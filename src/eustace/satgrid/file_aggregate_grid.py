"""Loop over satellite files and aggregate."""

__version__ = "$Revision: 472 $"
__author__ = "Joel R. Mitchelson"

import numpy
from grid import GridAxis
from grid import GridLatLon
from aggregate_grid import AggregateGridCountSumMinMax
from aggregate_grid import AggregateGridCount
from aggregate_grid import AggregateGridSumDevSq
from satellite_file import SatelliteFile

# This is a safe value to use to indicate invalid data:
# - Temperatures in K shouldn't go negative
# - Uncertainty or variance measures shouldn't go negative
FLAG_INVALID = -1.0


class LoadResult(object):
    """Represent result of loading a file."""

    STATUS_OK = 'ok'
    STATUS_MISSING = 'missing'
    STATUS_ERROR = 'error'

    def __init__(self, filename, status, message):
        self.filename = filename
        self.status = status
        self.message = message


class RunResult(object):
    """Represent results of running aggregation."""

    def __init__(self, fields, load):
        self.fields = fields
        self.load = load


class SatelliteFileAggregator(object):
    """Aggregate multiple satellite files onto one grid."""

    # Aggregator classes need many instance attributes to express all relevant information
    # pylint: disable=too-many-instance-attributes

    def __init__(self, axis_lat, axis_lon, observation_name, qc_mask_obs, qc_filter_obs, qc_mask_valid, qc_filter_valid, max_bin_obs):
        """Initialise with grid parameters and observation to aggregate."""

        # Allow lots of arguments here - could consider having data class for these in future
        # pylint: disable=too-many-arguments

        self.grid = GridLatLon(axis_lat, axis_lon)
        self.outputgrid = AggregateGridCountSumMinMax(self.grid.dimensions)
        self.obsgrid = AggregateGridCount(self.grid.dimensions)
        self.observation_name = observation_name
        self.qc_mask_obs = qc_mask_obs
        self.qc_filter_obs = qc_filter_obs
        self.qc_mask_valid = qc_mask_valid
        self.qc_filter_valid = qc_filter_valid

        # Grid for two-pass variance calculation
        self.devsqgrid = AggregateGridSumDevSq(self.grid.dimensions)

    def aggregate_file(self, satfile):
        """Aggregate the specified file."""

        # indices for valid data, variable at those points, and lat/lon coords
        # at those points
        valid_indices = satfile.get_indices_from_qc(self.qc_mask_valid, self.qc_filter_valid)
        observations = satfile.get_variable_at_indices(self.observation_name, valid_indices)
        uncertainty_random = satfile.get_variable_at_indices(self.observation_name + '_unc_ran', valid_indices)
        uncertainty_local_atm = satfile.get_variable_at_indices(self.observation_name + '_unc_loc_atm', valid_indices)
        uncertainty_local_sfc = satfile.get_variable_at_indices(self.observation_name + '_unc_loc_sfc', valid_indices)
        uncertainty_systematic = satfile.get_variable(self.observation_name + '_unc_sys')
        coord_lat = satfile.get_variable_at_indices('lat', valid_indices)
        coord_lon = satfile.get_variable_at_indices('lon', valid_indices)

        # ingest observations
        self.outputgrid.aggregate_at_indices(self.grid.compute_indices(coord_lat, coord_lon), observations, uncertainty_random, uncertainty_local_atm, uncertainty_local_sfc, uncertainty_systematic)

        # valid indices at all obs of interest (e.g. daytime), and
        # corresponding lat/lon
        obs_indices = satfile.get_indices_from_qc(self.qc_mask_obs, self.qc_filter_obs)
        obs_coord_lat = satfile.get_variable_at_indices('lat', obs_indices)
        obs_coord_lon = satfile.get_variable_at_indices('lon', obs_indices)

        # count obs of interest
        self.obsgrid.aggregate_at_indices(self.grid.compute_indices(obs_coord_lat, obs_coord_lon))

    def aggregate_file_second_order(self, fields0, satfile):
        """Aggregate for second-order fields, having already run once."""

        valid_indices = satfile.get_indices_from_qc(self.qc_mask_valid, self.qc_filter_valid)
        obs = satfile.get_variable_at_indices(self.observation_name, valid_indices)
        grid_bin_obs_mean = fields0['tsmean']
        coord_lat = satfile.get_variable_at_indices('lat', valid_indices)
        coord_lon = satfile.get_variable_at_indices('lon', valid_indices)
        bin_indices = self.grid.compute_indices(coord_lat, coord_lon)
        self.devsqgrid.aggregate_at_indices(grid_bin_obs_mean, bin_indices, obs)

    def compute_fields_first_order(self):
        """Compute the fields which can be computed without knowing the mean."""

        fields = {}
        fields['latitude'] = self.grid.axis_lat.get_centres()
        fields['longitude'] = self.grid.axis_lon.get_centres()
        fields.update(self.outputgrid.compute_fields(flag_invalid=FLAG_INVALID))
        fields['total_number_of_observations'] = self.obsgrid.count
        return fields

    def compute_fields_second_order(self):
        """Compute fields which require aggregate_file_second_order to have been run."""

        valid_indices = numpy.nonzero(self.outputgrid.count > numpy.int32(1))
        sum_devsq_valid = self.devsqgrid.sum_devsq[valid_indices]
        nvalid = self.outputgrid.count[valid_indices]
        ntotal = self.obsgrid.count[valid_indices]

        # compute variance over valid observations
        variance_valid = sum_devsq_valid / (nvalid - 1)
        sampling_uncertainty_valid = (
            (ntotal - nvalid) * numpy.sqrt(variance_valid)) / (ntotal - 1)

        # make arrays all flagged as invalid
        variance = numpy.empty(shape=self.outputgrid.dimensions, dtype=numpy.float32)
        variance.fill(FLAG_INVALID)
        sampling_uncertainty = numpy.empty(shape=self.outputgrid.dimensions, dtype=numpy.float32)
        sampling_uncertainty.fill(FLAG_INVALID)

        # populate with results
        variance[valid_indices] = variance_valid
        sampling_uncertainty[valid_indices] = sampling_uncertainty_valid

        # return as fields
        return {'tsvariance': variance, 'tsmean_unc_spl': sampling_uncertainty}

    def run(self, filenames):
        """Perform aggregation on files, specified by file name."""

        # load of attempts to load files and status
        loadresults = []

        # load all files
        satfiles = []
        for filename in filenames:

            # initialise status to missing
            status = LoadResult.STATUS_MISSING
            message = None

            # load if exists, otherwise log non-existence and continue
            if filename.check_readable():

                # Catching all exceptions is appropriate here as the intention is
                # to output any exception method to the log, no matter what it was.
                # pylint: disable=broad-except

                try:

                    # load data
                    satfiles.append(SatelliteFile(filename))

                    # flag as ok
                    status = LoadResult.STATUS_OK

                except Exception as detail:

                    # flag as error
                    status = LoadResult.STATUS_ERROR

                    # record error message
                    message = detail.__str__()

            # log results of this file
            loadresults.append(LoadResult(filename, status, message))

        # process if there are any files present
        if satfiles:

            # first pass
            for satfile in satfiles:
                self.aggregate_file(satfile)

            # compute first order fields
            fields = self.compute_fields_first_order()

            # second pass for squared deviation
            for satfile in satfiles:
                self.aggregate_file_second_order(fields, satfile)

            # second order fields
            fields.update(self.compute_fields_second_order())

        else:

            # no files present - empty fields item
            fields = None

        # build result
        return RunResult(fields, loadresults)

    @staticmethod
    def from_command_arguments(args):
        """ Construct aggregation instance using arguments from specified args namespace."""

        return SatelliteFileAggregator(
            axis_lat=GridAxis(args.y0, args.ys, args.yn, args.wrapcoords),
            axis_lon=GridAxis(args.x0, args.xs, args.xn, args.wrapcoords),
            qc_mask_obs=args.qc_mask_obs, qc_filter_obs=args.qc_filter_obs,
            qc_mask_valid=args.qc_mask_valid, qc_filter_valid=args.qc_filter_valid,
            observation_name=args.obsname,
            max_bin_obs=args.maxbinobs)
