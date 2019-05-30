
import os.path
import iris
from datetime import datetime
import numpy
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystemInputLoader
from eustace.analysis.advanced_standard.fileio.observation_structure_source_connector import ObservationStructureSourceConnector

iris.FUTURE.netcdf_promote = True
"""
Required to avoid this iris warning:

    IrisDeprecation: NetCDF default loading behaviour currently does not expose variables which
    define reference surfaces for dimensionless vertical coordinates as independent Cubes.
    This behaviour is deprecated in favour of automatic promotion to Cubes.
    To switch to the new behaviour, set iris.FUTURE.netcdf_promote to True.
"""

class AnalysisSystemInputLoaderHadCRUT4(AnalysisSystemInputLoader):
    """Load HadCRUT4 data at specified time indices on demand."""

    def __init__(self, filelist):
        """Construct for given base path and time index model."""

        # Keep track of file list
        self.filelist = filelist
        
        # Reference epoch
        timeunits = iris.load(self.filelist[0])[0].coord('time').units
        self.epoch = timeunits.num2date(0)

    def datetime_at_time_index(self, time_index):
        """Convert the given time index to a datetime.
           For HadCRUT4 this assumes time_index is a month number."""

        return ObservationSourceHadCRUT4.midmonth(self.epoch, time_index)

    def load_observation_source(self, time_index):
        """Load HadCRUT4 data at specified time index, as ObservationSource instance."""

        # Load as ObservationSource instance
        return ObservationSourceHadCRUT4(self.filelist, time_index)

    def load_observation_structure(self, observable, time_index, log=None):
        """Load HadCRUT4 data at specified time index, as ObservationStructure instance."""

        # Compute datetime for this time index
        current_date = self.datetime_at_time_index(time_index)

        # Get source
        source = self.load_observation_source(time_index)

        # Convert to observation structure instance
        return ObservationStructureSourceConnector(
            observationsource=source,
            observable=observable,
            corresponding_datetime=current_date)

class ObservationSourceHadCRUT4(ObservationSource):
    """ObservationSource interface to HadCRUT4."""

    @staticmethod
    def midmonth(epoch, monthnumber):
        """Compute the mid-point of given month number since epoch and return as datetime object.
           This is needed because HadCRUT4 uses midmonth values on the time axis."""

        # Convert month number to the mid-month point
        total = 12*epoch.year + monthnumber
        t0 = datetime(total / 12, 1 + (total % 12), 1)
        t1 = datetime((total + 1)/12, 1 + ((total + 1) % 12), 1)
        halfmonth = (t1 - t0) / 2
        return t0 + halfmonth

    @staticmethod
    def monthnumber_to_decimal(referencefile, monthnumber):

        # Reference epoch
        timeunits = iris.load(referencefile)[0].coord('time').units
        epoch = timeunits.num2date(0)

        # Midmonth point
        processdate = ObservationSourceHadCRUT4.midmonth(epoch, monthnumber)

        # Compute corresponding decimal time
        return timeunits.date2num(processdate)

    def __init__(self, filelist, monthnumber):
        """Initialise using data on given base path and specified integer month number."""

        # Store filelist
        self.filelist = filelist

        # Store month number
        self.monthnumber = monthnumber

        # Convert to decimal request
        requesttime = ObservationSourceHadCRUT4.monthnumber_to_decimal(filelist[0], monthnumber)

        # Constraints for Iris to select required data
        self.select_time = iris.Constraint(time=lambda cell: cell.point == requesttime)
        self.select_tas = iris.Constraint(cube_func=lambda cube:cube.var_name=='temperature_anomaly')
        self.select_unc_tas = iris.Constraint(cube_func=lambda cube:cube.var_name=='standard_error')

    def observables(self):
        """The names of variables estimated from this source."""

        return [ ObservationSource.TMEAN ]

    def observation_location_lookup(self):
        """NumPy array in which column number corresponds to location id and rows are latitude and longitude."""

        tascube = iris.load(self.filelist)[1]
        lon, lat = numpy.meshgrid(numpy.float64(tascube.coord('longitude').points), numpy.float64(tascube.coord('latitude').points))
        return numpy.vstack((lat.ravel(), lon.ravel()))

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""

        if observable == ObservationSource.TMEAN:

            # Load data
            tas_cube = iris.load(self.filelist, self.select_tas & self.select_time)
            unc_tas_cube = iris.load(self.filelist, self.select_unc_tas & self.select_time)

            # Get numerical arrays
            tas = tas_cube.concatenate_cube().data.ravel()
            tas_unc0 = unc_tas_cube[0].data.ravel()
            tas_unc1 = unc_tas_cube[1].data.ravel()

            # Valid indices for temperature obs
            valid_indices = numpy.nonzero(~tas.mask)

            # Check uncertainty information covers inputs
            if any(tas_unc0.mask[valid_indices]) or any(tas_unc1.mask[valid_indices]):
                raise ValueError('missing uncertainty information')

            # Ideally would compute this just by calling:
            #
            #     unc_tas_combined = numpy.sqrt(tas_unc0**2 + tas_unc1**2)
            #
            # but that raises a warning on the invalid values, so must specify the valid ones only.
            #
            unc_tas_combined = numpy.zeros(shape=tas.shape)
            unc_tas_combined[valid_indices] = numpy.sqrt(tas_unc0.data[valid_indices]**2 + tas_unc1.data[valid_indices]**2)
            
            # Check all unc_tas_combined are reasonable
            missing_uncertainty_indices = numpy.nonzero(~tas.mask & (unc_tas_combined < 0.001))[0]
            if len(missing_uncertainty_indices):
                message = 'missing uncertainties at {0}: {1}'.format(missing_uncertainty_indices, unc_tas_combined[missing_uncertainty_indices])
                print 'WARNING: ', message
                print '--> adjusting to sqrt(0.1) K'
                unc_tas_combined[missing_uncertainty_indices] = numpy.sqrt(0.1)
                
            # Return as observations class
            return Observations(
                mask=tas.mask.ravel(),
                time=self.monthnumber,
                location=numpy.array(range(tas.shape[0]), dtype=numpy.int64),
                measurement=tas.data,
                uncorrelatederror=unc_tas_combined,
                locallycorrelatederror=None)
        
    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty."""
        
        return [ ]
