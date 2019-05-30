"""Example catalogue use: define satellite gridding over land."""

import argparse
from eumopps.catalogue.operationcommand import OperationArguments
from eumopps.timeutils import datetime_numeric
from eustace.timeutils.epoch import days_since_epoch
from eustace.timeutils.epoch import epoch_plus_days
from eustace.satgrid.satellite_file import SatelliteFilename
from eustace.satgrid.file_aggregate_grid import SatelliteFileAggregator
from eustace.satgrid.filebuilder import FileBuilderAggregateField

INPUTNAME = 'Aqua-MODIS'


class LandSatGridArguments(OperationArguments):
    """Define configuration options for day and night subsets."""

    # Allow multiple two-character names as on command line
    # pylint: disable=invalid-name,too-many-instance-attributes,too-many-locals

    def __init__(self, prefix, qc_filter_obs, qc_filter_valid):
        """Build based on a few configuration parameters."""

        super(LandSatGridArguments, self).__init__()

        # the configurable bits which are different between day and night
        self.outputpattern = prefix + '.%Y%m%d.nc'
        self.qc_filter_obs = qc_filter_obs
        self.qc_filter_valid = qc_filter_valid

        # fixed parameters
        self.qc_mask_obs = 1
        self.qc_mask_valid = 7
        self.x0 = -180.0
        self.xs = 0.25
        self.xn = 1440
        self.y0 = -90.0
        self.ys = 0.25
        self.yn = 720
        self.obsname = 'LST'
        self.wrapcoords = True
        self.maxbinobs = 512
        self.complevel = 4
        self.source = ''

SUBSETSPECS = [LandSatGridArguments('satgrid.day', 0, 0), LandSatGridArguments('satgrid.night', 1, 1)]
"""Subsets to make."""


def specify(manager, options):
    """Return the list of operations and expected output given the input catalogue."""

    # Require many local variables here
    # pylint: disable=too-many-locals

    # parse input option list
    parser = argparse.ArgumentParser()
    parser.add_argument('--start', required=True)
    parser.add_argument('--end', required=True)
    args = parser.parse_args(options)

    # process a whole number of days
    startday = int(days_since_epoch(datetime_numeric.parse(args.start)))
    endday = int(days_since_epoch(datetime_numeric.parse(args.end)))

    # find items in catalogue grouped by day
    lstinput = manager.references_groupbyday(INPUTNAME, subsetindex=0)
    auxinput = manager.references_groupbyday(INPUTNAME, subsetindex=1)

    # the main output dataset
    dataset = manager.newdataset()

    # list of subsets (one for each command pattern)
    subsets = [dataset.newsubset([spec.outputpattern]) for spec in SUBSETSPECS]

    # iterate over subsets
    for subsetindex, subset in enumerate(subsets):

        # process each day
        for dayindex in range(startday, endday + 1):

            # convert to datetime object for formatting etc
            day = epoch_plus_days(dayindex)

            # loop over LST files for this day
            inputs = []
            for lstreference in lstinput[dayindex]:

                # find AUX file matching the LST file
                lsttime = manager.match(lstreference).time
                print 'day {0} time {1}'.format(dayindex, lsttime)
                auxreference = next(reference for reference in auxinput[dayindex] if manager.match(reference).time == lsttime)

                # append pair
                inputs.extend([lstreference, auxreference])

            # build output filename
            outputs = [subset.newfiletime(day)]

            # Append this operation to make it
            dataset.newoperation(inputs, outputs, SUBSETSPECS[subsetindex])


def run(inputpathnames, outputpathnames, args):
    """Run aggregation on specified input pathnames (expected LST/AUX file pairs) and output according to args."""

    # convert flat list of pathnames into satellite file objects
    inputs = []
    for inputindex in range(0, len(inputpathnames), 2):
        inputs.append(SatelliteFilename(
            inputpathnames[inputindex], inputpathnames[inputindex + 1]))

    # the aggregation class
    aggregator = SatelliteFileAggregator.from_command_arguments(args)

    # run it
    run_result = aggregator.run(inputs)

    # save result
    print 'Saving: ', outputpathnames[0]
    FileBuilderAggregateField.save(
        outputpathnames[0], run_result.fields, args.source, args.complevel)
