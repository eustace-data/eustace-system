"""Parse command line and run gridding accordingly."""

__version__ = "$Revision: 1335 $"
__author__ = "Joel R. Mitchelson"

import cProfile
import pstats
import StringIO
import argparse
import json
import sys
from datetime import datetime
from file_aggregate_grid import SatelliteFileAggregator
from satellite_file import SatelliteFilename
from filebuilder import FileBuilderAggregateField


class ProgramResult(object):
    """Represent results from program run."""

    OUTPUT_STATUS_OK = 'ok'
    OUTPUT_STATUS_NO_DATA = 'no_data'

    def __init__(self, args, time_mode, load, output_status):
        self.args = args
        self.time_mode = time_mode
        self.load = load
        self.output_status = output_status


class FilenameList(object):
    """Represent a set of filenames built from input specification."""

    TIME_MODE_SINGLE = 'single'
    TIME_MODE_DAY = 'day'

    def __init__(self, args):

        # initialise to empty
        times = []

        if len(args.date) == 12:
            # requested a specific time
            times.append(datetime.strptime(args.date, '%Y%m%d%H%M'))
            time_mode = FilenameList.TIME_MODE_SINGLE
        else:
            # requested all times on a given day
            time_mode = FilenameList.TIME_MODE_DAY
            selected_date = datetime.strptime(args.date, '%Y%m%d').date()
            for hour in range(0, 24):
                for minute in range(0, 60, 5):
                    times.append(datetime(
                        selected_date.year, selected_date.month, selected_date.day, hour, minute))

        # store mode used
        self.time_mode = time_mode

        # build files from times
        self.filenames = [SatelliteFilename.from_pattern_and_time(
            args.path, args.pattern_lst, args.pattern_aux, t) for t in times]


def run(args):
    """Run aggregation from given args."""

    # build the list of times required based on input spec
    filename_list = FilenameList(args)

    # build aggregator
    aggregator = SatelliteFileAggregator.from_command_arguments(args)

    # run
    run_result = aggregator.run(filename_list.filenames)

    # store result if produced
    if run_result.fields:
        FileBuilderAggregateField.save(
            args.o, run_result.fields, args.source, args.complevel)
        output_status = ProgramResult.OUTPUT_STATUS_OK
    else:
        output_status = ProgramResult.OUTPUT_STATUS_NO_DATA

    # return results to be printed to stdout
    return ProgramResult(args, filename_list.time_mode, run_result.load, output_status)


def main():
    """Main program entry point."""

    # command line options
    parser = argparse.ArgumentParser('satgrid')
    parser.add_argument('-x0', default=-180.00, type=float, help='Longitude axis start')
    parser.add_argument('-xs', default=0.25, type=float, help='Longitude grid resolution')
    parser.add_argument('-xn', default=1440, type=int, help='Longitude number of grid points')
    parser.add_argument('-y0', default=-90.00, type=float, help='Latitude axis start')
    parser.add_argument('-ys', default=0.25, type=float, help='Latitude grid resolution')
    parser.add_argument('-yn', default=720, type=int, help='Latitude number of grid points')
    parser.add_argument('-qc_mask_obs', default=1, type=int, help='Bitmask to select the QC bits used to determine the satellite observations of interest (e.g. daytime)')
    parser.add_argument('-qc_filter_obs', default=0, type=int, help='Value taken by QC bits indicated by QC_MASK_OBS, at all satellite observations of interest')
    parser.add_argument('-qc_mask_valid', default=7, type=int, help='Bitmask to select the QC bits used to determine the subset of observations which are of interest and thought to be valid (e.g. daytime, high confidence, non-cloudy)')
    parser.add_argument('-qc_filter_valid', default=0, type=int, help='Value taken by QC bits indicated by QC_MASK_VALID on the valid subset.')
    parser.add_argument('-obsname', default='LST', help='Observation to aggregate')
    parser.add_argument('-o', default='satgrid.nc', help='Output filename')
    parser.add_argument('-profile', default=False, action='store_true', help='Run profiler')
    parser.add_argument('-pattern_lst', default='%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc')
    parser.add_argument('-pattern_aux', default='%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc')
    parser.add_argument('-path', default='/gws/nopw/j04/eustace/data/incoming/MODIS')
    parser.add_argument('-wrapcoords', default=True, action='store_true', help='Coordinate wrap-around')
    parser.add_argument('-maxbinobs', default=512, help='Upper limit on observations per grid box')
    parser.add_argument('-source', default='', help='Identifier to place in source attribute of output netcdf file.')
    parser.add_argument('-complevel', default=4, type=int, help='Compression level (0 means none)')
    parser.add_argument('date')

    # parse
    args = parser.parse_args()

    # start profiler if requested
    profile = None
    if args.profile:
        profile = cProfile.Profile()
        profile.enable()

    # run
    program_result = run(args)

    # print result information
    json.dump(program_result, sys.stdout, default=lambda o: o.__dict__)

    # print profiling options if requested
    if args.profile:

        # stop profiler
        profile.disable()

        # print profile info
        profile_results = StringIO.StringIO()
        profile_stats = pstats.Stats(
            profile, stream=profile_results).sort_stats('cumulative')
        profile_stats.print_stats()
        print profile_results.getvalue()
