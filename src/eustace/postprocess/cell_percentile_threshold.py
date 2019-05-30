"""Value QC check on EUSTACE analysis products

Flags computed:

1. Climatology uncertainty threshold on this day
2. Large scale uncertainty threshold on this day
3. Local constraint threshold on this day

"""

import argparse
import os.path
import numpy
from dateutil.relativedelta import relativedelta
import datetime
from dateutil import parser

import iris

from eustace.outputformats.outputvariable import OutputVariable
from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eumopps.version.svn import get_revision_id_for_module
from netCDF4 import default_fillvals
from eustace.outputformats import definitions
import eustace.timeutils.epoch

from apply_mask import load_analysis, load_climatology, load_large_scale, load_local, save_flag_file

import outlier_check
from outlier_check import compute_window_mean
from dateutil import parser
from operation_count import operation_dates
import flags

# bit mask flags used in output file
FLAG_TYPE = flags.FLAG_TYPE
TYPE_NAME = flags.TYPE_NAME
FLAG_MAX_N = flags.FLAG_MAX_N
FLAG_MAX_USED = flags.FLAG_MAX_USED
NULL_FLAG = flags.NULL_FLAG

# time window checks
DAY_FLAG = flags.DAY_FLAG
PRIOR_WINDOW_FLAG = flags.PRIOR_WINDOW_FLAG
POST_WINDOW_FLAG = flags.POST_WINDOW_FLAG

# calendar day checks
CALENDAR_DAY_FLAG    = flags.CALENDAR_DAY_FLAG
PRIOR_CALENDAR_FLAG  = flags.PRIOR_CALENDAR_FLAG
POST_CALENDAR_FLAG   = flags.POST_CALENDAR_FLAG

# slow component uncertatinty thresholds
CLIMATOLOGY_UNC_FLAG  = flags.CLIMATOLOGY_UNC_FLAG
LARGE_SCALE_UNC_FLAG  = flags.LARGE_SCALE_UNC_FLAG

# extreme value checks
AREAL_LOW_FLAG      = flags.AREAL_LOW_FLAG
AREAL_HIGH_FLAG     = flags.AREAL_HIGH_FLAG
EXTREME_LOW_FLAG    = flags.EXTREME_LOW_FLAG
EXTREME_HIGH_FLAG   = flags.EXTREME_HIGH_FLAG

# omitted data source flags
MISSING_MARINE_FLAG  = flags.MISSING_MARINE_FLAG

#
LOCATION_LOW_FLAG = flags.LOCATION_LOW_FLAG
LOCATION_HIGH_FLAG = flags.LOCATION_HIGH_FLAG

# missing data indicator
MISSING_FLAG_FLAG    = flags.MISSING_FLAG_FLAG

FLAG_MEANINGS = flags.FLAG_MEANINGS

# parameters for uncertainty thresholding
QUARTILE_LIMIT = 8.0


"""

Core methods for thresholding

"""

def location_threshold_checks(analysis_directory, percentile_directory, output_directory, iteration, processdate):
    """check constraint in a local window of +/- n days"""
    
    # date time filename formated datestring for day on interest
    datestring = "{:04d}{:02d}{:02d}".format(processdate.year, processdate.month, processdate.day )
    print processdate
    
    # get the component field for this date
    try:
        analysis, uncertainty, influence, ensemble = load_analysis(analysis_directory, iteration, datestring)
    except:
        print "Error loading large scale file" 

    analysis_field = numpy.squeeze( analysis.data )

    percentile_file = os.path.join(percentile_directory, 'eustace_analysis_{:d}_percentiles_{:02d}{:02d}.npz'.format(iteration, processdate.month, processdate.day))
    percentiles = numpy.load( percentile_file )
    
    lower_quartile = percentiles['lower_quartile_block']
    median = percentiles['median_block']
    upper_quartile = percentiles['upper_quartile_block']

    # check analysis values exceed thresholds
    lower_threshold_exceedance = analysis_field < median - QUARTILE_LIMIT * (median-lower_quartile)
    upper_threshold_exceedance = analysis_field > median + QUARTILE_LIMIT * (upper_quartile-median)

    flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )    
    flag_values[lower_threshold_exceedance] = flag_values[lower_threshold_exceedance] | LOCATION_LOW_FLAG
    flag_values[upper_threshold_exceedance] = flag_values[upper_threshold_exceedance] | LOCATION_HIGH_FLAG
    
    return flag_values

"""

Calls setup to run as batches in lsf

"""

def flagging_operation(  reference_time_string,
                            operation_index,
                            analysis_directory,
                            iteration,
                            output_directory,
                            percentile_directory,
                            ):
    
    # get dates in the month to be processed in this operation
    reference_time = parser.parse(reference_time_string)
    processing_dates = operation_dates(reference_time, operation_index)
    
    # derive flags for each day
    for processdate in processing_dates:
        
        # check that directory exists and if not then make it
        if not os.path.exists(os.path.join(output_directory, str(processdate.year))):
            os.makedirs(os.path.join(output_directory, str(processdate.year)))
        
        # compute flags
        flag_values = location_threshold_checks(analysis_directory, percentile_directory, output_directory, iteration, processdate)
        
        #save
        outputfile = os.path.join(output_directory, '{:04d}'.format(processdate.year), 'eustace_analysis_{:d}_qc_flags_{:04d}{:02d}{:02d}.nc'.format(iteration, processdate.year, processdate.month, processdate.day))
        save_flag_file(flag_values, processdate, outputfile)

def main():

    print 'Submission of advanced standard analysis location based threshold jobs'
    
    parser = argparse.ArgumentParser(description='Monthly batches of masking operations')
    
    parser.add_argument('--reference_time_string', type=str, default="1880-01-01", help='reference first day to run formatted YYY-MM-DD. Should be first day of month.')
    parser.add_argument('--operation_index',  type=int, default=0, help='the number of months since the reference_time_string that this masking operation corresponds to')
    parser.add_argument('--analysis_directory', type=str, default="/work/scratch/cmorice/advanced_standard/", help='root directory of the analysis')
    parser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    parser.add_argument('--output_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='output location')
    parser.add_argument('--percentile_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='percentile_directory')

    args = parser.parse_args()

    flagging_operation( args.reference_time_string,
                        args.operation_index,
                        args.analysis_directory,
                        args.iteration,
                        args.output_directory,
                        args.percentile_directory )

if __name__ == '__main__':
    
    main()
    
