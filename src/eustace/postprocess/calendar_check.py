"""Calendar day checks in climatology and large scale threshold flags for EUSTACE analysis products"""

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

from apply_mask import load_flags, save_flag_file

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

# missing data indicator
MISSING_FLAG_FLAG    = flags.MISSING_FLAG_FLAG

FLAG_MEANINGS = flags.FLAG_MEANINGS

# useful derived combinations
CALENDAR_DAY_FLAG    = flags.CALENDAR_DAY_FLAG # helper combination that say there is no observation nearby in time

# will be flagged when the number of constrained days is less than COUNT_THRESHOLD 
COUNT_THRESHOLD = 1

# will be flagged when the number of constrained days in whole period is less than TOTAL_COUNT_THRESHOLD 
TOTAL_COUNT_THRESHOLD = 3

"""

File I/O

"""             

def calendar_constraint_checks(input_directory, output_directory, iteration, reference_time, start_year, end_year, target_flag):
    """check constraint in a local window of +/- n days"""
    
    # initialise contraint flags to None
    prior_constraint_count = None
    post_constraint_count = None
    total_constraint_count = None
    # loop over window to check the constraint
    print "looping over:", start_year-reference_time.year, end_year-reference_time.year+1
    for yearshift in range(start_year-reference_time.year, end_year-reference_time.year+1):
    
        # get the day's analysis data
        day_to_load = reference_time + relativedelta(years=yearshift)
        datestring = "{:04d}{:02d}{:02d}".format(day_to_load.year, day_to_load.month, day_to_load.day )
        try:
            reference_flags = load_flags(input_directory, iteration, datestring)
            reference_flags = numpy.squeeze(reference_flags.data[:,:,:])
        except:
            print "failed to load", datestring
            print input_directory
            continue
    
        isconstrained = ~(numpy.bitwise_and(reference_flags, target_flag) == target_flag)
        
        if total_constraint_count is None:
            total_constraint_count = numpy.zeros(isconstrained.shape, numpy.int64)
        total_constraint_count[isconstrained]+=1
        
        if day_to_load.year == reference_time.year:
            # get and keep the flag values for the centre day
            flag_values = reference_flags
            
        else:
            # check to see if constrained and accumulate prior/post count

            if prior_constraint_count is None:
                prior_constraint_count = numpy.zeros(isconstrained.shape, numpy.int64)
            if post_constraint_count is None:
                post_constraint_count = numpy.zeros(isconstrained.shape, numpy.int64)
            
            # accumulate number of previous days constrained in window
            if day_to_load.year < reference_time.year:
                prior_constraint_count[isconstrained]+=1

            # accumulate number of following days constrained in window
            if day_to_load.year > reference_time.year:
                post_constraint_count[isconstrained]+=1

    # compute summary flags          
    prior_constrained = prior_constraint_count < COUNT_THRESHOLD
    post_constrained = post_constraint_count < COUNT_THRESHOLD

    climate_constrained = total_constraint_count < TOTAL_COUNT_THRESHOLD

    # set flag values for output
    flag_values[prior_constrained] = flag_values[prior_constrained] | PRIOR_CALENDAR_FLAG
    flag_values[post_constrained] = flag_values[post_constrained] | POST_CALENDAR_FLAG
    flag_values[climate_constrained] = flag_values[climate_constrained] | MISSING_FLAG_FLAG
    
    return flag_values

"""

Calls setup to run as batches in lsf

"""

def flagging_operation( reference_time_string,
                        operation_index,
                        input_directory,
                        iteration,
                        output_directory,
                        start_year,
                        end_year ):

    # get dates in the month to be processed in this operation
    reference_time = parser.parse(reference_time_string)
    processing_dates = operation_dates(reference_time, operation_index)
    
    # derive flags for each day
    for processdate in processing_dates:
        
        # check that directory exists and if not then make it
        if not os.path.exists(os.path.join(output_directory, str(processdate.year))):
            os.makedirs(os.path.join(output_directory, str(processdate.year)))

        flag_values = calendar_constraint_checks( input_directory,
                                                  output_directory,
                                                  iteration,
                                                  processdate,
                                                  start_year,
                                                  end_year,
                                                  CALENDAR_DAY_FLAG )

        # save to NetCDF
        outputfile = os.path.join(output_directory, '{:04d}'.format(processdate.year), 'eustace_analysis_{:d}_qc_flags_{:04d}{:02d}{:02d}.nc'.format(iteration, processdate.year, processdate.month, processdate.day))
        save_flag_file(flag_values, processdate, outputfile)


def main():

    print 'Submission of advanced standard analysis jobs'
    
    parser = argparse.ArgumentParser(description='Monthly batches of masking operations')
    
    parser.add_argument('--reference_time_string', type=str, default="1880-01-01", help='reference first day to run formatted YYY-MM-DD. Should be first day of month.')
    parser.add_argument('--operation_index',  type=int, default=0, help='the number of months since the reference_time_string that this masking operation corresponds to')
    parser.add_argument('--input_directory', type=str, default="/work/scratch/cmorice/advanced_standard/", help='directory in which masking information can be found')
    parser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    parser.add_argument('--output_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='operation index at which the analysis grid produced')
    parser.add_argument('--start_year',  type=int, default=1880, help='first analysis year to load for constraint tests')
    parser.add_argument('--end_year',  type=int, default=2015, help='last analysis year to load for constraint tests')
    
    args = parser.parse_args()

    flagging_operation( args.reference_time_string,
                        args.operation_index,
                        args.input_directory,
                        args.iteration,
                        args.output_directory,
                        args.start_year,
                        args.end_year
                       )

if __name__ == '__main__':
    
    main()
    
