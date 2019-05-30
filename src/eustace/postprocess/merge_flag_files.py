"""

Merge the calendar and threshold mask files

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

from operation_count import operation_dates
from apply_mask import load_flags, save_flag_file

import flags

def merged_flags(input_directory_1, input_directory_2, output_directory, iteration, processdate):
    """merge two flag files"""
    
    datestring = "{:04d}{:02d}{:02d}".format(processdate.year, processdate.month, processdate.day )

    flags_1 = load_flags(input_directory_1, iteration, datestring)
    flags_1 = numpy.squeeze(flags_1.data[:,:,:])
    
    flags_2 = load_flags(input_directory_2, iteration, datestring)
    flags_2 = numpy.squeeze(flags_2.data[:,:,:])
    
    # work around for incorrect flag usage
    incorrect_flag = numpy.bitwise_and(flags_2, flags.AREAL_LOW_FLAG) == flags.AREAL_LOW_FLAG
    #print incorrect_flag
    flags_2[incorrect_flag] = flags_2[incorrect_flag] | flags.MISSING_MARINE_FLAG # set the correct flag
    flags_2 = flags_2 & ~flags.AREAL_LOW_FLAG # negate the erroneous flag
    
    merged_flags = flags_1 | flags_2
    
    return merged_flags


"""

Calls setup to run as batches in lsf

"""

def merging_operation( reference_time_string,
                        operation_index,
                        iteration,
                        input_directory_1,
                        input_directory_2,
                        output_directory,
                        ):

    # get dates in the month to be processed in this operation
    reference_time = parser.parse(reference_time_string)
    processing_dates = operation_dates(reference_time, operation_index)
    
    # derive flags for each day
    for processdate in processing_dates:
        
        # check that directory exists and if not then make it
        if not os.path.exists(os.path.join(output_directory, str(processdate.year))):
            os.makedirs(os.path.join(output_directory, str(processdate.year)))

        flag_values = merged_flags(input_directory_1, input_directory_2, output_directory, iteration, processdate)
        
        # save to NetCDF
        outputfile = os.path.join(output_directory, '{:04d}'.format(processdate.year), 'eustace_analysis_{:d}_qc_flags_{:04d}{:02d}{:02d}.nc'.format(iteration, processdate.year, processdate.month, processdate.day))
        save_flag_file(flag_values, processdate, outputfile)
 
    
def main():

    print 'Submission of flag merging jobs'
    
    parser = argparse.ArgumentParser(description='Monthly batches of masking operations')
    
    parser.add_argument('--reference_time_string', type=str, default="1880-01-01", help='reference first day to run formatted YYY-MM-DD. Should be first day of month.')
    parser.add_argument('--operation_index',  type=int, default=0, help='the number of months since the reference_time_string that this masking operation corresponds to')
    parser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    parser.add_argument('--input_directory_1', type=str, default="/work/scratch/cmorice/masking/", help='first directory containing flag files to be merged')
    parser.add_argument('--input_directory_2',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='second directory containing flag files to be merged')
    parser.add_argument('--output_directory', type=str, default="/work/scratch/cmorice/masking/", help='directory in which threshold masking for largescale, climatology and extremes can be found')
    args = parser.parse_args()

    merging_operation(  args.reference_time_string,
                        args.operation_index,
                        args.iteration,
                        args.input_directory_1,
                        args.input_directory_2,
                        args.output_directory
                         )

if __name__ == '__main__':
    
    main()
    

