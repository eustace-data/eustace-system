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

# missing data indicator
MISSING_FLAG_FLAG    = flags.MISSING_FLAG_FLAG

FLAG_MEANINGS = flags.FLAG_MEANINGS

# parameters for uncertainty thresholding
CLIMATOLOGY_UNC_SMOOTH_WINDOW = 13  # window for smoothing climatology uncertainties before thresholding in number of grid cells
LARGE_SCALE_UNC_SMOOTH_WINDOW = 17  # window for smoothing large-scale uncertainties before thresholding in number of grid cells
LOCAL_INFLUENCE_SMOOTH_WINDOW = 5  # window for smoothing local constraints before thresholding in number of grid cells
PADMODE_LATITUDE='reflect'
PADMODE_LONGITUDE='wrap'

"""

Land sea mask fix

"""

def get_marine_flag():

    LAND_VALUE = 0.0
    SEA_VALUE  = 100.0
    COAST_FILE='/gws/nopw/j04/eustace/data/internal/climatology_covariates/coastal_influence.test.0.25_0.25.nc'
    COAST_FILE_VAR = 'coastal_influence'
    variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == COAST_FILE_VAR))
    coastal_influence = iris.load_cube( COAST_FILE, variable_constraint ).data[:,:]

    coastal_influence[0,:] = LAND_VALUE
    coastal_influence[-1,:] = SEA_VALUE
    coastal_influence[:,0] = coastal_influence[:,1] 
    coastal_influence[:,-1] = coastal_influence[:,-2]
    
    marine_locations = coastal_influence == 100.0
    
    flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )  
    flag_values[marine_locations] = flag_values[marine_locations] | MISSING_MARINE_FLAG
    
    return flag_values
    

"""

Core methods for thresholding on climatology and large scale uncertainty influence

"""

def climatology_uncertainty_checks(analysis_directory, output_directory, iteration, processdate, uncertainty_threshold):
    """check constraint in a local window of +/- n days"""
    
    # date time filename formated datestring for day on interest
    datestring = "{:04d}{:02d}{:02d}".format(processdate.year, processdate.month, processdate.day )
    print processdate
    
    # get the component field for this date
    try:
        analysis, uncertainty, influence = load_climatology(analysis_directory, iteration, datestring)
    except:
        print "Error loading climatology file" 

    climatology_uncertainty_field = numpy.squeeze( uncertainty.data )
    
    # check whether smoothed uncertainties are greater than the uncertainty_threshold
    climatology_uncertainty_exceedance = compute_window_mean(climatology_uncertainty_field,
                                                             CLIMATOLOGY_UNC_SMOOTH_WINDOW, 
                                                             CLIMATOLOGY_UNC_SMOOTH_WINDOW, 
                                                             padmode_latitude=PADMODE_LATITUDE, 
                                                             padmode_longitude=PADMODE_LONGITUDE) > uncertainty_threshold

    flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )    
    flag_values[climatology_uncertainty_exceedance] = flag_values[climatology_uncertainty_exceedance] | CLIMATOLOGY_UNC_FLAG
    
    return flag_values

def large_scale_uncertainty_checks(analysis_directory, output_directory, iteration, processdate, uncertainty_threshold):
    """check constraint in a local window of +/- n days"""
    
    # date time filename formated datestring for day on interest
    datestring = "{:04d}{:02d}{:02d}".format(processdate.year, processdate.month, processdate.day )
    print processdate
    
    # get the component field for this date
    try:
        analysis, uncertainty, influence = load_large_scale(analysis_directory, iteration, datestring)
    except:
        print "Error loading large scale file" 

    large_scale_uncertainty_field = numpy.squeeze( uncertainty.data )
    
    # check whether smoothed uncertainties are greater than the uncertainty_threshold
    large_scale_uncertainty_exceedance = compute_window_mean(large_scale_uncertainty_field,
                                                             LARGE_SCALE_UNC_SMOOTH_WINDOW, 
                                                             LARGE_SCALE_UNC_SMOOTH_WINDOW, 
                                                             padmode_latitude=PADMODE_LATITUDE, 
                                                             padmode_longitude=PADMODE_LONGITUDE) > uncertainty_threshold

    flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )    
    flag_values[large_scale_uncertainty_exceedance] = flag_values[large_scale_uncertainty_exceedance] | LARGE_SCALE_UNC_FLAG
    
    return flag_values

def local_influence_checks(analysis_directory, output_directory, iteration, processdate, constraint_threshold):
    """check constraint in a local window of +/- n days"""
    
    # date time filename formated datestring for day on interest
    datestring = "{:04d}{:02d}{:02d}".format(processdate.year, processdate.month, processdate.day )
    print processdate
    
    # get the component field for this date
    try:
        analysis, uncertainty, influence, ensemble = load_analysis(analysis_directory, iteration, datestring)
    except:
        print "Error loading large scale file" 

    local_influence_field = numpy.squeeze( influence.data )
    
    # check whether smoothed uncertainties are greater than the uncertainty_threshold
    local_influence_exceedance = compute_window_mean(local_influence_field,
                                                             LOCAL_INFLUENCE_SMOOTH_WINDOW, 
                                                             LOCAL_INFLUENCE_SMOOTH_WINDOW, 
                                                             padmode_latitude=PADMODE_LATITUDE, 
                                                             padmode_longitude=PADMODE_LONGITUDE) < constraint_threshold

    flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )    
    flag_values[local_influence_exceedance] = flag_values[local_influence_exceedance] | DAY_FLAG
    
    return flag_values

"""

Calls setup to run as batches in lsf

"""

def flagging_operation(  reference_time_string,
                            operation_index,
                            analysis_directory,
                            iteration,
                            output_directory,
                            climatology_limit,
                            largescale_limit,
                            constraint_threshold ):
    
    # get dates in the month to be processed in this operation
    reference_time = parser.parse(reference_time_string)
    processing_dates = operation_dates(reference_time, operation_index)
    
    # derive flags for each day
    for processdate in processing_dates:
        
        # check that directory exists and if not then make it
        if not os.path.exists(os.path.join(output_directory, str(processdate.year))):
            os.makedirs(os.path.join(output_directory, str(processdate.year)))
        
        # compute flags
        climatology_flag_values = climatology_uncertainty_checks(analysis_directory, output_directory, iteration, processdate, climatology_limit)
        large_scale_flag_values = large_scale_uncertainty_checks(analysis_directory, output_directory, iteration, processdate, largescale_limit)
        local_influence_flag_values = local_influence_checks(analysis_directory, output_directory, iteration, processdate, constraint_threshold)
        marine_flag = get_marine_flag()
        
        # join flags
        flag_values = climatology_flag_values | large_scale_flag_values | local_influence_flag_values | marine_flag
        
        #save
        outputfile = os.path.join(output_directory, '{:04d}'.format(processdate.year), 'eustace_analysis_{:d}_qc_flags_{:04d}{:02d}{:02d}.nc'.format(iteration, processdate.year, processdate.month, processdate.day))
        save_flag_file(flag_values, processdate, outputfile)

def main():

    print 'Submission of advanced standard analysis observation constraint flagging jobs'
    
    parser = argparse.ArgumentParser(description='Monthly batches of masking operations')
    
    parser.add_argument('--reference_time_string', type=str, default="1880-01-01", help='reference first day to run formatted YYY-MM-DD. Should be first day of month.')
    parser.add_argument('--operation_index',  type=int, default=0, help='the number of months since the reference_time_string that this masking operation corresponds to')
    parser.add_argument('--analysis_directory', type=str, default="/work/scratch/cmorice/advanced_standard/", help='root directory of the analysis')
    parser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    parser.add_argument('--output_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='operation index at which the analysis grid produced')
    parser.add_argument('--constraint_threshold',  type=float, default=0.6, help='threshold at which the observation influence indicates that the analysis is constrained')
    parser.add_argument('--climatology_limit',  type=float, default=0.5, help='threshold for masking on climatology uncertainty')
    parser.add_argument('--largescale_limit',  type=float, default=0.5, help='threshold for masking on largescale uncertainty')
    args = parser.parse_args()

    flagging_operation( args.reference_time_string,
                        args.operation_index,
                        args.analysis_directory,
                        args.iteration,
                        args.output_directory,
                        args.climatology_limit,
                        args.largescale_limit,
                        args.constraint_threshold )

if __name__ == '__main__':
    
    main()
    
