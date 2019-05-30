"""Checks for constraint on latent variables that are projected onto the grid"""


import argparse
import os.path
import numpy
from dateutil.relativedelta import relativedelta
import datetime
from dateutil import parser

import iris

from eustace.outputformats.outputvariable import OutputVariable
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure
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
MISSING_MARINE_FLAG  = flags.AREAL_LOW_FLAG

# missing data indicator
MISSING_FLAG_FLAG    = flags.MISSING_FLAG_FLAG

FLAG_MEANINGS = flags.FLAG_MEANINGS

# parameters for uncertainty thresholding
CLIMATOLOGY_LATENT_FLAG  = flags.CLIMATOLOGY_LATENT_FLAG
LARGE_SCALE_LATENT_FLAG   = flags.LARGE_SCALE_LATENT_FLAG
CONSTRAINT_THRESHOLD = 0.2#6

"""

Analysis gridder like operation to get projection of latent variable constraints onto each grid cell

"""

from eustace.analysis.advanced_standard.examples.moving_climatology import AnalysisSystem_EUSTACE
import eustace.analysis.advanced_standard.examples.model_structure
from eustace.analysis.advanced_standard.components.storage_files_batch import SpaceTimeComponentSolutionStorageBatched_Files
import eustace.timeutils.epoch
from eustace.analysis.advanced_standard.fileio.output_projector import Projector
from eustace.outputformats import definitions

import eustace.timeutils.epoch

def latent_variable_flag(input_directory, output_directory, iteration, processing_dates):
    
    # manually setup the analysis model for the R1413 run - Warning: the eustace svn revision must be correct for the global bias model interpretation to that run analysis
    
    storage_climatology = SpaceTimeComponentSolutionStorageBatched_Files( statefilename_read='/work/scratch/cmorice/advanced_standard/climatology_solution_9/climatology_solution_9.pickle',
                                                                          sample_filename_read='/work/scratch/cmorice/advanced_standard/climatology_solution_sample_9/climatology_solution_sample_9.pickle',
                                                                          prior_sample_filename_read='/work/scratch/cmorice/advanced_standard/climatology_solution_prior_sample_9/climatology_solution_prior_sample_9.pickle',
                                                                          keep_in_memory = True )
    
    storage_large_scale = SpaceTimeComponentSolutionStorageBatched_Files( statefilename_read='/work/scratch/cmorice/advanced_standard/large_scale_solution_9/large_scale_solution_9.pickle',
                                                                          sample_filename_read='/work/scratch/cmorice/advanced_standard/large_scale_solution_sample_9/large_scale_solution_sample_9.pickle',
                                                                          prior_sample_filename_read='/work/scratch/cmorice/advanced_standard/large_scale_solution_prior_sample_9/large_scale_solution_prior_sample_9.pickle',
                                                                          keep_in_memory = True )
                                                                          
    storage_local = eustace.analysis.advanced_standard.components.storage_files_batch.SpatialComponentSolutionStorageIndexed_Files()
    covariates_descriptor = "/gws/nopw/j04/eustace/data/internal/climatology_covariates/covariates.json"
    insitu_biases = True
    breakpoints_file = "/gws/nopw/j04/eustace/data/internal/D1.7/daily/eustace_stations_global_R001127_daily_status.nc"
    global_biases = True
    global_biases_group_list = ["surfaceairmodel_ice_global" , "surfaceairmodel_land_global", "surfaceairmodel_ocean_global"]
    compute_uncertainties = False
    method = 'EXACT'
    compute_sample = False
    sample_size = definitions.GLOBAL_SAMPLE_SHAPE[3]
    compute_prior_sample = False


    print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method)
    
    
    grid_resolution = [180. / definitions.GLOBAL_FIELD_SHAPE[1], 360. / definitions.GLOBAL_FIELD_SHAPE[2]]
    
    latitudes=numpy.linspace(-90.+grid_resolution[0]/2., 90.-grid_resolution[0]/2, num=definitions.GLOBAL_FIELD_SHAPE[1])
    longitudes=numpy.linspace(-180.+grid_resolution[1]/2., 180.-grid_resolution[1]/2, num=definitions.GLOBAL_FIELD_SHAPE[2])
    
    #timebase = TimeBaseDays(eustace.timeutils.epoch.EPOCH)
    #processdates = [timebase.number_to_datetime(daynumber) for daynumber in time_indices]
    
    # get times as understood by the analysis sustem
    time_indices =[eustace.timeutils.epoch.days_since_epoch(t) for t in processing_dates]
    
    cell_sampling   = [1, 1]
    blocking = 10

    # thinned set of sample indices for inclusion in output product
    sample_indices = range(definitions.GLOBAL_SAMPLE_SHAPE[3])
    
    climatology_projector = None
    large_scale_projector = None
    local_projector = None

    
    for ( inner_index, time_index, processdate ) in zip( range(len(time_indices)), time_indices, processing_dates ):
        print time_index
        
        # initialise flags
        flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )
        
        # Configure output grid
        outputstructure = OutputRectilinearGridStructure(time_index, processdate,
                                                     latitudes=latitudes,
                                                     longitudes=longitudes)
        
        # climatology component
        print 'Evaluating: climatology'
        if climatology_projector is None:
            climatology_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            climatology_projector.set_component(analysissystem.components[0])
            
            latent_climatology_constraint = evaluate_latent_variable_constraint(climatology_projector)
        
        climatology_projector.update_time_index(time_index, keep_design = False)
        climatology_projector.evaluate_design_matrix()
        
        climatology_statistic = evaluate_constraint_statistic(climatology_projector, latent_climatology_constraint, CONSTRAINT_THRESHOLD).reshape(definitions.GLOBAL_FIELD_SHAPE[1:])

         

        flag_values[climatology_statistic] = flag_values[climatology_statistic] | CLIMATOLOGY_LATENT_FLAG
        
        # large scale component
        print 'Evaluating: large-scale'
        if large_scale_projector is None:
            large_scale_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            large_scale_projector.set_component(analysissystem.components[1])
            
            latent_large_scale_constraint = evaluate_latent_variable_constraint(large_scale_projector)
            
        large_scale_projector.update_time_index(time_index, keep_design = False)
        large_scale_projector.evaluate_design_matrix()
        
        large_scale_statistic = evaluate_constraint_statistic(large_scale_projector, latent_large_scale_constraint, CONSTRAINT_THRESHOLD).reshape(definitions.GLOBAL_FIELD_SHAPE[1:])

        flag_values[large_scale_statistic] = flag_values[large_scale_statistic] | LARGE_SCALE_LATENT_FLAG
        
        outputfile = os.path.join(output_directory, '{:04d}'.format(processdate.year), 'eustace_analysis_{:d}_qc_flags_{:04d}{:02d}{:02d}.nc'.format(iteration, processdate.year, processdate.month, processdate.day))
        save_flag_file(flag_values, processdate, outputfile)

def evaluate_latent_variable_constraint(projector):
    # compute the ratio of posterior and prior uncertainty variances
    
    VARIANCE_RATIO_UPPER_BOUND = 1.0
    
    post_std = numpy.std( projector.state_samples, axis=1 )
    prior_std = numpy.std( projector.prior_samples, axis=1 )
    print "deviations post and prior:"
    print post_std
    print prior_std
    
    latent_variable_constraint = 1.0 - numpy.minimum( post_std**2 / prior_std**2, VARIANCE_RATIO_UPPER_BOUND)
    print latent_variable_constraint
    return latent_variable_constraint

def evaluate_constraint_statistic(projector, latent_variable_constraint, constraint_threshold):
    # Indicates where unconstrained variables project onto an output
    
    unconstraint_variable_indicator = numpy.zeros(latent_variable_constraint.shape)
    unconstraint_variable_indicator[latent_variable_constraint < constraint_threshold]  = 1.0
    
    print projector.design_matrix.shape
    print unconstraint_variable_indicator.shape
    projection = projector.design_matrix.dot(unconstraint_variable_indicator)
    print projection
    constraint_statistic = projection > 0.0
    
    return constraint_statistic
    
        
"""

Calls setup to run as batches in lsf

"""

def flagging_operation(  reference_time_string,
                            operation_index,
                            analysis_directory,
                            iteration,
                            output_directory ):
    
    # get dates in the month to be processed in this operation
    reference_time = parser.parse(reference_time_string)
    processing_dates = operation_dates(reference_time, operation_index)
    
    # check that directory exists and if not then make it
    if not os.path.exists(os.path.join(output_directory, str(processing_dates[0].year))):
        os.makedirs(os.path.join(output_directory, str(processing_dates[0].year)))
    
    latent_variable_flag(analysis_directory, output_directory, iteration, processing_dates)
    

def main():

    print 'Submission of advanced standard analysis observation constraint flagging jobs'
    
    parser = argparse.ArgumentParser(description='Monthly batches of masking operations')
    
    parser.add_argument('--reference_time_string', type=str, default="1880-01-01", help='reference first day to run formatted YYY-MM-DD. Should be first day of month.')
    parser.add_argument('--operation_index',  type=int, default=0, help='the number of months since the reference_time_string that this masking operation corresponds to')
    parser.add_argument('--analysis_directory', type=str, default="/work/scratch/cmorice/advanced_standard/", help='root directory of the analysis')
    parser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    parser.add_argument('--output_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='operation index at which the analysis grid produced')

    args = parser.parse_args()

    flagging_operation( args.reference_time_string,
                        args.operation_index,
                        args.analysis_directory,
                        args.iteration,
                        args.output_directory, )

if __name__ == '__main__':
    
    main()
