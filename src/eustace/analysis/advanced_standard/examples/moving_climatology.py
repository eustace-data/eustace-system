
import os
import psutil
import sys

from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.outputformats import definitions

from eumopps.version.svn import get_revision_id_for_module

from model_structure import ClimatologyDefinition, LargeScaleDefinition, LocalDefinition, NonStationaryLocalDefinition

from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.elements.bias_insitu_land import InsituLandBiasElement

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.observationsource import ObservationSource

from eustace.timeutils.decimaltime import datetime_to_decimal_year
from datetime import datetime
from dateutil import relativedelta
from dateutil.rrule import rrule, MONTHLY
import eustace.timeutils.epoch

import numpy
import copy

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

from eustace.analysis.observationsource import ObservationSource

from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eustace.outputformats import definitions

from inputloader_rawbinary import AnalysisSystemInputLoaderRawBinary_OneDay
from eustace.analysis.fileio.observationsource_missing import ObservationSourceMissing

from eumopps.timeutils.timebase import TimeBaseDays
from eumopps.timeutils import datetime_numeric


"""

AnalysisSystem classes for the advanced standard analysis and optimisation system

"""

class AnalysisSystem_EUSTACE(AnalysisSystem):
    """Analysis system in which the storage objects for each component are specified during construction."""
    
    def __init__(self, storage_climatology, storage_large_scale, storage_local, covariates_descriptor, 
               insitu_biases = False, breakpoints_file = None, global_biases = False, global_biases_group_list = [], 
               compute_uncertainties = False, method = 'EXACT',
               compute_sample = False, sample_size = definitions.GLOBAL_SAMPLE_SHAPE[3], compute_prior_sample = False):

        super(AnalysisSystem_EUSTACE, self).__init__(

            components = [ 
                SpaceTimeComponent(ClimatologyDefinition(covariates_descriptor), storage_climatology, True, compute_uncertainties, method, compute_sample, sample_size, compute_prior_sample),
                SpaceTimeComponent(LargeScaleDefinition(insitu_biases, breakpoints_file), storage_large_scale, True, compute_uncertainties, method, compute_sample, sample_size, compute_prior_sample),
                SpatialComponent(LocalDefinition(global_biases, global_biases_group_list), storage_local, compute_uncertainties, method, compute_sample, sample_size, compute_prior_sample)],

            observable = ObservationSource.TMEAN)

"""

Functions to be called from EUMOPPS for observation processing, analysis fitting and projecting onto the output grid

"""

def process_inputs(storage_climatology, storage_large_scale, storage_local, 
           inputsources, time_index, 
           component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
           compute_uncertainties, method,
           compute_sample, sample_size, compute_prior_sample):
    
    # Build inputloaders from list of sources
    inputloaders = [ AnalysisSystemInputLoaderRawBinary_OneDay(time_index=time_index, **source) for source in inputsources ]
    
    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method,
                                    compute_sample, sample_size, compute_prior_sample)
    
    # Build and store measurement systems
    analysissystem.update_component_time(inputloaders, component_index, time_index)

def process_inputs_batch(storage_climatology, storage_large_scale, storage_local, 
           inputsources, time_indices, 
           component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
           compute_uncertainties, method,
           compute_sample, sample_size, compute_prior_sample):

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method,
                                    compute_sample, sample_size, compute_prior_sample)
    
    # for each time step construct inputloaders for raw binary reading and then pass to the measurement update accumulator
    #time_indices = time_index
  
    # need to construct inputloaders = [ [inputloader_0,...,inputloader_n]_{t=0}, ..., [inputloader_0,...,inputloader_n]_{t=T} ]
    inputloaders = []
    
    
    for source in inputsources:
        print source
        #for variable in source:
            #for filename in source[variable]:
                #print filename
    
    print time_indices
    for i in range(len(time_indices)):
        
        # construct list of inputloaders for time index i
        substep_sources = []
        for source in inputsources:
            
            substep_source = {}
            
            for input_key in source:
                if input_key == 'name':
                    continue
                elif input_key == 'fixed_location_lookup_filename':
                    substep_source[input_key] = source[input_key][i]
                    print substep_source[input_key]
                else:
                    substep_source[input_key] = {}
                    for variable in source[input_key]:
                        substep_source[input_key][variable] = source[input_key][variable][i]
                        print substep_source[input_key][variable]
                        
            substep_sources.append( substep_source )
            
        inputloaders.append( [ AnalysisSystemInputLoaderRawBinary_OneDay(time_index=time_indices[i], **source_input) for source_input in substep_sources ] )
        
    # accumulate the measurement updates and store
    analysissystem.update_component_times(inputloaders, component_index, time_indices)
    
def solve(storage_climatology, storage_large_scale, storage_local, 
      component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
      compute_uncertainties, method,
          compute_sample, sample_size, compute_prior_sample):

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                        compute_uncertainties, method,
                        compute_sample, sample_size, compute_prior_sample)

    # Solve
    analysissystem.update_component_solution(component_index)

def measurement_merge(storage_climatology, storage_large_scale, storage_local, 
      component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
      compute_uncertainties, method, output_file):

    import pickle

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                        compute_uncertainties, method,
          compute_sample, sample_size)

    # Get the merged update object
    merged_updates = analysissystem.components[component_index].component_solution().merge_updates()
    
    # Save to disk
    with open(output_file, 'wb') as f:
        pickle.dump(merged_updates, f)

#build_output_projectors
def output_grid(storage_climatology, storage_large_scale, storage_local, 
            outputfile, climatologyfile, largescalefile, localfile,
            processdate, time_index, 
            covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
            compute_uncertainties, method,
            compute_sample, sample_size, compute_prior_sample):

    from eustace.analysis.advanced_standard.fileio.output_projector import Projector

    print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method)
    
    
    grid_resolution = [180. / definitions.GLOBAL_FIELD_SHAPE[1], 360. / definitions.GLOBAL_FIELD_SHAPE[2]]
    
    latitudes=numpy.linspace(-90.+grid_resolution[0]/2., 90.-grid_resolution[0]/2, num=definitions.GLOBAL_FIELD_SHAPE[1])
    longitudes=numpy.linspace(-180.+grid_resolution[1]/2., 180.-grid_resolution[1]/2, num=definitions.GLOBAL_FIELD_SHAPE[2])
    
    cell_sampling   = [1, 1]
    blocking = 10

    # Configure output grid
    outputstructure = OutputRectilinearGridStructure(time_index, processdate,
                                                     latitudes=latitudes,
                                                     longitudes=longitudes)


    # thinned set of sample indices for inclusion in output product
    sample_indices = range(definitions.GLOBAL_SAMPLE_SHAPE[3])
    
    # climatology component
    print 'Evaluating: climatology'
    climatology_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
    climatology_projector.set_component(analysissystem.components[0])
    climatology_projector.evaluate_design_matrix()

    climatology_expected_value = climatology_projector.project_expected_value().reshape((-1,1))
    climatology_uncertainties = climatology_projector.project_sample_deviation()
    climatology_samples = climatology_projector.project_sample_values(sample_indices = sample_indices) + climatology_expected_value
    climatology_unconstraint = climatology_uncertainties**2 / climatology_projector.project_sample_deviation(prior = True)**2
    
    climatology_projector = None # clear projector from memory
    
    print climatology_expected_value.shape, climatology_uncertainties.shape, climatology_samples.shape
    
    # large scale component
    print 'Evaluating: large-scale'
    large_scale_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
    large_scale_projector.set_component(analysissystem.components[1])
    large_scale_projector.evaluate_design_matrix()
    
    large_scale_expected_value = large_scale_projector.project_expected_value().reshape((-1,1))
    large_scale_uncertainties = large_scale_projector.project_sample_deviation()
    large_scale_samples = large_scale_projector.project_sample_values(sample_indices = sample_indices) + large_scale_expected_value
    large_scale_unconstraint = large_scale_uncertainties**2 / large_scale_projector.project_sample_deviation(prior = True)**2
    
    large_scale_projector = None # clear projector from memory
    
    print large_scale_expected_value.shape, large_scale_uncertainties.shape, large_scale_samples.shape
    
    # local component
    print 'Evaluating: local'
    local_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
    local_projector.set_component(analysissystem.components[2])
    local_projector.evaluate_design_matrix()

    local_expected_value = local_projector.project_expected_value().reshape((-1,1))
    local_uncertainties = local_projector.project_sample_deviation()
    local_samples = local_projector.project_sample_values(sample_indices = sample_indices) + local_expected_value    
    local_unconstraint = local_uncertainties**2 / local_projector.project_sample_deviation(prior = True)**2
    
    local_projector = None # clear projector from memory
    
    print local_expected_value.shape, local_uncertainties.shape, local_samples.shape
    
    # Save results
    print outputfile
    # main merged product output files
    filebuilder = FileBuilderGlobalField(
        outputfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional output',
        __name__, 
        '')
    
    
    climatology_fraction = local_unconstraint # defined as ratio of posterior to prior variance in local component
    
    result_expected_value = climatology_expected_value + large_scale_expected_value + local_expected_value
    result_expected_uncertainties = numpy.sqrt(climatology_uncertainties**2 + large_scale_uncertainties**2 + local_uncertainties**2)
    
    filebuilder.add_global_field(definitions.TAS, result_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, result_expected_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_CLIMATOLOGY_FRACTION, climatology_fraction.reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    for index in range(definitions.GLOBAL_SAMPLE_SHAPE[3]):
        variable = copy.deepcopy(definitions.TASENSEMBLE)
        variable.name = variable.name + '_' + str(index)
        selected_sample = (climatology_samples[:,index] + large_scale_samples[:,index] + local_samples[:,index]).ravel()
        filebuilder.add_global_field(variable, selected_sample.reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    filebuilder.save_and_close()
    
    # climatology only output
    filebuilder = FileBuilderGlobalField(
        climatologyfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional component output - climatology',
        __name__, 
        '')
    
    filebuilder.add_global_field(definitions.TAS, climatology_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, climatology_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_CLIMATOLOGY_FRACTION, climatology_unconstraint.reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    # large scale only output
    filebuilder = FileBuilderGlobalField(
        largescalefile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional component output - large scale',
        __name__, 
        '')
    
    filebuilder.add_global_field(definitions.TASPERTURBATION, large_scale_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, large_scale_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_CLIMATOLOGY_FRACTION, large_scale_unconstraint.reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    # local only output
    filebuilder = FileBuilderGlobalField(
        localfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional component output - local',
        __name__, 
        '')

    filebuilder.add_global_field(definitions.TASPERTURBATION, local_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, local_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_CLIMATOLOGY_FRACTION, local_unconstraint.reshape(definitions.GLOBAL_FIELD_SHAPE))


#build_output_projectors
def output_grid_batch(storage_climatology, storage_large_scale, storage_local, 
                        outputfiles, climatologyfiles, largescalefiles, localfiles,
                        time_indices, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method,
                        compute_sample, sample_size, compute_prior_sample):

    from eustace.analysis.advanced_standard.fileio.output_projector import Projector

    variance_ratio_upper_bound = 1.0

    print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method)
    
    
    grid_resolution = [180. / definitions.GLOBAL_FIELD_SHAPE[1], 360. / definitions.GLOBAL_FIELD_SHAPE[2]]
    
    latitudes=numpy.linspace(-90.+grid_resolution[0]/2., 90.-grid_resolution[0]/2, num=definitions.GLOBAL_FIELD_SHAPE[1])
    longitudes=numpy.linspace(-180.+grid_resolution[1]/2., 180.-grid_resolution[1]/2, num=definitions.GLOBAL_FIELD_SHAPE[2])
    
    timebase = TimeBaseDays(eustace.timeutils.epoch.EPOCH)
    #processdates = [datetime_numeric.build( timebase.number_to_datetime(daynumber) ) for daynumber in time_indices]
    processdates = [timebase.number_to_datetime(daynumber) for daynumber in time_indices]
    
    cell_sampling   = [1, 1]
    blocking = 10

    # thinned set of sample indices for inclusion in output product
    sample_indices = range(definitions.GLOBAL_SAMPLE_SHAPE[3])
    
    climatology_projector = None
    large_scale_projector = None
    local_projector = None
    
    for ( inner_index, time_index, processdate ) in zip( range(len(time_indices)), time_indices, processdates ):
        print time_index
        # Configure output grid
        outputstructure = OutputRectilinearGridStructure(time_index, processdate,
                                                     latitudes=latitudes,
                                                     longitudes=longitudes)
        
        # climatology component
        print 'Evaluating: climatology'
        if climatology_projector is None:
            climatology_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            climatology_projector.set_component(analysissystem.components[0])
        
        climatology_projector.update_time_index(time_index, keep_design = False)
        climatology_projector.evaluate_design_matrix()

        climatology_expected_value = climatology_projector.project_expected_value().reshape((-1,1))
        climatology_uncertainties = climatology_projector.project_sample_deviation()
        climatology_samples = climatology_projector.project_sample_values(sample_indices = sample_indices) + climatology_expected_value
        climatology_unconstraint = numpy.minimum( climatology_uncertainties**2 / climatology_projector.project_sample_deviation(prior = True)**2, variance_ratio_upper_bound )
        
        # large scale component
        print 'Evaluating: large-scale'
        if large_scale_projector is None:
            large_scale_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            large_scale_projector.set_component(analysissystem.components[1])
            
        large_scale_projector.update_time_index(time_index, keep_design = False)
        large_scale_projector.evaluate_design_matrix()
        
        large_scale_expected_value = large_scale_projector.project_expected_value().reshape((-1,1))
        large_scale_uncertainties = large_scale_projector.project_sample_deviation()
        large_scale_samples = large_scale_projector.project_sample_values(sample_indices = sample_indices) + large_scale_expected_value
        large_scale_unconstraint = numpy.minimum( large_scale_uncertainties**2 / large_scale_projector.project_sample_deviation(prior = True)**2, variance_ratio_upper_bound)
        
        
        
        # local component - time handling updates state to new time but does not recompute the design matrix
        print 'Evaluating: local'
        if local_projector is None:
            local_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            local_projector.set_component(analysissystem.components[2])
            local_projector.evaluate_design_matrix()
        else:
            local_projector.update_time_index(time_index, keep_design = True)
            local_projector.set_component(analysissystem.components[2], keep_design = True)
        
        print analysissystem.components
        
        local_expected_value = local_projector.project_expected_value().reshape((-1,1))
        local_uncertainties = local_projector.project_sample_deviation()
        local_samples = local_projector.project_sample_values(sample_indices = sample_indices) + local_expected_value  

        local_unconstraint = numpy.minimum( local_uncertainties**2 / local_projector.project_sample_deviation(prior = True)**2, variance_ratio_upper_bound )
        
        
        # Save results
        outputfile = outputfiles[inner_index]
        print outputfile
        # main merged product output files
        filebuilder = FileBuilderGlobalField(
            outputfile, 
            eustace.timeutils.epoch.days_since_epoch(processdate),
            'EUSTACE Analysis',
            get_revision_id_for_module(eustace),
            definitions.TAS.name,
            '',
            'Provisional output',
            __name__, 
            '')
        
        
        result_expected_value = climatology_expected_value + large_scale_expected_value + local_expected_value
        result_expected_uncertainties = numpy.sqrt(climatology_uncertainties**2 + large_scale_uncertainties**2 + local_uncertainties**2)
        climatology_fraction = local_unconstraint # defined as ratio of posterior to prior variance in local component
        
        filebuilder.add_global_field(definitions.TAS, result_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TASUNCERTAINTY, result_expected_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, 1.0 - climatology_fraction.reshape(definitions.GLOBAL_FIELD_SHAPE))
        
        for index in range(definitions.GLOBAL_SAMPLE_SHAPE[3]):
            variable = copy.deepcopy(definitions.TASENSEMBLE)
            variable.name = variable.name + '_' + str(index)
            selected_sample = (climatology_samples[:,index] + large_scale_samples[:,index] + local_samples[:,index]).ravel()
            filebuilder.add_global_field(variable, selected_sample.reshape(definitions.GLOBAL_FIELD_SHAPE))
        
        filebuilder.save_and_close()
        
        # climatology only output
        climatologyfile = climatologyfiles[inner_index]
        filebuilder = FileBuilderGlobalField(
            climatologyfile, 
            eustace.timeutils.epoch.days_since_epoch(processdate),
            'EUSTACE Analysis',
            get_revision_id_for_module(eustace),
            definitions.TAS.name,
            '',
            'Provisional component output - climatology',
            __name__, 
            '')
        
        filebuilder.add_global_field(definitions.TAS, climatology_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TASUNCERTAINTY, climatology_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, 1.0 - climatology_unconstraint.reshape(definitions.GLOBAL_FIELD_SHAPE))
        
        filebuilder.save_and_close()
        
        # large scale only output
        largescalefile = largescalefiles[inner_index]
        filebuilder = FileBuilderGlobalField(
            largescalefile, 
            eustace.timeutils.epoch.days_since_epoch(processdate),
            'EUSTACE Analysis',
            get_revision_id_for_module(eustace),
            definitions.TAS.name,
            '',
            'Provisional component output - large scale',
            __name__, 
            '')
        
        filebuilder.add_global_field(definitions.TASPERTURBATION, large_scale_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TASUNCERTAINTY, large_scale_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, 1.0 - large_scale_unconstraint.reshape(definitions.GLOBAL_FIELD_SHAPE))
        
        filebuilder.save_and_close()
        
        # local only output
        localfile = localfiles[inner_index]
        filebuilder = FileBuilderGlobalField(
            localfile, 
            eustace.timeutils.epoch.days_since_epoch(processdate),
            'EUSTACE Analysis',
            get_revision_id_for_module(eustace),
            definitions.TAS.name,
            '',
            'Provisional component output - local',
            __name__, 
            '')

        filebuilder.add_global_field(definitions.TASPERTURBATION, local_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TASUNCERTAINTY, local_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, 1.0 - local_unconstraint.reshape(definitions.GLOBAL_FIELD_SHAPE))
        
        filebuilder.save_and_close()
        
        print "Memory usage (MB):", psutil.Process(os.getpid()).memory_info().rss / (1024*1024)

def early_look_grid_batch(storage_climatology, storage_large_scale, storage_local, 
                        outputfiles,
                        time_indices, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method,
                        compute_sample, sample_size, compute_prior_sample):
    """Produce 'early look' NetCDF output files without loading or gridding uncertainty information
    
    For inspection of analysis output prior to the final gridding step.
    
    """
    
    from eustace.analysis.advanced_standard.fileio.output_projector import Projector

    print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method)
    
    
    grid_resolution = [180. / definitions.GLOBAL_FIELD_SHAPE[1], 360. / definitions.GLOBAL_FIELD_SHAPE[2]]
    
    latitudes=numpy.linspace(-90.+grid_resolution[0]/2., 90.-grid_resolution[0]/2, num=definitions.GLOBAL_FIELD_SHAPE[1])
    longitudes=numpy.linspace(-180.+grid_resolution[1]/2., 180.-grid_resolution[1]/2, num=definitions.GLOBAL_FIELD_SHAPE[2])
    
    timebase = TimeBaseDays(eustace.timeutils.epoch.EPOCH)
    
    processdates = [timebase.number_to_datetime(daynumber) for daynumber in time_indices]
    
    cell_sampling   = [1, 1]
    blocking = 10

    # thinned set of sample indices for inclusion in output product
    
    climatology_projector = None
    large_scale_projector = None
    local_projector = None
    
    for ( inner_index, time_index, processdate ) in zip( range(len(time_indices)), time_indices, processdates ):
        print time_index
        # Configure output grid
        outputstructure = OutputRectilinearGridStructure(time_index, processdate,
                                                     latitudes=latitudes,
                                                     longitudes=longitudes)
        
        # climatology component
        print 'Evaluating: climatology'
        if climatology_projector is None:
            climatology_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            climatology_projector.set_component(analysissystem.components[0])
        
        climatology_projector.update_time_index(time_index, keep_design = False)
        climatology_projector.evaluate_design_matrix()

        climatology_expected_value = climatology_projector.project_expected_value().reshape((-1,1))    
        
        # large scale component
        print 'Evaluating: large-scale'
        if large_scale_projector is None:
            large_scale_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            large_scale_projector.set_component(analysissystem.components[1])
            
        large_scale_projector.update_time_index(time_index, keep_design = False)
        large_scale_projector.evaluate_design_matrix()
        
        large_scale_expected_value = large_scale_projector.project_expected_value().reshape((-1,1))        
        
        
        # local component - time handling updates state to new time but does not recompute the design matrix
        print 'Evaluating: local'
        if local_projector is None:
            local_projector = Projector(latitudes, longitudes, grid_resolution, time_index, cell_sampling, blocking)
            local_projector.set_component(analysissystem.components[2])
            local_projector.evaluate_design_matrix()
        else:
            local_projector.update_time_index(time_index, keep_design = True)
            local_projector.set_component(analysissystem.components[2], keep_design = True)
        
        print analysissystem.components
        
        local_expected_value = local_projector.project_expected_value().reshape((-1,1))
        
        
        # Save results
        outputfile = outputfiles[inner_index]
        print outputfile
        # main merged product output files
        filebuilder = FileBuilderGlobalField(
            outputfile, 
            eustace.timeutils.epoch.days_since_epoch(processdate),
            'EUSTACE Analysis',
            get_revision_id_for_module(eustace),
            definitions.TAS.name,
            '',
            'Provisional output',
            __name__, 
            '')
        
        
        field_definition_tas = definitions.OutputVariable.from_template(definitions.TEMPLATE_TEMPERATURE, 'tas', quantity='average', cell_methods='time: mean')
        field_definition_tas_climatology = definitions.OutputVariable.from_template(definitions.TEMPLATE_TEMPERATURE, 'tas_climatology', quantity='average', cell_methods='time: mean')
        field_definition_tas_large_scale = definitions.OutputVariable.from_template(definitions.TEMPLATE_PERTURBATION, 'tas_large_scale', quantity='average', cell_methods='time: mean')
        field_definition_tas_daily_local = definitions.OutputVariable.from_template(definitions.TEMPLATE_PERTURBATION, 'tas_daily_local', quantity='average', cell_methods='time: mean')
        
        result_expected_value = climatology_expected_value + large_scale_expected_value + local_expected_value
        filebuilder.add_global_field(field_definition_tas, result_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(field_definition_tas_climatology, climatology_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(field_definition_tas_large_scale, large_scale_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(field_definition_tas_daily_local, local_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.save_and_close()
        
        print "Memory usage (MB):", psutil.Process(os.getpid()).memory_info().rss / (1024*1024)

def output_homog(storage_climatology, storage_large_scale, storage_local, 
            outputfile, climatologyfile, largescalefile, localfile,
            processdate, time_index, 
            covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
            compute_uncertainties, method,
            compute_sample, sample_size, compute_prior_sample):

    from eustace.preprocess.fileio.insitu_land_breakpoints import ObservationBreakPointSourceInsituLand
    from eustace.analysis.advanced_standard.fileio.output_projector import Projector

    print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method)


    component = analysissystem.components[1]
    element = component.storage.element_read()
    hyperparameters = component.storage.hyperparameters_read()
    
    currentstate = component.solutionstorage.state_read()
    element_states = component.storage.element_read().element_prior(hyperparameters).element_states(currentstate)
    
    breakpoints_reader = ObservationBreakPointSourceInsituLand(breakpoints_file)
    station_locations = breakpoints_reader.observation_location_lookup()
    #station_indices = breakpoints_reader.observation_location_lookup()
    breakpoints_reader.dataset.close()
    
    for i, subelement in enumerate(element.combination):

        if isinstance(subelement, BiasElement):
            mystate = element_states[i]
            
            if hasattr(subelement, 'observed_breakpoints'):
                import pickle
                
                outdict = {'breakpoints': subelement.observed_breakpoints,
                           'locations':   station_locations.T,
                           'biases':      mystate}
                
                with open('/work/scratch/cmorice/advanced_standard/observed_breakpoints.pickle', 'w') as f:
                    pickle.dump(outdict, f)
