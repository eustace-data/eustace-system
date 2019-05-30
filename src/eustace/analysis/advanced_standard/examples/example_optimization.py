import numpy
import os.path

from model_structure import ShortScaleSetup
from moving_climatology import AnalysisSystem_EUSTACE, ClimatologyDefinition, LargeScaleDefinition
from eustace.analysis.advanced_standard.optimizationsystem import OptimizationSystem
from eustace.analysis.advanced_standard.elements.local_view import LocalSubRegion, LocalSuperTriangle, ExtendedCombinationHyperparameters, SphereMeshViewGlobal
from eustace.analysis.advanced_standard.elements.local_view import NonStationaryLocal, ExpandedLocalHyperparameters
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.elements.combination import CombinationElement, CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement, LocalHyperparameters
from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.spatialdelayed import DelayedSpatialComponent
from eustace.analysis.observationsource import ObservationSource
from eustace.outputformats import definitions

"""

AnalysisSystem setups required for the optimization

"""

class GlobalOptimizationSystem_EUSTACE(OptimizationSystem):
    """An OptimizationSystem with a global scale local component
    
    Used as an intermediate step towards construction of OptimizationSystem 
    for regional optimisation.
    
    """
    
    def __init__(self, storage_climatology, storage_large_scale, storage_local_bias, storage_local_spde,
               covariates_descriptor, insitu_biases = False, breakpoints_file = None, global_biases = False, global_biases_group_list = [], 
               compute_uncertainties = False, method = 'EXACT',
               compute_sample = False, sample_size = definitions.GLOBAL_SAMPLE_SHAPE[3]):
                   
        # initialise the OptimizationSystem
        super(GlobalOptimizationSystem_EUSTACE, self).__init__(
                        components = [ 
                                        SpaceTimeComponent(ClimatologyDefinition(covariates_descriptor), storage_climatology, True, compute_uncertainties, method, compute_sample, sample_size),
                                        SpaceTimeComponent(LargeScaleDefinition(insitu_biases, breakpoints_file), storage_large_scale, True, compute_uncertainties, method, compute_sample, sample_size),
                                        SpatialComponent(PureBiasComponentDefinition(global_biases, global_biases_group_list), storage_local_bias, compute_uncertainties, method, compute_sample, sample_size),
                                        DelayedSpatialComponent(PureLocalComponentDefinition(), storage_local_spde, compute_uncertainties, method, compute_sample, sample_size)
                                    ],
                        observable = ObservationSource.TMEAN )

class RegionOptimizationSystem_EUSTACE(OptimizationSystem):
    """Analysis system in which the storage objects for each component are specified during construction."""
    
    def __init__(self, storage_climatology, storage_large_scale, storage_local_bias, storage_region_spde,
               covariates_descriptor, insitu_biases = False, breakpoints_file = None, global_biases = False, global_biases_group_list = [], 
               compute_uncertainties = False, method = 'EXACT',
               compute_sample = False, sample_size = definitions.GLOBAL_SAMPLE_SHAPE[3],
               neighbourhood_level = 0, region_index = 0, regionspec = 'LocalSubRegion'):

        # initialise the OptimizationSystem
        super(RegionOptimizationSystem_EUSTACE, self).__init__(
                        components = [ 
                                        SpaceTimeComponent(ClimatologyDefinition(covariates_descriptor), storage_climatology, True, compute_uncertainties, method, compute_sample, sample_size),
                                        SpaceTimeComponent(LargeScaleDefinition(insitu_biases, breakpoints_file), storage_large_scale, True, compute_uncertainties, method, compute_sample, sample_size),
                                        SpatialComponent(PureBiasComponentDefinition(global_biases, global_biases_group_list), storage_local_bias, compute_uncertainties, method, compute_sample, sample_size),
                                        DelayedSpatialComponent(LocalViewDefinition( neighbourhood_level, region_index, regionspec ), storage_region_spde, compute_uncertainties, method, compute_sample, sample_size)
                                    ],
                        observable = ObservationSource.TMEAN )

"""

Optimisation Component Definitions

"""
                                                            
class LocalViewDefinition(ComponentStorage_InMemory):
    """Define local component."""
    
    def __init__(self, neighbourhood_level, region_index, regionspec = 'LocalSubRegion'):         
    
        setup = ShortScaleSetup()

        if regionspec == 'LocalSubRegion':
            elementclass = LocalSubRegion
        elif regionspec == 'LocalSuperTriangle':
            elementclass = LocalSuperTriangle

        super(LocalViewDefinition, self).__init__(
            CombinationElement([elementclass(setup.local_settings.n_triangulation_divisions, neighbourhood_level, region_index)]),
            ExtendedCombinationHyperparameters([LocalHyperparameters(numpy.log(setup.local_settings.amplitude), 
                                                            numpy.log(numpy.radians(setup.local_settings.space_length_scale)))]))

class PureLocalComponentDefinition(ComponentStorage_InMemory):
    """Define a component that only contains the daily local random field model."""
    
    def __init__(self):         
    
        setup = ShortScaleSetup()

        super(PureLocalComponentDefinition, self).__init__(
            CombinationElement([LocalElement(setup.local_settings.n_triangulation_divisions)]),     
            ExtendedCombinationHyperparameters([LocalHyperparameters(numpy.log(setup.local_settings.amplitude), 
                                                            numpy.log(numpy.radians(setup.local_settings.space_length_scale)))]))

class PureBiasComponentDefinition(ComponentStorage_InMemory):
    """Define a component that only contains the daily bias model."""
    
    def __init__(self, bias_terms = False, global_biases_group_list = []):         
    
        setup = ShortScaleSetup()
        
        if bias_terms:
            bias_elements = [BiasElement(groupname, 1) for groupname in global_biases_group_list]
            bias_hyperparameters = [CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude)) for index in range(len(global_biases_group_list))]
        else:
            bias_elements, bias_hyperparameters = [], []
        
        super(PureBiasComponentDefinition, self).__init__(
            CombinationElement(bias_elements),         
            CombinationHyperparameters(bias_hyperparameters))

    
"""

Extracting states for component splits and local spde views

"""
                                                            
def split_states_time( full_component, optimisation_component, non_optimisation_component, element_optimisation_flags, time_key ):
    """Assign the state vector for a full component into optimisation and non-optimisation components"""
    
    full_state = full_component.solutionstorage.partial_state_read(time_key)
    element_states = full_component.storage.element_read().element_prior(full_component.storage.hyperparameters_read()).element_states(full_state)
    
    assert(len(element_states) == len(element_optimisation_flags))
    
    optimisation_states = []
    non_optimisation_states = []
    
    for element_index, optimisation_flag in enumerate(element_optimisation_flags):
        
        if optimisation_flag:
            optimisation_states.append(element_states[element_index])
        else:
            non_optimisation_states.append(element_states[element_index])
    print optimisation_states, non_optimisation_states
    
    # handling for un-initialised states
    if any(elem is None for elem in optimisation_states):
        optimisation_states = None
    else:
        optimisation_states = numpy.concatenate( optimisation_states )
        
    if any(elem is None for elem in non_optimisation_states):
        non_optimisation_states = None
    else:
        non_optimisation_states =numpy.concatenate( non_optimisation_states )
        
    optimisation_component.solutionstorage.partial_state_write( optimisation_states, time_key )
    non_optimisation_component.solutionstorage.partial_state_write( non_optimisation_states, time_key )
    
    
def extract_local_view_states_time( full_component, view_component, view_flags, time_key ):
    """Extract state from global component for a local view where view_flags indicates which sub-elements are local views"""
    
    full_state = full_component.solutionstorage.partial_state_read(time_key)
    element_states = full_component.storage.element_read().element_prior(full_component.storage.hyperparameters_read()).element_states(full_state)
    view_states = []
    for element_index, view_flag in enumerate(view_flags):
        if view_flag:
            view_state_indices = view_component.storage.element_read().combination[element_index].spde.get_active_vertex_indices()
            view_states.append( element_states[element_index][view_state_indices] )
        else:
            view_states.append(element_states[element_index])
            
    view_component.solutionstorage.partial_state_write(  numpy.concatenate( view_states ), time_key )
    

"""

Splitting analysis system into optimisable components and running the optimisation

"""

def run_split_states(storage_climatology, storage_large_scale, storage_local,
                     storage_local_bias, storage_local_spde,
                     covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                     compute_uncertainties, method,
                     compute_sample, sample_size,
                     time_index):
    
    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method, compute_sample, sample_size)
   
    # Build split system
    splitsystem = GlobalOptimizationSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local_bias, storage_local_spde,
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method, compute_sample, sample_size)

    # Only want to optimise the spde sub-element in analysissystem's local component - this should be the first sub element of component index 2
    element_optimisation_flags = [True,] + [False,]*(len(analysissystem.components[2].storage.element.combination) - 1)

    # Define components to be split by element according to element_optimisation_flags
    full_component = analysissystem.components[2]

    non_optimisation_component = splitsystem.components[2]
    optimisation_component = splitsystem.components[3]

    # Split state between optimisation and non optimisation components at specified time key
    split_states_time( full_component, optimisation_component, non_optimisation_component, element_optimisation_flags, time_index )


def optimize_region(storage_climatology, storage_large_scale, storage_local_bias, storage_region_spde,
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                        compute_uncertainties, method,
                        compute_sample, sample_size,
                        neighbourhood_level, region_index, regionspec,
                        inputdescriptor, time_keys,
                        hyperparameter_storage_file):
                        
    # Build regional system for optimisation
    regionsystem = RegionOptimizationSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local_bias, storage_region_spde,
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method, compute_sample, sample_size, neighbourhood_level, region_index, regionspec)
    
    optimization_component_index = 3
    
    # Build projected observation information and precision matrices as 
    # preprocessing step and solve to get optimization set point.
    
    regionsystem.update_component(inputdescriptor, optimization_component_index, time_keys)
    
    # run the optimization of the regional component hyperparameters
    regionsystem.optimize_component(optimization_component_index, hyperparameter_storage_file)

"""

Construction of non-stationary hyperparameter files

"""

def merge_local_parameterisations(full_resolution_level, neighbourhood_level, hyperparameter_filenames, output_filename):
    """Merge local hyperparameter fits by averaging over local fits at each latent variable"""
    
    sigma_accumulator = None
    rho_accumulator = None
    contribution_counter = None
    
    # step over regions summing local (sigma, rho) values
    for region_index, hyperparameter_file in enumerate(hyperparameter_filenames):
        
        if not os.path.isfile(hyperparameter_file):
            continue
        
        local_region = LocalViewDefinition(neighbourhood_level, region_index, regionspec = 'LocalSuperTriangle')
        local_region.hyperparameters.values_from_npy_savefile(hyperparameter_file)
        
        
        local_spde = local_region.element.combination[0].spde
        local_hyperparameters = local_region.hyperparameters.get_array()
        
        accumulators = SphereMeshViewGlobal.accumulate_local_parameterisations(sigma_accumulator,
                                                                                rho_accumulator,
                                                                                contribution_counter,
                                                                                local_spde,
                                                                                local_hyperparameters)
                                                                                                                                        
        sigma_accumulator, rho_accumulator, contribution_counter = accumulators
        

    # compute the hyperparameters as the logarithm of the average of local sigma and rho estimates at each node
    log_sigmas, log_rhos = SphereMeshViewGlobal.finalise_local_parameterisation_sigma_rho(sigma_accumulator,
                                                                                            rho_accumulator,
                                                                                            contribution_counter)

    # store the merged hyperparameter values to disk
    setup = ShortScaleSetup()    
    nonstationary_element = NonStationaryLocal(setup.local_settings.n_triangulation_divisions)
    
    nonstationary_hyperparameters =  ExtendedCombinationHyperparameters([ExpandedLocalHyperparameters(log_sigma = log_sigmas, log_rho   = log_rhos)])
    nonstationary_hyperparameters.values_to_npy_savefile(output_filename)
    
def initialise_local_hyperparameter_file(filename):
    """Write a non-stationary hyperparameter file for initial stationary parameters"""
    
    setup = ShortScaleSetup()
    
    local_element = NonStationaryLocal(setup.local_settings.n_triangulation_divisions)
    n_latent_variables = local_element.spde.n_latent_variables()
    
    hyperparameters =  ExtendedCombinationHyperparameters([ExpandedLocalHyperparameters(log_sigma = numpy.ones(n_latent_variables) * numpy.log(setup.local_settings.amplitude),
                                                                                        log_rho   = numpy.ones(n_latent_variables) * numpy.log(numpy.radians(setup.local_settings.space_length_scale)))])
    hyperparameters.values_to_npy_savefile(filename)
    
