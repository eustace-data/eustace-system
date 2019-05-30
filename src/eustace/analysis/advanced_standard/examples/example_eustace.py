"""Example infilled analysis using EUSTACE data from EUMOPPS catalogue."""

import copy
import eustace.timeutils.epoch
import numpy
import os.path

from datetime import datetime
from dateutil import relativedelta
from dateutil.rrule import rrule, MONTHLY

from eumopps.version.svn import get_revision_id_for_module

from eustaceconfig import WORKSPACE_PATH

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.elements.bias_insitu_land import InsituLandBiasElement
from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import LoadCovariateElement
from eustace.analysis.advanced_standard.elements.kronecker import SpaceTimeKroneckerElement
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanElement
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.local_view import NonStationaryLocal, ExpandedLocalHyperparameters, ExtendedCombinationHyperparameters
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

from eustace.analysis.observationsource import ObservationSource

from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eustace.outputformats import definitions

from inputloader_rawbinary import AnalysisSystemInputLoaderRawBinary_OneDay

class ClimatologySetup():
  """Define the setup for the climatology component"""
  
  def __init__(self):
    self.n_triangulation_divisions=6
    self.n_harmonics=4
    self.n_spatial_components=5
    self.amplitude=2.
    self.space_length_scale=5. #length scale in units of degrees
    self.grandmean_amplitude=15.0

class LargeScaleSetup():
  """Define the setup for the large scale component"""
  
  def __init__(self):
    
    # Note: node spacing is constant in seconds, not months.
    # Months are not the same length in seconds so nodes
    # are not exactly placed at month starts.
    time_extension = relativedelta.relativedelta(years = 1)
    model_startdate = datetime(2006, 1, 1)-time_extension
    model_enddate = datetime(2013, 1, 1)+time_extension
    avg_months_per_node = 3
    
    self.n_triangulation_divisions=3
    self.alpha=2
    self.starttime=eustace.timeutils.epoch.days_since_epoch( model_startdate )
    self.endtime=eustace.timeutils.epoch.days_since_epoch( model_enddate )
    self.n_nodes = (len( list(rrule(MONTHLY, dtstart=model_startdate, until=model_enddate)) ) - 1 ) / avg_months_per_node +1
    self.overlap_factor=2.5
    self.H = 1
    self.amplitude=1.
    self.space_length_scale=15.0  # length scale in units of degrees
    self.time_length_scale=15.0   # length scal in units of days
    self.bias_amplitude=.9



class LocalSetup():
  """Define the setup for the local component"""
  
  def __init__(self):
    self.n_triangulation_divisions=4#7
    self.amplitude=2.
    self.space_length_scale=2.0   #length scale in units of degrees
    self.bias_amplitude=15.0

class NonStationaryLocalSetup():
  """Define the setup for the non-stational local component - hyperparameter file specified in catalogue"""
  
  def __init__(self):
    self.n_triangulation_divisions=4#7
    self.bias_amplitude=15.0

class ClimatologyDefinition(ComponentStorage_InMemory):
    """Define climatology component."""

    def __init__(self, covariates_descriptor):
        if covariates_descriptor is not None:
            loader = LoadCovariateElement(covariates_descriptor)
            loader.check_keys()
            covariate_elements, covariate_hyperparameters = loader.load_covariates_and_hyperparameters()
            print('The following fields have been added as covariates of the climatology model')
            print(loader.data.keys())
        else:
            covariate_elements, covariate_hyperparameters = [], []

        setup = ClimatologySetup()
        
        super(ClimatologyDefinition, self).__init__(
            CombinationElement( [SeasonalElement(setup.n_triangulation_divisions, 
                                                 setup.n_harmonics, 
                                                 include_local_mean=True), GrandMeanElement()] + covariate_elements),
            CombinationHyperparameters( [SeasonalHyperparameters(setup.n_spatial_components, 
                                                                 numpy.log(setup.amplitude), 
                                                                 numpy.log(numpy.radians(setup.space_length_scale))), 
                                         CovariateHyperparameters(numpy.log(setup.grandmean_amplitude))] + covariate_hyperparameters))

class LargeScaleDefinition(ComponentStorage_InMemory):
    """Define large scale component."""

    def __init__(self, bias_terms = False, breakpoints_file = None):
      
	setup = LargeScaleSetup()
	
	if bias_terms:
	  bias_element = [InsituLandBiasElement(breakpoints_file, apply_policy=True, cut_value=3)]
	  bias_hyperparameters = [CovariateHyperparameters(numpy.log(setup.bias_amplitude))]
	else:
	  bias_element, bias_hyperparameters = [], []

        super(LargeScaleDefinition, self).__init__(
            CombinationElement([SpaceTimeKroneckerElement(setup.n_triangulation_divisions, 
                                                          setup.alpha, 
                                                          setup.starttime, 
                                                          setup.endtime, 
                                                          setup.n_nodes, 
                                                          setup.overlap_factor, 
                                                          setup.H)] + bias_element),
            CombinationHyperparameters([SpaceTimeSPDEHyperparameters(numpy.log(setup.amplitude), 
                                                                     numpy.log(numpy.radians(setup.space_length_scale)), 
                                                                     numpy.log(setup.time_length_scale))] + bias_hyperparameters))
        
class LocalDefinition(ComponentStorage_InMemory):
    """Define local component."""
    
    def __init__(self, bias_terms = False, global_biases_group_list = []):         
    
        setup = LocalSetup()
        
        if bias_terms:
            bias_elements = [BiasElement(groupname, 1) for groupname in global_biases_group_list]
            bias_hyperparameters = [CovariateHyperparameters(numpy.log(setup.bias_amplitude)) for index in range(len(global_biases_group_list))]
        else:
            bias_elements, bias_hyperparameters = [], []

        super(LocalDefinition, self).__init__(
            CombinationElement([LocalElement(setup.n_triangulation_divisions)] + bias_elements),         
            CombinationHyperparameters([LocalHyperparameters(numpy.log(setup.amplitude), 
                                                             numpy.log(numpy.radians(setup.space_length_scale)))] + bias_hyperparameters))

class NonStationaryLocalDefinition(ComponentStorage_InMemory):
    """Define non stational with stored hyperparameters local component."""
    
    def __init__(self, bias_terms = False, global_biases_group_list = [], local_hyperparameter_file = None):         
    
        setup = NonStationaryLocalSetup()
        
        if bias_terms:
            bias_elements = [BiasElement(groupname, 1) for groupname in global_biases_group_list]
            bias_hyperparameters = [CovariateHyperparameters(numpy.log(setup.bias_amplitude)) for index in range(len(global_biases_group_list))]
        else:
            bias_elements, bias_hyperparameters = [], []

        local_hyperparameters = ExpandedLocalHyperparameters(log_sigma = None, log_rho = None)
        local_hyperparameters.values_from_npy_savefile(local_hyperparameter_file)
        
        super(NonStationaryLocalDefinition, self).__init__(
            CombinationElement([NonStationaryLocal(setup.n_triangulation_divisions)] + bias_elements),         
            CombinationHyperparameters([local_hyperparameters]+bias_hyperparameters))

class AnalysisSystem_EUSTACE(AnalysisSystem):
    """Analysis system in which the storage objects for each component are specified during construction."""
    
    def __init__(self, storage_climatology, storage_large_scale, storage_local, covariates_descriptor, 
                        insitu_biases = False, breakpoints_file = None, global_biases = False, global_biases_group_list = [], 
                        compute_uncertainties = False, method = 'EXACT',
                        compute_sample = False, sample_size = definitions.GLOBAL_SAMPLE_SHAPE[3]):

        super(AnalysisSystem_EUSTACE, self).__init__(

            components = [ 
                SpaceTimeComponent(ClimatologyDefinition(covariates_descriptor), storage_climatology, True, compute_uncertainties, method, compute_sample, sample_size),
                SpaceTimeComponent(LargeScaleDefinition(insitu_biases, breakpoints_file), storage_large_scale, True, compute_uncertainties, method, compute_sample, sample_size),
                SpatialComponent(LocalDefinition(global_biases, global_biases_group_list), storage_local, compute_uncertainties, method, compute_sample, sample_size)],

            observable = ObservationSource.TMEAN)

def process_inputs(storage_climatology, storage_large_scale, storage_local, 
                    inputsources, time_index, 
                    component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                    compute_uncertainties, method,
                    compute_sample, sample_size):

    # Build inputloaders from list of sources
    inputloaders = [ AnalysisSystemInputLoaderRawBinary_OneDay(time_index=time_index, **source) for source in inputsources ]

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                                            covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                                            compute_uncertainties, method,
                                            compute_sample, sample_size)
    
    # Build and store measurement systems
    analysissystem.update_component_time(inputloaders, component_index, time_index)

def solve(storage_climatology, storage_large_scale, storage_local, 
            component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
            compute_uncertainties, method,
            compute_sample, sample_size):

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                                            covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                                            compute_uncertainties, method,
                                            compute_sample, sample_size)

    # Solve
    analysissystem.update_component_solution(component_index)

def measurement_merge(storage_climatology, storage_large_scale, storage_local, 
      component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
      compute_uncertainties, method, output_file):

    import pickle

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                        compute_uncertainties, method)

    # Get the merged update object
    merged_updates = analysissystem.components[component_index].component_solution().merge_updates()
    
    # Save to disk
    with open(output_file, 'wb') as f:
        pickle.dump(merged_updates, f)

def output_grid(storage_climatology, storage_large_scale, storage_local, 
                outputfile, processdate, time_index, 
                covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                compute_uncertainties, method,
                compute_sample, sample_size):

    print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                                            covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                                            compute_uncertainties, method,
                                        compute_sample, sample_size)

    #Configure output grid
    outputstructure = OutputRectilinearGridStructure(
	time_index, processdate,
	latitudes=numpy.linspace(-90.+definitions.GLOBAL_FIELD_RESOLUTION/2., 90.-definitions.GLOBAL_FIELD_RESOLUTION/2., num=definitions.GLOBAL_FIELD_SHAPE[1]),
	longitudes=numpy.linspace(-180.+definitions.GLOBAL_FIELD_RESOLUTION/2., 180.-definitions.GLOBAL_FIELD_RESOLUTION/2., num=definitions.GLOBAL_FIELD_SHAPE[2]))

    # Evaluate expected value at these locations
    for field in ['MAP', 'post_STD']:
      print 'Evaluating: ',field
      result_expected_value = analysissystem.evaluate_expected_value('MAP', outputstructure, 'GRID_CELL_AREA_AVERAGE', [1,1], 1000)
      result_expected_uncertainties = analysissystem.evaluate_expected_value('post_STD', outputstructure, 'GRID_CELL_AREA_AVERAGE', [1,1], 1000)
      
    print 'Evaluating: climatology fraction'
    climatology_fraction = analysissystem.evaluate_climatology_fraction(outputstructure, [1,1], 1000)
   
    print 'Evaluating: the sample'
    sample = analysissystem.evaluate_projected_sample(outputstructure)
   
    # Save results
    filebuilder = FileBuilderGlobalField(
        outputfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'Infilling Example',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Example data only',
        __name__, 
        '')
    filebuilder.add_global_field(definitions.TAS, result_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, result_expected_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_CLIMATOLOGY_FRACTION, climatology_fraction.reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    for index in range(definitions.GLOBAL_SAMPLE_SHAPE[3]):
      variable = copy.deepcopy(definitions.TASENSEMBLE)
      variable.name = variable.name + '_' + str(index)
      selected_sample = sample[:,index].ravel()+result_expected_value
      filebuilder.add_global_field(variable, selected_sample.reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    filebuilder.save_and_close()

def output_grid_component(storage_climatology, storage_large_scale, storage_local, 
        outputfile, processdate, time_index, 
        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
        compute_uncertainties, method):

    print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method)

    # Configure output grid
    outputstructure = OutputRectilinearGridStructure(
        time_index, processdate,
        latitudes=numpy.linspace(-89.875, 89.875, num=definitions.GLOBAL_FIELD_SHAPE[1]),
        longitudes=numpy.linspace(-179.875, 179.875, num=definitions.GLOBAL_FIELD_SHAPE[2]))

    # Evaluate expected value at these locations
    for field in ['MAP', 'post_STD']:
      print 'Evaluating: ',field
      result_expected_value = analysissystem.evaluate_expected_value('MAP', outputstructure, 'GRID_CELL_AREA_AVERAGE', [1,1], 1000)
      result_expected_uncertainties = analysissystem.evaluate_expected_value('post_STD', outputstructure, 'GRID_CELL_AREA_AVERAGE', [1,1], 1000)
      
    # Save results
    filebuilder = FileBuilderGlobalField(
        outputfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'Infilling Example',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Example data only',
        __name__, 
        '')
    filebuilder.add_global_field(definitions.TAS, result_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, result_expected_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    filebuilder.save_and_close()