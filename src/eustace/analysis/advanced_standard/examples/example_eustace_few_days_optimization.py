"""Example infilled analysis using a small amount of EUSTACE daily data."""

import argparse
import eustaceconfig
import numpy
import os.path
import time 

from datetime import datetime
from dateutil.relativedelta import relativedelta

from eustaceconfig import WORKSPACE_PATH

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.analysissystem import memory_usage_megabytes

from eustace.analysis.advanced_standard.optimizationsystem import OptimizationSystem

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.spatialdelayed import DelayedSpatialComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.elements.bias_insitu_land import InsituLandBiasElement
from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import LoadCovariateElement
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanElement
from eustace.analysis.advanced_standard.elements.kronecker import SpaceTimeKroneckerElement
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.local_view import LocalSubRegion, LocalSuperTriangle, NonStationaryLocal
from eustace.analysis.advanced_standard.elements.local_view import ExtendedLocalHyperparameters, ExtendedCombinationHyperparameters, ExpandedLocalHyperparameters
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters

from eustace.analysis.advanced_standard.examples.example_optimization import split_states_time, extract_local_view_states_time

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

from eustace.analysis.observationsource import ObservationSource

from eustace.outputformats import definitions

from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField

from eustace.timeutils.epoch import days_since_epoch

from eustaceconfig import WORKSPACE_PATH

from inputloader_rawbinary import AnalysisSystemInputLoaderRawBinary_Sources

BREAKPOINTS_FILE = os.path.join(eustaceconfig.WORKSPACE_PATH,'data/internal/D1.7/daily/eustace_stations_global_R001127_daily_status.nc')
GLOBAL_BIASES_GROUP_LIST = ['surfaceairmodel_ice_global' , 'surfaceairmodel_land_global', 'surfaceairmodel_ocean_global']

def main():

    print 'Advanced standard example using a few days of EUSTACE data'
    parser = argparse.ArgumentParser(description='Advanced standard example using a few days of EUSTACE data')
    parser.add_argument('outpath', help='directory where the output should be redirected')
    parser.add_argument('--json_descriptor', default = None, help='a json descriptor containing the covariates to include in the climatology model')
    parser.add_argument('--land_biases', action='store_true', help='include insitu land homogenization bias terms')
    parser.add_argument('--global_biases', action='store_true', help='include global satellite bias terms')
    parser.add_argument('--n_iterations', type=int, default=5, help='number of solving iterations')
    args = parser.parse_args()

    # Input data path
    basepath = os.path.join('/work/scratch/eustace/rawbinary3')

    # Days to process
    #time_indices = range(int(days_since_epoch(datetime(2006, 2, 1))), int(days_since_epoch(datetime(2006, 2, 2))))
    #time_indices = range(int(days_since_epoch(datetime(1906, 2, 1))), int(days_since_epoch(datetime(1906, 2, 2))))
    
    date_list = [datetime(2006, 1, 1) + relativedelta(days=k) for k in range(3)]
    
    #backwards_list = [date_list[i] for i in range(11, -1, -1)]
    #date_list = backwards_list
    
    time_indices = [int(days_since_epoch(date)) for date in date_list]
    
    # Sources to use
    sources = [ 'surfaceairmodel_land', 'surfaceairmodel_ocean', 'surfaceairmodel_ice', 'insitu_land', 'insitu_ocean' ] 
    sources = [ 'insitu_land', 'insitu_ocean' ]
    #sources = [ 'surfaceairmodel_land' ]
    # CLIMATOLOGY COMPONENT: combining the seasonal core along with latitude harmonics, altitude and coastal effects
    
    if args.json_descriptor is not None:
      loader = LoadCovariateElement(args.json_descriptor)
      loader.check_keys()
      covariate_elements, covariate_hyperparameters = loader.load_covariates_and_hyperparameters()
      print('The following fields have been added as covariates of the climatology model')
      print(loader.data.keys())
    else:
      covariate_elements, covariate_hyperparameters = [], []

    #climatology_element = CombinationElement( [SeasonalElement(n_triangulation_divisions=2, n_harmonics=2, include_local_mean=False), GrandMeanElement()]+covariate_elements)
    #climatology_hyperparameters = CombinationHyperparameters( [SeasonalHyperparameters(n_spatial_components=2, common_log_sigma=0.0, common_log_rho=0.0), CovariateHyperparameters(numpy.log(15.0))] + covariate_hyperparameters )
    
    climatology_element = CombinationElement( [GrandMeanElement(),]+covariate_elements)
    climatology_hyperparameters = CombinationHyperparameters( [CovariateHyperparameters(numpy.log(15.0)),] + covariate_hyperparameters )
    
    #climatology_element =SeasonalElement(n_triangulation_divisions=2, n_harmonics=2, include_local_mean=False)
    #climatology_hyperparameters = SeasonalHyperparameters(n_spatial_components=2, common_log_sigma=0.0, common_log_rho=0.0)
    
    climatology_component = SpaceTimeComponent(ComponentStorage_InMemory(climatology_element, climatology_hyperparameters), SpaceTimeComponentSolutionStorage_InMemory(), compute_uncertainties=True, method='APPROXIMATED')


    # LARGE SCALE (kronecker product) COMPONENT: combining large scale trends with bias terms accounting for homogeneization effects

    if args.land_biases:
        bias_element, bias_hyperparameters = [InsituLandBiasElement(BREAKPOINTS_FILE)], [CovariateHyperparameters(numpy.log(.9))]
        print('Adding bias terms for insitu land homogenization')
    else:
        bias_element, bias_hyperparameters = [], []
        
    large_scale_element = CombinationElement( [SpaceTimeKroneckerElement(n_triangulation_divisions=2, alpha=2, starttime=-30, endtime=365*1+30, n_nodes=12*1+2 , overlap_factor=2.5, H=1)] + bias_element)
    large_scale_hyperparameters = CombinationHyperparameters( [SpaceTimeSPDEHyperparameters(space_log_sigma=0.0, space_log_rho=numpy.log(numpy.radians(15.0)), time_log_rho=numpy.log(15.0))] + bias_hyperparameters) 
    large_scale_component =  SpaceTimeComponent(ComponentStorage_InMemory(large_scale_element, large_scale_hyperparameters), SpaceTimeComponentSolutionStorage_InMemory(), compute_uncertainties=True, method='APPROXIMATED')
                             
    
    # LOCAL COMPONENT: combining local scale variations with global satellite bias terms

    if args.global_biases:
        bias_elements = [BiasElement(groupname, 1) for groupname in GLOBAL_BIASES_GROUP_LIST]
        bias_hyperparameters = [CovariateHyperparameters(numpy.log(15.0)) for index in range(3)]
        print('Adding global bias terms for all the surfaces')
    else:
        bias_elements, bias_hyperparameters = [], []

    n_triangulation_divisions_local = 7
    local_log_sigma = numpy.log(5)
    local_log_rho   = numpy.log(numpy.radians(5.0))
    local_element = NonStationaryLocal(n_triangulation_divisions=n_triangulation_divisions_local)
    n_local_nodes = local_element.spde.n_latent_variables()
    local_scale_element = CombinationElement([local_element] + bias_elements)
    local_hyperparameters = ExpandedLocalHyperparameters(log_sigma=numpy.repeat(local_log_sigma, n_local_nodes), log_rho=numpy.repeat(local_log_rho, n_local_nodes))
    local_scale_hyperparameters = CombinationHyperparameters([local_hyperparameters] + bias_hyperparameters)
    local_component = DelayedSpatialComponent(ComponentStorage_InMemory(local_scale_element, local_scale_hyperparameters), SpatialComponentSolutionStorage_InMemory(), compute_uncertainties=True, method='APPROXIMATED')
    print "hyperparameter storage:", local_component.storage.hyperparameters
    print 'Analysing inputs'

    # Analysis system using the specified components, for the Tmean observable
    ##analysis_system = AnalysisSystem(
    ##    [ climatology_component, large_scale_component, local_component ],
    ##    ObservationSource.TMEAN)

    analysis_system = OptimizationSystem([ climatology_component, local_component], 
                                      ObservationSource.TMEAN)


    # Object to load raw binary inputs at time indices
    inputloaders = [ AnalysisSystemInputLoaderRawBinary_Sources(basepath, source, time_indices) for source in sources ]

    for iteration in range(args.n_iterations):
	
	message = 'Iteration {}'.format(iteration)
	print(message)
	
	# Update with data
	analysis_system.update(inputloaders, time_indices)
    
    
    
    ##################################################
    
    # Optimize local model hyperparameters
    
    # Loop over local regions, generate optimization systems, fit hyperparameters and save
    
    # split spde and bias models for local component into two components
    global_spde_sub_component_definition = ComponentStorage_InMemory(CombinationElement([local_element]), CombinationHyperparameters([local_hyperparameters]))
    global_spde_sub_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
    global_spde_sub_component = DelayedSpatialComponent(global_spde_sub_component_definition, global_spde_sub_component_storage_solution)
    
    bias_sub_component_definition = ComponentStorage_InMemory(CombinationElement(bias_elements), CombinationHyperparameters(bias_hyperparameters))
    bias_sub_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
    bias_sub_component = DelayedSpatialComponent(bias_sub_component_definition, bias_sub_component_storage_solution)
    
    element_optimisation_flags = [True, False, False, False] # one spde, three biases
    
    for time_key in time_indices:
        split_states_time( local_component, global_spde_sub_component, bias_sub_component, element_optimisation_flags, time_key )
    
    # Define subregions and extract their states
    neighbourhood_level = 1
        
    n_subregions = global_spde_sub_component.storage.element_read().combination[0].spde.n_triangles_at_level(neighbourhood_level)
    hyperparameter_file_template = "local_hyperparameters.%i.%i.%i.npy"
    
    
    
    fit_hyperparameters = True
    optimization_component_index = 2
    if fit_hyperparameters:
        for region_index in range(n_subregions):
            # Setup model for local subregion of neighours with super triangle
            view_flags = [True,]
            region_element = CombinationElement([LocalSubRegion(n_triangulation_divisions_local, neighbourhood_level, region_index)])
            region_hyperparameters = ExtendedCombinationHyperparameters( [LocalHyperparameters( log_sigma = local_log_sigma, log_rho = local_log_rho )] )
            region_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
            region_sub_component = DelayedSpatialComponent(ComponentStorage_InMemory(region_element, region_hyperparameters), region_component_storage_solution)
            
            for time_key in time_indices:
                print "region_index, time_key:", region_index, time_key
                extract_local_view_states_time( global_spde_sub_component, region_sub_component, view_flags, time_key )
        
            print "running optimization for region:", region_index

            region_optimization_system = OptimizationSystem([ climatology_component, bias_sub_component, region_sub_component], 
                                                        ObservationSource.TMEAN)
        
            for time_key in time_indices:                                                
                region_optimization_system.update_component_time(inputloaders, optimization_component_index, time_key)
        
            # commented version that works for few days inputs
            #region_optimization_system.components[optimization_component_index].component_solution().optimize()
            #region_optimization_system.components[optimization_component_index].storage.hyperparameters.get_array()
            #hyperparameter_file = os.path.join(args.outpath, hyperparameter_file_template % (n_triangulation_divisions_local, neighbourhood_level, region_index) )
            #region_sub_component.storage.hyperparameters.values_to_npy_savefile( hyperparameter_file )
            
            # replaced with version for full processing based json dump of input files - need to generate the input_descriptor dict
            hyperparameter_file = os.path.join(args.outpath, hyperparameter_file_template % (n_triangulation_divisions_local, neighbourhood_level, region_index) )
            region_optimization_system.process_inputs(input_descriptor, optimization_component_index, time_indices)
            region_optimization_system.optimize_component(optimization_component_index, hyperparameter_storage_file = hyperparameter_file)
            
            fitted_hyperparameters_converted = region_sub_component.storage.hyperparameters.get_array()
            fitted_hyperparameters_converted[0] = numpy.exp( fitted_hyperparameters_converted[0] )
            fitted_hyperparameters_converted[1] = numpy.exp( fitted_hyperparameters_converted[1] ) * 180.0/ numpy.pi
            print 'fitted_hyperparameters_converted:', fitted_hyperparameters_converted
        
    # Setup model for the super triangle without neighbours for hyperparameter merging
    region_spdes = []
    region_hyperparameter_values = []
    for region_index in range(n_subregions):
        # Redefine the region sub component as a supertriangle rather than a neighbourhood
        region_element = CombinationElement([LocalSuperTriangle(n_triangulation_divisions_local, neighbourhood_level, region_index)])
        region_hyperparameters = ExtendedCombinationHyperparameters( [LocalHyperparameters( log_sigma = local_log_sigma, log_rho = local_log_rho )] )
        region_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
        region_sub_component = DelayedSpatialComponent(ComponentStorage_InMemory(region_element, region_hyperparameters), region_component_storage_solution)
        
        # Read the optimized hyperparameters
        hyperparameter_file = os.path.join(args.outpath, hyperparameter_file_template % (n_triangulation_divisions_local, neighbourhood_level, region_index) )
        region_sub_component.storage.hyperparameters.values_from_npy_savefile( hyperparameter_file )

        # Append the spde model and hyperparameters to their lists for merging
        region_spdes.append(region_element.combination[0].spde)
        region_hyperparameter_values.append(region_sub_component.storage.hyperparameters.get_array())

            
    # merge and save hyperparameters
    full_spde = local_element.spde
    new_hyperparameter_values, global_sigma_design, global_rho_design = full_spde.merge_local_parameterisations( region_spdes, region_hyperparameter_values, merge_method = 'exp_average' )
    
    local_hyperparameters.set_array(new_hyperparameter_values)
    hyperparameter_file_merged = "merged_hyperparameters.%i.%i.npy" % (n_triangulation_divisions_local, neighbourhood_level)    
    local_hyperparameters.values_to_npy_savefile(  os.path.join(args.outpath, hyperparameter_file_merged) )
    
    # Refit local model with the optimized hyperparameters    
    analysis_system.update_component(inputloaders, 1, time_indices)
    
    ##################################################
    
    print 'Computing outputs'
    
    # Produce an output for each time index
    for time_index in time_indices:

        # Get date for output
        outputdate = inputloaders[0].datetime_at_time_index(time_index)
        print 'Evaluating output grid: ', outputdate

        #Configure output grid
        outputstructure = OutputRectilinearGridStructure(
            time_index, outputdate,
            latitudes=numpy.linspace(-89.875, 89.875, num=definitions.GLOBAL_FIELD_SHAPE[1]),
            longitudes=numpy.linspace(-179.875, 179.875, num=definitions.GLOBAL_FIELD_SHAPE[2]))

        # print 'Size of grid : ', outputstructure.number_of_observations()

        # Evaluate expected value at these locations
        result_expected_value = analysis_system.evaluate_expected_value('MAP', outputstructure, 'POINTWISE')
	result_expected_uncertainties = analysis_system.evaluate_expected_value('post_STD', outputstructure, 'POINTWISE')

	# Make output filename
        pathname = 'eustace_example_output_{0:04d}{1:02d}{2:02d}.nc'.format(outputdate.year, outputdate.month, outputdate.day)
	pathname = os.path.join(args.outpath, pathname)
        print 'Saving: ', pathname

        # Save results
        filebuilder = FileBuilderGlobalField(
            pathname, 
            time_index,
            'Infilling Example',
            'UNVERSIONED',
            definitions.TAS.name,
            '',
            'Example data only',
            'eustace.analysis.advanced_standard.examples.example_eustace_few_days', 
            '')
        filebuilder.add_global_field(definitions.TAS, result_expected_value.reshape(definitions.GLOBAL_FIELD_SHAPE))
        filebuilder.add_global_field(definitions.TASUNCERTAINTY, result_expected_uncertainties.reshape(definitions.GLOBAL_FIELD_SHAPE))
	filebuilder.save_and_close()

    print 'Complete'

if __name__=='__main__':

    # Call main method
    main()
