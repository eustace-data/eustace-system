"""Example infilled analysis using a small amount of EUSTACE daily data."""

import argparse
import copy
import eustaceconfig
import numpy
import os.path
import time 

from datetime import datetime

from eustaceconfig import WORKSPACE_PATH

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.analysissystem import memory_usage_megabytes

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
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanElement
from eustace.analysis.advanced_standard.elements.kronecker import SpaceTimeKroneckerElement
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters

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
    time_indices = range(int(days_since_epoch(datetime(2006, 2, 1))), int(days_since_epoch(datetime(2006, 2, 2))))

    # Sources to use
    sources = [ 'surfaceairmodel_land', 'surfaceairmodel_ocean', 'surfaceairmodel_ice', 'insitu_land', 'insitu_ocean' ]    

    #SETUP
    # setup for the seasonal core: climatology covariates setup read from file
    seasonal_setup = {'n_triangulation_divisions':5,
		      'n_harmonics':4,
		      'n_spatial_components':6,
		      'amplitude':2.,
		      'space_length_scale':5., # length scale in units of degrees
		     }
    grandmean_amplitude = 15.0
    
    # setup for the large scale component
    spacetime_setup = {'n_triangulation_divisions':2,
		       'alpha':2,
		       'starttime':0,
		       'endtime':10.,
		       'n_nodes':2,
		       'overlap_factor':2.5,
		       'H':1,
		       'amplitude':1.,
		       'space_lenght_scale':15.0, # length scale in units of degrees
		       'time_length_scale':15.0   # length scal in units of days
		      }
    bias_amplitude = .9

    # setup for the local component
    local_setup = {'n_triangulation_divisions':6,
                   'amplitude':2.,
                   'space_length_scale':2. # length scale in units of degrees
                  }
    globalbias_amplitude = 15.0

    # CLIMATOLOGY COMPONENT: combining the seasonal core along with latitude harmonics, altitude and coastal effects    
    if args.json_descriptor is not None:
      loader = LoadCovariateElement(args.json_descriptor)
      loader.check_keys()
      covariate_elements, covariate_hyperparameters = loader.load_covariates_and_hyperparameters()
      print('The following fields have been added as covariates of the climatology model')
      print(loader.data.keys())
    else:
      covariate_elements, covariate_hyperparameters = [], []

    climatology_element = CombinationElement( [SeasonalElement(n_triangulation_divisions=seasonal_setup['n_triangulation_divisions'], 
							       n_harmonics=seasonal_setup['n_harmonics'], 
							       include_local_mean=True), 
					       GrandMeanElement()]+covariate_elements)       
    climatology_hyperparameters = CombinationHyperparameters( [SeasonalHyperparameters(n_spatial_components=seasonal_setup['n_spatial_components'], 
										       common_log_sigma=numpy.log(seasonal_setup['amplitude']), 
										       common_log_rho=numpy.log(numpy.radians(seasonal_setup['space_length_scale']))), 
							       CovariateHyperparameters(numpy.log(grandmean_amplitude))] + covariate_hyperparameters )
    climatology_component = SpaceTimeComponent(ComponentStorage_InMemory(climatology_element, climatology_hyperparameters), SpaceTimeComponentSolutionStorage_InMemory(), 
                                                                         compute_uncertainties=True, method='APPROXIMATED',
                                                                         compute_sample=True, sample_size=definitions.GLOBAL_SAMPLE_SHAPE[3])

    # LARGE SCALE (kronecker product) COMPONENT: combining large scale trends with bias terms accounting for homogeneization effects    
    if args.land_biases:
	bias_element, bias_hyperparameters = [InsituLandBiasElement(BREAKPOINTS_FILE)], [CovariateHyperparameters(numpy.log(bias_amplitude))]
	print('Adding bias terms for insitu land homogenization')
    else:
	bias_element, bias_hyperparameters = [], []

    large_scale_element = CombinationElement( [SpaceTimeKroneckerElement(n_triangulation_divisions=spacetime_setup['n_triangulation_divisions'], 
                                                                         alpha=spacetime_setup['alpha'], 
                                                                         starttime=spacetime_setup['starttime'], 
                                                                         endtime=spacetime_setup['endtime'], 
                                                                         n_nodes=spacetime_setup['n_nodes'], 
                                                                         overlap_factor=spacetime_setup['overlap_factor'], 
                                                                         H=spacetime_setup['H'])] + bias_element)
    large_scale_hyperparameters = CombinationHyperparameters( [SpaceTimeSPDEHyperparameters(space_log_sigma=numpy.log(spacetime_setup['amplitude']),
                                                                                            space_log_rho=numpy.log(numpy.radians(spacetime_setup['space_lenght_scale'])), 
                                                                                            time_log_rho=numpy.log(spacetime_setup['time_length_scale']))] + bias_hyperparameters) 
    large_scale_component =  SpaceTimeComponent(ComponentStorage_InMemory(large_scale_element, large_scale_hyperparameters), SpaceTimeComponentSolutionStorage_InMemory(), 
                                                                          compute_uncertainties=True, method='APPROXIMATED',
                                                                          compute_sample=True, sample_size=definitions.GLOBAL_SAMPLE_SHAPE[3])
                                 
    # LOCAL COMPONENT: combining local scale variations with global satellite bias terms    
    if args.global_biases:
	bias_elements = [BiasElement(groupname, 1) for groupname in GLOBAL_BIASES_GROUP_LIST]
	bias_hyperparameters = [CovariateHyperparameters(numpy.log(globalbias_amplitude)) for index in range(len(GLOBAL_BIASES_GROUP_LIST))]
	print('Adding global bias terms for all the surfaces')
    else:
	bias_elements, bias_hyperparameters = [], []

    local_scale_element = CombinationElement([LocalElement(n_triangulation_divisions=local_setup['n_triangulation_divisions'])] + bias_elements)
    local_scale_hyperparameters = CombinationHyperparameters([LocalHyperparameters(log_sigma=numpy.log(local_setup['amplitude']), 
                                                                                   log_rho=numpy.log(numpy.radians(local_setup['space_length_scale'])))] + bias_hyperparameters)
    local_component = SpatialComponent(ComponentStorage_InMemory(local_scale_element, local_scale_hyperparameters), SpatialComponentSolutionStorage_InMemory(), 
                                                                 compute_uncertainties=True, method='APPROXIMATED',
                                                                 compute_sample=True, sample_size=definitions.GLOBAL_SAMPLE_SHAPE[3])

    # Analysis system using the specified components, for the Tmean observable
    print 'Analysing inputs'

    analysis_system = AnalysisSystem(
        [ climatology_component, large_scale_component, local_component ],
        ObservationSource.TMEAN)

    # Object to load raw binary inputs at time indices
    inputloaders = [ AnalysisSystemInputLoaderRawBinary_Sources(basepath, source, time_indices) for source in sources ]

    for iteration in range(args.n_iterations):
	
	message = 'Iteration {}'.format(iteration)
	print(message)
	
	# Update with data
	analysis_system.update(inputloaders, time_indices)

    print 'Computing outputs'

    # Produce an output for each time index
    for time_index in time_indices:

        # Get date for output
        outputdate = inputloaders[0].datetime_at_time_index(time_index)
        print 'Evaluating output grid: ', outputdate

        #Configure output grid
        outputstructure = OutputRectilinearGridStructure(
            time_index, outputdate,
            latitudes=numpy.linspace(-90.+definitions.GLOBAL_FIELD_RESOLUTION/2., 90.-definitions.GLOBAL_FIELD_RESOLUTION/2., num=definitions.GLOBAL_FIELD_SHAPE[1]),
            longitudes=numpy.linspace(-180.+definitions.GLOBAL_FIELD_RESOLUTION/2., 180.-definitions.GLOBAL_FIELD_RESOLUTION/2., num=definitions.GLOBAL_FIELD_SHAPE[2]))

        # Evaluate expected value at these locations
        for field in ['MAP', 'post_STD']:
	  print 'Evaluating: ',field
	  result_expected_value = analysis_system.evaluate_expected_value('MAP', outputstructure, 'GRID_CELL_AREA_AVERAGE', [1,1], 1000)
	  result_expected_uncertainties = analysis_system.evaluate_expected_value('post_STD', outputstructure, 'GRID_CELL_AREA_AVERAGE', [1,1], 1000)
	  
	print 'Evaluating: climatology fraction'
	climatology_fraction = analysis_system.evaluate_climatology_fraction(outputstructure, [1,1], 1000)

	print 'Evaluating: the sample'
	sample = analysis_system.evaluate_projected_sample(outputstructure)

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
        filebuilder.add_global_field(definitions.TAS_CLIMATOLOGY_FRACTION, climatology_fraction.reshape(definitions.GLOBAL_FIELD_SHAPE))

	for index in range(definitions.GLOBAL_SAMPLE_SHAPE[3]):
	  variable = copy.deepcopy(definitions.TASENSEMBLE)
	  variable.name = variable.name + '_' + str(index)
	  selected_sample = sample[:,index].ravel()+result_expected_value
	  filebuilder.add_global_field(variable, selected_sample.reshape(definitions.GLOBAL_FIELD_SHAPE))
	  
	filebuilder.save_and_close()

    print 'Complete'

if __name__=='__main__':

    # Call main method
    main()
