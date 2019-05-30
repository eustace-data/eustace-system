
import eumopps.version.svn
import eumopps.catalogue.catalogue
import eumopps.catalogue.dataset
import eumopps.catalogue.operation
import eumopps.catalogue.placeholder
import eumopps.catalogue.step
from stepmidmonth import StepMidMonth
from eumopps.catalogue.fileio.formatnetcdf import CatalogueWriterNetCDF, CatalogueReaderNetCDF
import numpy
import os.path
import iris

from eustaceconfig import WORKSPACE_PATH
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters
from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeFactorElement
from eustace.analysis.advanced_standard.elements.spacetimespde import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure
from inputloader_hadcrut4 import AnalysisSystemInputLoaderHadCRUT4
from outputformat_hadcrut4 import FileBuilderHadCRUT4ExampleOutput
from outputformat_hadcrut4 import TAS_ANOMALY

class ClimatologyDefinition(ComponentStorage_InMemory):
    """Define climatology component."""

    def __init__(self):

        super(ClimatologyDefinition, self).__init__(
            SeasonalElement(n_triangulation_divisions=5, n_harmonics=5, include_local_mean=True),
            SeasonalHyperparameters(n_spatial_components=6, common_log_sigma=1.0, common_log_rho=0.0))

class LargeScaleDefinition(ComponentStorage_InMemory):
    """Define large scale component."""

    def __init__(self):

        # Number of factors for large scale (factor analysis) component and initial hyperparameters
        n_factors = 5
        factors = [ ]
        factor_hyperparameters = [ ]
        for factor_index in range(n_factors):

            factor_hyperparameters.append( SpaceTimeSPDEHyperparameters(
                    space_log_sigma=0.0,
                    space_log_rho=numpy.log(10.0 * numpy.pi/180 + 25.0 * numpy.pi/180 *(n_factors - factor_index) / n_factors),
                    time_log_rho=numpy.log(1/12.0 + 6/12.0*(n_factors - factor_index) / n_factors)) )
            
            factors.append( SpaceTimeFactorElement(n_triangulation_divisions=5, alpha=2, starttime=0, endtime=36, overlap_factor=2.5, H=1) )
            
        super(LargeScaleDefinition, self).__init__(
            CombinationElement(factors),
            CombinationHyperparameters(factor_hyperparameters))
        
class LocalDefinition(ComponentStorage_InMemory):
    """Define local component."""
    
    def __init__(self):                      
 
        super(LocalDefinition, self).__init__(
            LocalElement(n_triangulation_divisions=4),
            LocalHyperparameters(log_sigma=0.0, log_rho=numpy.log(10.0 * numpy.pi/180)))

class HadCRUT4_InMemory(object):
    """Class for purpose of holding the global in-memory variables."""

    definition_climatology = ClimatologyDefinition()
    definition_large_scale = LargeScaleDefinition()
    definition_local = LocalDefinition()
    solution_climatology = SpaceTimeComponentSolutionStorage_InMemory()
    solution_large_scale = SpaceTimeComponentSolutionStorage_InMemory()
    solution_local = SpatialComponentSolutionStorage_InMemory()

class AnalysisSystem_HadCRUT4_InMemory(AnalysisSystem):

    def __init__(self):

        # Build analysis system by referencing the global variables
        # In distributed computing we would construct these using classes that load from file

        super(AnalysisSystem_HadCRUT4_InMemory, self).__init__(

            components = [ 
                SpaceTimeComponent(HadCRUT4_InMemory.definition_climatology, HadCRUT4_InMemory.solution_climatology),
                SpaceTimeComponent(HadCRUT4_InMemory.definition_large_scale, HadCRUT4_InMemory.solution_large_scale),
                SpatialComponent(HadCRUT4_InMemory.definition_local, HadCRUT4_InMemory.solution_local) ],

            observable = ObservationSource.TMEAN)

def process_month(inputfilelist, time_index, component_index, processdate):
    
    inputloaders = [ AnalysisSystemInputLoaderHadCRUT4(inputfilelist) ]
    AnalysisSystem_HadCRUT4_InMemory().update_component_time(inputloaders, component_index, time_index)

def solve(component_index):

    AnalysisSystem_HadCRUT4_InMemory().update_component_solution(component_index)

def output_month(outputfile, time_index, processdate):

    print 'Saving: ', processdate
    print 'Output: ', outputfile

    # Configure output grid
    outputstructure = OutputRectilinearGridStructure(
        time_index, processdate,
        latitudes=numpy.linspace(-87.5, 87.5, num=36),
        longitudes=numpy.linspace(-177.5, 177.5, num=72))

    # Evaluate expected value at these locations
    result_expected_value = AnalysisSystem_HadCRUT4_InMemory().evaluate_expected_value(outputstructure)

    # Save results
    filebuilder = FileBuilderHadCRUT4ExampleOutput(outputfile, outputstructure)
    filebuilder.add_global_field(TAS_ANOMALY, result_expected_value.reshape(1,36,72))
    filebuilder.save_and_close()

def main():

    print 'EUSTACE example using HadCRUT4 monthly data'

    eumopps.version.svn.set_allow_unversioned_code(True)

    # Catalogue to use
    cataloguefilename = 'catalogue.nc'

    # Output path assumed current directory
    outputpath = '.'

    # Input data path
    inputpath = os.path.join(WORKSPACE_PATH, 'data/incoming/HadCRUT4.5.0.0')

    # Number of months to process
    nummonths = 1

    # Build catalogue of inputs
    catalogue = eumopps.catalogue.catalogue.Catalogue(

        datasets = [

            # Input data set
            eumopps.catalogue.dataset.CatalogueDataSet(
                name='hadcrut4',
                path=inputpath,
                subsets = [
                    eumopps.catalogue.dataset.CatalogueDataSubset(layout=eumopps.catalogue.storage.DataStorageFiles(patterns=[ 'hadcrut4_median_netcdf.nc' ])),
                    eumopps.catalogue.dataset.CatalogueDataSubset(layout=eumopps.catalogue.storage.DataStorageFiles(patterns=[ 'hadcrut4_uncorrelated_supplementary.nc' ])),
                    eumopps.catalogue.dataset.CatalogueDataSubset(layout=eumopps.catalogue.storage.DataStorageFiles(patterns=[ 'hadcrut4_blended_uncorrelated.nc' ]))
                    ]),

            # The data set to make
            eumopps.catalogue.dataset.CatalogueDataSet(
                name='infilled',
                path=outputpath,
                subsets = [
                    eumopps.catalogue.dataset.CatalogueDataSubset(layout=eumopps.catalogue.storage.DataStorageFiles(patterns=[ 'example_hadcrut4_%Y%m.nc' ])),
                    ])
            ],

        operations = [

            eumopps.catalogue.operation.Operation(
                name = 'climatology_component_input',
                runmodule = {
                    'python_function': 'eustace.analysis.advanced_standard.examples.example_hadcrut4_eumopps.process_month',
                    'inputfilelist':  [ 
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 0),
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 1),
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 2) 
                        ],
                    'time_index': eumopps.catalogue.placeholder.StepIndex(),
                    'component_index': 0,
                    'processdate': eumopps.catalogue.placeholder.StepTime(),
                    },
                step = StepMidMonth('18500101000000', nummonths)),

            eumopps.catalogue.operation.Operation(
                name = 'climatology_component_solve',
                runmodule = {
                    'python_function': 'eustace.analysis.advanced_standard.examples.example_hadcrut4_eumopps.solve',
                    'component_index': 0
                    },
                step = eumopps.catalogue.step.StepOnce()),

            eumopps.catalogue.operation.Operation(
                name = 'large_scale_component_input',
                runmodule = {
                    'python_function': 'eustace.analysis.advanced_standard.examples.example_hadcrut4_eumopps.process_month',
                    'inputfilelist':  [ 
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 0),
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 1),
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 2) 
                        ],
                    'time_index': eumopps.catalogue.placeholder.StepIndex(),
                    'component_index': 1,
                    'processdate': eumopps.catalogue.placeholder.StepTime(),
                    },
                step = StepMidMonth('18500101000000', nummonths)),

            eumopps.catalogue.operation.Operation(
                name = 'large_scale_component_solve',
                runmodule = {
                    'python_function': 'eustace.analysis.advanced_standard.examples.example_hadcrut4_eumopps.solve',
                    'component_index': 1
                    },
                step = eumopps.catalogue.step.StepOnce()),

            eumopps.catalogue.operation.Operation(
                name = 'local_component_input',
                runmodule = {
                    'python_function': 'eustace.analysis.advanced_standard.examples.example_hadcrut4_eumopps.process_month',
                    'inputfilelist':  [ 
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 0),
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 1),
                        eumopps.catalogue.placeholder.InputFile('hadcrut4', 2) 
                        ],
                    'time_index': eumopps.catalogue.placeholder.StepIndex(),
                    'component_index': 2,
                    'processdate': eumopps.catalogue.placeholder.StepTime(),
                    },
                step = StepMidMonth('18500101000000', nummonths)),

            eumopps.catalogue.operation.Operation(
                name = 'local_component_solve',
                runmodule = {
                    'python_function': 'eustace.analysis.advanced_standard.examples.example_hadcrut4_eumopps.solve',
                    'component_index': 2
                    },
                step = eumopps.catalogue.step.StepOnce()),

            eumopps.catalogue.operation.Operation(
                name = 'output_grid',
                runmodule = {
                    'python_function': 'eustace.analysis.advanced_standard.examples.example_hadcrut4_eumopps.output_month',
                    'outputfile':  eumopps.catalogue.placeholder.OutputFile('infilled'),
                    'time_index': eumopps.catalogue.placeholder.StepIndex(),
                    'processdate': eumopps.catalogue.placeholder.StepTime(),
                    },
                step = StepMidMonth('18500101000000', nummonths))
            ])
         
    print 'Analysing inputs'

    # Populate inputs
    catalogue.datasets[0].search()
    catalogue.datasets[0].checksum()
    print 'Found input files: ', [ str(subset.matches[0].name) for subset in catalogue.datasets[0].subsets ]

    # Resolve operation references
    for operation in catalogue.operations:
        operation.resolve_operation_references(catalogue, outputpath)

    # Store catalogue
    CatalogueWriterNetCDF().save(cataloguefilename, catalogue)

    # Do operations in catalogue
    for operationgroup in [
        'climatology_component_input',
        'climatology_component_solve',
        'large_scale_component_input',
        'large_scale_component_solve',
        'local_component_input',
        'local_component_solve',
        'output_grid'  ]:

        operation_count = CatalogueReaderNetCDF().operationcount(cataloguefilename, operationgroup)
        print 'Processing: {0}   Operations to run: {1}'.format(operationgroup, operation_count)

        # Load and run input-processing operations
        for operation_index in range(operation_count):

            operation = CatalogueReaderNetCDF().load_operation(cataloguefilename, operationgroup)
            operation.resolve_single_operation(cataloguefilename, operation_index, 'now')
            operation.run()
        
    # Catalogue final checksums and store
    catalogue.dataset('infilled').checksum()
    CatalogueWriterNetCDF().save(cataloguefilename, catalogue)
    
if __name__=="__main__":
    main()
