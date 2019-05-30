"""Example infilled analysis using HadCRUT4 input data."""

import numpy
import os.path
from eustaceconfig import WORKSPACE_PATH
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
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
from eustace.analysis.observationsource import ObservationSource
from inputloader_hadcrut4 import AnalysisSystemInputLoaderHadCRUT4
from outputformat_hadcrut4 import FileBuilderHadCRUT4ExampleOutput
from outputformat_hadcrut4 import TAS_ANOMALY
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

def main():

    print 'EUSTACE example using HadCRUT4 monthly data'

    # Input data path
    input_basepath = os.path.join(WORKSPACE_PATH, 'data/incoming/HadCRUT4.5.0.0')

    # Input filenames
    input_filenames = [
        'hadcrut4_median_netcdf.nc',
        'hadcrut4_uncorrelated_supplementary.nc',
        'hadcrut4_blended_uncorrelated.nc' ]

    # Months to process
    time_indices = range(2)

    # Climatology component
    climatology_component = SpaceTimeComponent(ComponentStorage_InMemory(SeasonalElement(n_triangulation_divisions=5, n_harmonics=5, include_local_mean=True),
                                                                         SeasonalHyperparameters(n_spatial_components=6, common_log_sigma=1.0, common_log_rho=0.0)),
                                               SpaceTimeComponentSolutionStorage_InMemory())

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

    # Large scale (factor analysis) component
    large_scale_component = SpaceTimeComponent(ComponentStorage_InMemory(CombinationElement(factors), CombinationHyperparameters(factor_hyperparameters)),
                                               SpaceTimeComponentSolutionStorage_InMemory())

    # Local component
    local_component = SpatialComponent(ComponentStorage_InMemory(LocalElement(n_triangulation_divisions=4), 
                                                                 LocalHyperparameters(log_sigma=0.0, log_rho=numpy.log(10.0 * numpy.pi/180))),
                                       SpatialComponentSolutionStorage_InMemory())

    print 'Analysing inputs'

    # Analysis system using the specified components, for the Tmean observable
    analysis_system = AnalysisSystem(
        [ climatology_component, large_scale_component, local_component ],
        ObservationSource.TMEAN)

    # Make filelist
    input_filelist = [ os.path.join(input_basepath, filename) for filename in input_filenames ]

    # Object to load HadCRUT4 inputs at time indices
    inputloader = AnalysisSystemInputLoaderHadCRUT4(input_filelist)

    # Update with data
    analysis_system.update([ inputloader ], time_indices)

    print 'Computing outputs'

    # Produce an output for each time index
    for time_index in time_indices:

        # Make output filename
        outputdate = inputloader.datetime_at_time_index(time_index)
        pathname = 'example_output_{0:04d}{1:02d}.nc'.format(outputdate.year, outputdate.month)
        print 'Saving: ', pathname

        # Configure output grid
        outputstructure = OutputRectilinearGridStructure(
            time_index, outputdate,
            latitudes=numpy.linspace(-87.5, 87.5, num=36),
            longitudes=numpy.linspace(-177.5, 177.5, num=72))

        # Evaluate expected value at these locations
        result_expected_value = analysis_system.evaluate_expected_value(outputstructure)

        # Save results
        filebuilder = FileBuilderHadCRUT4ExampleOutput(pathname, outputstructure)
        filebuilder.add_global_field(TAS_ANOMALY, result_expected_value.reshape(1,36,72))
        filebuilder.save_and_close()

    print 'Complete'

if __name__=='__main__':

    # Call main method
    main()
