"""Test of example using simulated data, and local component only."""

import unittest
import numpy
import os
import scipy.sparse
import  tempfile
from netCDF4 import Dataset
from StringIO import StringIO

from eustace.timeutils.epoch import epoch_plus_days

from eustace.analysis.observationsource import ObservationSource

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystemInputLoader

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeFunction
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.geography_based import GeographyBasedElement

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

from eustaceconfig import WORKSPACE_PATH

from netCDF4 import Dataset

ALTITUDE_RELATIVE_PATH = 'data/internal/climatology_covariates/DEM_global_0.25_0.25.nc'

class SimulatedObservationStructure(ObservationStructure):

    def __init__(self, time_index, locations, measurement, uncorrelated_error_magnitude):

        # Store time index
        self.t = time_index
        self.locations = locations
        self.measurement = measurement
        self.uncorrelated_error_magnitude  = uncorrelated_error_magnitude

    def time_index(self):

        return self.t
        
    def time_datetime(self):

        return epoch_plus_days(self.t)

    def number_of_observations(self):

        return self.locations.shape[0]
        
    def location_polar_coordinates(self):

        return self.locations
                
    def observation_vector(self):
        """Array of observations."""

        return self.measurement

    def observation_precision(self):
        """Observation precision matrix (sparse)."""

        return scipy.sparse.diags(1.0 / (self.uncorrelated_error_magnitude ** 2), format='csc')


class SimulatedInputLoader(AnalysisSystemInputLoader):
    """Simulated loading of input data."""

    def __init__(self, fixed_locations, measurement, uncorrelated_error_magnitude):

        self.fixed_locations = fixed_locations
        self.measurement = measurement
        self.uncorrelated_error_magnitude = uncorrelated_error_magnitude


    def datetime_at_time_index(self, time_index):
        """Use EUSTACE epoch."""

        return epoch_plus_days(time_index)
        
    def load_observation_structure(self, observable, time_index, log=None):
        """Load simulated data at specified time index."""

        if observable == ObservationSource.TMEAN:

            return SimulatedObservationStructure(time_index, self.fixed_locations, self.measurement, self.uncorrelated_error_magnitude)


class TestMiniWorldGeographyBased(unittest.TestCase):

    def setUp(self):
	"""Generate a geography-based covariate miniworld file"""
	self.covariate_file = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.elements.test_geography_based', suffix='.nc',delete=True)

	nx, ny = (5, 3)
        x = numpy.linspace(0, 1, nx) 
        y = numpy.linspace(0, 0.5, ny)
        self.longitude, self.latitude = numpy.meshgrid(x, y)
	self.covariate=numpy.ones([nx, ny])
	dataset=Dataset(self.covariate_file.name,'w','NETCDF4')
	dataset.createDimension('x',3)
	dataset.createDimension('y',5)
	field_latitude=dataset.createVariable('lat',numpy.float32,dimensions=('x','y'))
	field_longitude=dataset.createVariable('lon',numpy.float32,dimensions=('x','y'))
	field_covariate=dataset.createVariable('covariate',numpy.float32,dimensions=('x','y'))
	field_latitude[:]=self.latitude
	field_longitude[:]=self.longitude
	field_covariate[:]=self.covariate
	dataset.close()

	self.altitude_datafile = os.path.join(WORKSPACE_PATH, ALTITUDE_RELATIVE_PATH)

    def tearDown(self):
	self.latitude = None
	self.longitude = None
	self.covariate = None

    def test_mini_world_geography_based_mock_data(self):
	"""Testing on a simple mock data file, with mock covariate values"""
	
	# GENERATING OBSERVATIONS
	# Simulated locations: they will exactly sits on the grid points of the covariate datafile
        locations = numpy.array( [ [ 0.0, 0.0 ],
                                   [ 0.0, 0.5 ],
                                   [ 0.05, 0.0 ] ] )

        # Simulated measurements: simple linear relation of type: y = 2*x
        measurement = numpy.array( [ 2., 2., 2. ] )

        # Simulated errors
        uncorrelatederror = 0.1*numpy.ones(measurement.shape)

        # Simulated inputs
        simulated_input_loader = SimulatedInputLoader(locations, measurement, uncorrelatederror)

        # Simulate evaluation of this time index
        simulated_time_indices = [ 0 ]


        # GENERATING THE MODEL
        # Local component
	geography_covariate_element = GeographyBasedElement(self.covariate_file.name,'lat','lon','covariate', 1.0)
	geography_covariate_element.load()
        geography_based_component = SpatialComponent(ComponentStorage_InMemory(geography_covariate_element, CovariateHyperparameters(-0.5*numpy.log(10.))), SpatialComponentSolutionStorage_InMemory())

	# GENERATING THE ANALYSIS
        # Analysis system using the specified components, for the Tmean observable
        analysis_system = AnalysisSystem([ geography_based_component ], ObservationSource.TMEAN, log=StringIO())

        # Update with data
        analysis_system.update([ simulated_input_loader ], simulated_time_indices)

        # Check state vector directly
        statevector = analysis_system.components[0].solutionstorage.partial_state_read(0).ravel()
	
        # These are the nodes where observations were put (see SimulatedObservationSource above)
        # - check they correspond to within 3 times the stated noise level
        self.assertAlmostEqual( 2., statevector[0], delta=0.3)
        
	# Also check entire state vector within outer bounds set by obs
        self.assertTrue(all(statevector < 2.0))
	
        # And check output corresponds too
        # (evaluate result on output structure same as input)
        simulated_output_structure = SimulatedObservationStructure(0, locations, None, None)
        result = analysis_system.evaluate_expected_value('MAP', simulated_output_structure, flag='POINTWISE')
        numpy.testing.assert_almost_equal(statevector[ 0 ]*numpy.ones(len(measurement)), result)

    def test_mini_world_latitude_harmonics(self):
	"""Testing on a simple mock data file using latitude harmonics"""

	# GENERATING OBSERVATIONS
	# Simulated locations: they will exactly sits on the grid points of the covariate datafile
        locations = numpy.array( [ [ 0.0, 0.0 ],
                                   [ 0.25, 0.5 ],
                                   [ 0.5, 0.0 ] ] )
	# Simulated model is y = a*cos(2x) + c*cos(4*x) + b*sin(2x) + d*sin(4*x) with x = latitude, so we expect a=c=1, c=d=0
        measurement = LatitudeFunction(numpy.cos, 2.0).compute(locations[:,0]).ravel() + LatitudeFunction(numpy.cos, 4.0).compute(locations[:,0]).ravel()

        # Simulated errors
        uncorrelatederror = 0.1*numpy.ones(measurement.shape)

        # Simulated inputs
        simulated_input_loader = SimulatedInputLoader(locations, measurement, uncorrelatederror)

        # Simulate evaluation of this time index
        simulated_time_indices = [ 0 ]

        latitude_harmonics_component = SpatialComponent(ComponentStorage_InMemory(LatitudeHarmonicsElement(), 
					      CombinationHyperparameters([ CovariateHyperparameters(-0.5*numpy.log(p)) for p in [ 10.0, 10.0, 10.0, 10.0 ] ])), 
                                              SpatialComponentSolutionStorage_InMemory()) 

        # Analysis system using the specified components, for the Tmean observable
        analysis_system = AnalysisSystem([ latitude_harmonics_component ], ObservationSource.TMEAN, log=StringIO())

	# GENERATING THE ANALYSIS

        # Update with data
        analysis_system.update([ simulated_input_loader ], simulated_time_indices)

        # Check state vector directly
        statevector = analysis_system.components[0].solutionstorage.partial_state_read(0).ravel()

        # These are the nodes where observations were put (see SimulatedObservationSource above)
        # - check they correspond to within 3 times the stated noise level
        self.assertAlmostEqual( 1., statevector[0], delta=0.3)
        self.assertAlmostEqual( 1., statevector[2], delta=0.3)
        self.assertAlmostEqual( 0., statevector[1], delta=0.3)
        self.assertAlmostEqual( 0., statevector[3], delta=0.3)
	# Also check entire state vector within outer bounds set by obs
        self.assertTrue(all(statevector < 1.0))
	
        # And check output corresponds too
        # (evaluate result on output structure same as input)
        simulated_output_structure = SimulatedObservationStructure(0, locations, None, None)
        result = analysis_system.evaluate_expected_value('MAP', simulated_output_structure, flag='POINTWISE')
	expected = statevector[0]*LatitudeFunction(numpy.cos, 2.0).compute(locations[:,0]).ravel() + statevector[1]*LatitudeFunction(numpy.sin, 2.0).compute(locations[:,0]).ravel()\
                  + statevector[2] *LatitudeFunction(numpy.cos, 4.0).compute(locations[:,0]).ravel()+ statevector[3]*LatitudeFunction(numpy.sin, 4.0).compute(locations[:,0]).ravel()
        numpy.testing.assert_almost_equal(expected, result)

    def test_mini_world_altitude(self):
	"""Testing using altitude as a covariate"""
	
	# GENERATING OBSERVATIONS
	# Simulated locations: they will exactly sits on the grid points of the covariate datafile
	DEM = Dataset(self.altitude_datafile)
	latitude = DEM.variables['lat'][:]
	longitude = DEM.variables['lon'][:]
	altitude = DEM.variables['dem'][:]

	indices = numpy.stack((numpy.array([1, 3, 267, 80, 10, 215, 17, 120]), numpy.array([2, 256, 9, 110, 290, 154, 34, 151])), axis=1)

	selected_location = []
	altitude_observations = []
	for couple in indices:
	  selected_location.append([latitude[couple[0],couple[1]], longitude[couple[0], couple[1]]])
	  altitude_observations.append(altitude[couple[0],couple[1]])
	DEM.close()

	locations = numpy.array(selected_location)
        # Simulated measurements: simple linear relation of type: y = PI*x
	measurement = numpy.pi*numpy.array(altitude_observations)

        # Simulated errors
        uncorrelatederror = 0.1*numpy.ones(measurement.shape)

        # Simulated inputs
        simulated_input_loader = SimulatedInputLoader(locations, measurement, uncorrelatederror)

        # Simulate evaluation of this time index
        simulated_time_indices = [ 0 ]

        # GENERATING THE MODEL
        # Local component
	geography_covariate_element = GeographyBasedElement(self.altitude_datafile,'lat','lon','dem', 1.0)
	geography_covariate_element.load()
        geography_based_component = SpatialComponent(ComponentStorage_InMemory(geography_covariate_element, CovariateHyperparameters(-0.5*numpy.log(10.))), SpatialComponentSolutionStorage_InMemory())

	# GENERATING THE ANALYSIS
        # Analysis system using the specified components, for the Tmean observable
        analysis_system = AnalysisSystem([ geography_based_component ], ObservationSource.TMEAN, log=StringIO())

        # Update with data
        analysis_system.update([ simulated_input_loader ], simulated_time_indices)

        # Check state vector directly
        statevector = analysis_system.components[0].solutionstorage.partial_state_read(0).ravel()

        # These are the nodes where observations were put (see SimulatedObservationSource above)
        # - check they correspond to within 3 times the stated noise level
        self.assertAlmostEqual( numpy.pi, statevector[0], delta=0.3)
        
	# Also check entire state vector within outer bounds set by obs
        self.assertTrue(all(statevector < numpy.pi))

        # And check output corresponds too
        # (evaluate result on output structure same as input)
        simulated_output_structure = SimulatedObservationStructure(0, locations, None, None)
        result = analysis_system.evaluate_expected_value('MAP', simulated_output_structure, flag='POINTWISE')
        numpy.testing.assert_almost_equal(statevector[0]*numpy.array(altitude_observations), result)

    def test_mini_world_altitude_with_latitude(self):
	"""Testing using altitude as a covariate"""
	
	# GENERATING OBSERVATIONS
	# Simulated locations: they will exactly sits on the grid points of the covariate datafile
	DEM = Dataset(self.altitude_datafile)
	latitude = DEM.variables['lat'][:]
	longitude = DEM.variables['lon'][:]
	altitude = DEM.variables['dem'][:]

	indices = numpy.stack((numpy.array([1, 3, 5, 7, 8, 9, 10, 11]), numpy.array([0, 0, 0, 0, 0, 0, 0, 0])), axis=1)

	selected_location = []
	altitude_observations = []
	for couple in indices:
	  selected_location.append([latitude[couple[0],couple[1]], longitude[couple[0], couple[1]]])
	  altitude_observations.append(altitude[couple[0],couple[1]])
	DEM.close()

	locations = numpy.array(selected_location)
	# Simulated model is y = z + a*cos(2x) + c*cos(4*x) + b*sin(2x) + d*sin(4*x), with z = altitude, x = latitude, a=b=c=d=0
	slope = 1e-3
	measurement = slope*numpy.array(altitude_observations) 

        # Simulated errors
        uncorrelatederror = 0.1*numpy.ones(measurement.shape)

        # Simulated inputs
        simulated_input_loader = SimulatedInputLoader(locations, measurement, uncorrelatederror)

        # Simulate evaluation of this time index
        simulated_time_indices = [ 0 ]

        # GENERATING THE MODEL
        # Local component
	geography_covariate_element = GeographyBasedElement(self.altitude_datafile,'lat','lon','dem', 1.0)
	geography_covariate_element.load()
	combined_element = CombinationElement([geography_covariate_element, LatitudeHarmonicsElement()])
	combined_hyperparamters = CombinationHyperparameters([ CovariateHyperparameters(-0.5*numpy.log(10.)), CombinationHyperparameters([ CovariateHyperparameters(-0.5*numpy.log(p)) for p in [ 10.0, 10.0, 10.0, 10.0 ] ])])
	combined_component = SpatialComponent(ComponentStorage_InMemory(combined_element, combined_hyperparamters), SpatialComponentSolutionStorage_InMemory())

	# GENERATING THE ANALYSIS
        # Analysis system using the specified components, for the Tmean observable
        analysis_system = AnalysisSystem([ combined_component ], ObservationSource.TMEAN, log=StringIO())

        # Update with data
        analysis_system.update([ simulated_input_loader ], simulated_time_indices)

        # Check state vector directly
        statevector = analysis_system.components[0].solutionstorage.partial_state_read(0).ravel()

        # These are the nodes where observations were put (see SimulatedObservationSource above)
        # - check they correspond to within 3 times the stated noise level
        self.assertAlmostEqual( slope, statevector[0], delta=0.3)
        self.assertAlmostEqual( 0., statevector[1], delta=0.3)
        self.assertAlmostEqual( 0., statevector[2], delta=0.3)
        self.assertAlmostEqual( 0., statevector[3], delta=0.3)
        self.assertAlmostEqual( 0., statevector[4], delta=0.3)

        # And check output corresponds too
        # (evaluate result on output structure same as input)
        simulated_output_structure = SimulatedObservationStructure(0, locations, None, None)
        result = analysis_system.evaluate_expected_value('MAP', simulated_output_structure, flag='POINTWISE')
	expected = statevector[0]*numpy.array(altitude_observations)\
                 + statevector[1]*LatitudeFunction(numpy.cos, 2.0).compute(locations[:,0]).ravel()\
                 + statevector[2]*LatitudeFunction(numpy.sin, 2.0).compute(locations[:,0]).ravel()\
                 + statevector[3]*LatitudeFunction(numpy.cos, 4.0).compute(locations[:,0]).ravel()\
                 + statevector[4]*LatitudeFunction(numpy.sin, 2.0).compute(locations[:,0]).ravel()
        numpy.testing.assert_almost_equal(expected, result)
        
        # test output gridding, pointwise limit
	outputstructure = OutputRectilinearGridStructure(2,  epoch_plus_days(2), latitudes=numpy.linspace(-60., 60., num=5), longitudes=numpy.linspace(-90., 90, num=10))
	pointwise_result = analysis_system.evaluate_expected_value('MAP', outputstructure, 'POINTWISE')
	pointwise_limit_result = analysis_system.evaluate_expected_value('MAP', outputstructure,  'GRID_CELL_AREA_AVERAGE', [1, 1], 10)
	numpy.testing.assert_array_almost_equal(pointwise_result, pointwise_limit_result)
