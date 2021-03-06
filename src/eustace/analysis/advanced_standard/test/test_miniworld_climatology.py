"""Test of example using simulated data, low-resolution mesh, and climatology model."""

import unittest
import numpy
import scipy.sparse
from StringIO import StringIO

from eustace.timeutils.epoch import epoch_plus_days

from eustace.analysis.mesh.mesh import MeshIcosahedronSubdivision
from eustace.analysis.mesh.geometry import cartesian_to_polar2d
from eustace.analysis.observationsource import ObservationSource

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystemInputLoader

from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters
from eustace.analysis.advanced_standard.elements.seasonal import datetime_to_decimal_year
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

class SimulatedObservationStructure(ObservationStructure):
    """Simulated input data structure."""

    def __init__(self, time_index, locations, measurement, uncorrelated_error_magnitude):

        self.t = time_index
        self.locations = locations
        self.measurement = measurement
        self.uncorrelated_error_magnitude  = uncorrelated_error_magnitude
	self.time_index_number = time_index
	self.corresponding_datetime = None
	
    def time_index(self):

        return self.t
        
    def time_datetime(self):

        return epoch_plus_days(self.t)

    def number_of_observations(self):

        return self.locations.shape[0]
        
    def location_polar_coordinates(self):

        return self.locations
                
    def observation_vector(self):

        return self.measurement

    def observation_precision(self):

        # A diagonal matrix with diagonal equal to 1 / sigma^2
        return scipy.sparse.eye(self.locations.shape[0], format='csc') * (1.0 / (self.uncorrelated_error_magnitude ** 2))


class SimulatedInputLoader(AnalysisSystemInputLoader):
    """Simulated loading of input data."""

    def __init__(self, fixed_locations, measurement, uncorrelated_error_magnitude):

        self.fixed_locations = fixed_locations
        self.measurement = measurement
        self.uncorrelated_error_magnitude = uncorrelated_error_magnitude

    def datetime_at_time_index(self, time_index):
        return epoch_plus_days(time_index)
        
    def load_observation_structure(self, observable, time_index, log=None):

        if observable == ObservationSource.TMEAN:

            return SimulatedObservationStructure(time_index, self.fixed_locations, self.measurement[time_index,:], self.uncorrelated_error_magnitude)
        
class TestMiniWorldClimatology(unittest.TestCase):

    def test_mini_world_noiseless(self):


	number_of_simulated_time_steps = 1

        # Build system
	element = SeasonalElement(n_triangulation_divisions=3, n_harmonics=5, include_local_mean=True) 
        hyperparameters = SeasonalHyperparameters(n_spatial_components=6, common_log_sigma=0.0, common_log_rho=0.0)
 
        component = SpaceTimeComponent(ComponentStorage_InMemory(element, hyperparameters), SpaceTimeComponentSolutionStorage_InMemory())

        analysis_system = AnalysisSystem([ component ], ObservationSource.TMEAN, log=StringIO())

        # use fixed locations from icosahedron
        fixed_locations = cartesian_to_polar2d( MeshIcosahedronSubdivision.build(3).points )

        # random measurement at each location
        numpy.random.seed(8976)
        field_basis = numpy.random.randn(fixed_locations.shape[0])
	#print(field_basis.shape)
        #time_basis = numpy.array(harmonics_list)
        # some time function that varies over a year
	#decimal_years = numpy.array([datetime_to_decimal_year(epoch_plus_days(step)) for step in range(number_of_simulated_time_steps)])
        time_basis = numpy.cos(numpy.linspace(0.1, 1.75*numpy.pi,  number_of_simulated_time_steps))
        # kronecker product of the two
	#print(numpy.expand_dims(time_basis, 1))
        measurement = numpy.kron(field_basis, numpy.expand_dims(time_basis, 1))#numpy.expand_dims(time_basis, 1))

	#print(measurement.shape)
        # Simulated inputs
        simulated_input_loader = SimulatedInputLoader(fixed_locations, measurement, 0.0001 )

        # Simulate evaluation of this time index
        simulated_time_indices = range(number_of_simulated_time_steps)

        # Iterate
        for iteration in range(5):
            analysis_system.update([ simulated_input_loader ], simulated_time_indices)
        
        # Get all results
        result = numpy.zeros(measurement.shape)
        for t in range(number_of_simulated_time_steps):
            result[t,:] = analysis_system.evaluate_expected_value('MAP', SimulatedObservationStructure(t, fixed_locations, None, None), flag='POINTWISE')
        # Should be very close to original because specified noise is low
        numpy.testing.assert_almost_equal(result, measurement)
        max_disparity = (numpy.abs(result - measurement)).ravel().max()
        self.assertTrue(max_disparity < 1E-5)

	# test output gridding, pointwise limit
	outputstructure = OutputRectilinearGridStructure(2,  epoch_plus_days(2), latitudes=numpy.linspace(-60., 60., num=5), longitudes=numpy.linspace(-90., 90, num=10))
	pointwise_result = analysis_system.evaluate_expected_value('MAP',outputstructure, 'POINTWISE')
	pointwise_limit_result = analysis_system.evaluate_expected_value('MAP', outputstructure,  'GRID_CELL_AREA_AVERAGE', [1, 1], 10)
	numpy.testing.assert_array_almost_equal(pointwise_result, pointwise_limit_result)
	
"""
    def test_mini_world_large_and_local(self):

        # Use a number of time steps
        number_of_simulated_time_steps = 30

        # Large-scale spatial variability
        simulated_large_variation = 10.0

        # Local variability
        simulated_local_variation = 1.0

        # Iterations to use
        number_of_solution_iterations = 5

        # Build system

        # Large-scale factor
        element_large = SpaceTimeFactorElement(n_triangulation_divisions=3, alpha=2, starttime=0, endtime=number_of_simulated_time_steps+1, overlap_factor=2.5, H=1)
        hyperparameters_large = SpaceTimeSPDEHyperparameters(
                space_log_sigma=0.0,
                space_log_rho=numpy.log(numpy.radians(5.0)),
                time_log_rho=numpy.log(1.0 / 365.0))
        component_large = SpaceTimeComponent(
            ComponentStorage_InMemory(element_large, hyperparameters_large),
            SpaceTimeComponentSolutionStorage_InMemory())

        # And a local process
        component_local = SpatialComponent(
            ComponentStorage_InMemory(LocalElement(n_triangulation_divisions=3), 
                                      LocalHyperparameters(log_sigma=0.0, log_rho=numpy.log(numpy.radians(5.0)))),
            SpatialComponentSolutionStorage_InMemory())

        analysis_system = AnalysisSystem([ component_large, component_local ], ObservationSource.TMEAN, log=StringIO())
        # analysis_system = AnalysisSystem([ component_large ], ObservationSource.TMEAN)
        # analysis_system = AnalysisSystem([ component_local ], ObservationSource.TMEAN)

        # use fixed locations from icosahedron
        fixed_locations = cartesian_to_polar2d( MeshIcosahedronSubdivision.build(3).points )

        # random measurement at each location
        numpy.random.seed(8976)
        field_basis = simulated_large_variation*numpy.random.randn(fixed_locations.shape[0])

        # some time function that varies over a year
        time_basis = numpy.cos(numpy.linspace(0.1, 1.75*numpy.pi,  number_of_simulated_time_steps))

        # kronecker product of the two
        large_scale_process = numpy.kron(field_basis, numpy.expand_dims(time_basis, 1))

        # Random local changes where mean change at each time is zero
        # local_process = simulated_local_variation * numpy.random.randn(large_scale_process.shape[0], large_scale_process.shape[1])
        # local_process -= numpy.tile(local_process.mean(axis=1), (local_process.shape[1], 1)).T

        local_process = numpy.zeros(large_scale_process.shape)
        somefield = simulated_local_variation * numpy.random.randn(1, large_scale_process.shape[1])
        somefield -= somefield.ravel().mean()
        local_process[10,:] = somefield
        local_process[11,:] = -somefield

        # Add the two processes
        measurement = large_scale_process + local_process

        # Simulated inputs
        simulated_input_loader = SimulatedInputLoader(fixed_locations, measurement, 0.001)

        # Simulate evaluation of this time index
        simulated_time_indices = range(number_of_simulated_time_steps)

        # Iterate over all component solutions
        #for iteration in range(number_of_solution_iterations):    
        #    analysis_system.update([ simulated_input_loader ], simulated_time_indices)

        # Iterate for non-linear large scale component
        for iteration in range(number_of_solution_iterations):
          analysis_system.update_component([ simulated_input_loader ], 0, simulated_time_indices)

        # Single solution for local component
        analysis_system.update_component([ simulated_input_loader ], 1, simulated_time_indices)
        
        # Get all results
        result = numpy.zeros(measurement.shape)
        for t in range(number_of_simulated_time_steps):
            result[t,:] = analysis_system.evaluate_expected_value(SimulatedObservationStructure(t, fixed_locations, None, None))

        disparity_large_scale = (numpy.abs(result - large_scale_process)).ravel().max()
        # print 'large scale disparity: ', disparity_large_scale

        disparity_overall = (numpy.abs(result - measurement)).ravel().max()
        # print 'overall disparity: ', disparity_overall

        numpy.testing.assert_almost_equal(result, measurement, decimal=4)
        self.assertTrue(disparity_overall < 1E-5)        

        # numpy.testing.assert_almost_equal(result, large_scale_process, decimal=4)
"""
