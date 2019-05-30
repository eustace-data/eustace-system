"""Test of example using simulated data, low-resolution mesh, with kronecker and local elements."""

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
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.kronecker import SpaceTimeKroneckerElement
from eustace.analysis.advanced_standard.elements.spacetimespde import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory


class SimulatedObservationStructure(ObservationStructure):
    """Simulated input data structure."""

    def __init__(self, time_index, locations, measurement, uncorrelated_error_magnitude):

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
        
class TestMiniWorldKronecker(unittest.TestCase):

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
        element_large = SpaceTimeKroneckerElement(n_triangulation_divisions=1, alpha=2, starttime=0, endtime=number_of_simulated_time_steps+1, n_nodes=number_of_simulated_time_steps+2,overlap_factor=2.5, H=1)
        initial_hyperparameters_large = SpaceTimeSPDEHyperparameters(
                space_log_sigma=0.0,
                space_log_rho=numpy.log(numpy.radians(5.0)),
                time_log_rho=numpy.log(1.0 / 365.0))

        component_large = SpaceTimeComponent(
            ComponentStorage_InMemory(element_large, initial_hyperparameters_large),
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

        # All systems linear so single update should be ok
        analysis_system.update([ simulated_input_loader ], simulated_time_indices)
        
        # Get all results
        result = numpy.zeros(measurement.shape)
        for t in range(number_of_simulated_time_steps):
            result[t,:] = analysis_system.evaluate_expected_value('MAP', SimulatedObservationStructure(t, fixed_locations, None, None), flag='POINTWISE')

        disparity_large_scale = (numpy.abs(result - large_scale_process)).ravel().max()
        # print 'large scale disparity: ', disparity_large_scale

        disparity_overall = (numpy.abs(result - measurement)).ravel().max()
        # print 'overall disparity: ', disparity_overall

        numpy.testing.assert_almost_equal(result, measurement, decimal=4)
        self.assertTrue(disparity_overall < 1E-4)

        # numpy.testing.assert_almost_equal(result, large_scale_process, decimal=4)
