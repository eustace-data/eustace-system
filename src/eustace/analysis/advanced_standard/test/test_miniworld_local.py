"""Test of example using simulated data, low-resolution mesh, and local component only."""

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
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

class SimulatedObservationStructure(ObservationStructure):

    def __init__(self, time_index):

        # Store time index
        self.t = time_index

        # Inspect the same mesh as used for output
        # mesh = MeshIcosahedronSubdivision.build(number_of_subdivisions=1)
        # print cartesian_to_polar2d(mesh.points)
        #
        # Index 12: [   0,    0  ]
        # Index 17: [   0,  180  ]
        # Index 41: [   0,   -90 ]
        
        # Simulated locations
        self.locations = numpy.array( [ [ 0.0,   0.0 ],
                                        [ 0.0, 180.0 ],
                                        [ 0.0, -90.0 ] ] )

        # Simulated measurements
        self.measurement = numpy.array( [ 20.0, -15.0, 5.0 ] )

        # Simulated errors
        self.uncorrelatederror = 0.1*numpy.ones(self.measurement.shape)


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

        return scipy.sparse.diags(1.0 / (self.uncorrelatederror ** 2), format='csc')


class SimulatedInputLoader(AnalysisSystemInputLoader):
    """Simulated loading of input data."""

    def __init__(self):

        pass

    def datetime_at_time_index(self, time_index):
        """Use EUSTACE epoch."""

        return epoch_plus_days(time_index)
        
    def load_observation_structure(self, observable, time_index, log=None):
        """Load simulated data at specified time index."""

        if observable == ObservationSource.TMEAN:

            return SimulatedObservationStructure(time_index)
        
class TestMiniWorldLocal(unittest.TestCase):

    def test_mini_world_local(self):

        # Local component
        local_component = SpatialComponent(ComponentStorage_InMemory(LocalElement(n_triangulation_divisions=1), 
                                                                     LocalHyperparameters(log_sigma=0.0, log_rho=numpy.log(1.0))),
                                           SpatialComponentSolutionStorage_InMemory(), compute_uncertainties=True, method='APPROXIMATED')

        # Analysis system using the specified components, for the Tmean observable
        analysis_system = AnalysisSystem([ local_component ], ObservationSource.TMEAN, log=StringIO())

        # Simulated inputs
        simulated_input_loader = SimulatedInputLoader()

        # Simulate evaluation of this time index
        simulated_time_indices = [ 0 ]

        # Update with data
        analysis_system.update([ simulated_input_loader ], simulated_time_indices)

        # Check state vector directly
        statevector = analysis_system.components[0].solutionstorage.partial_state_read(0).ravel()
        # These are the nodes where observations were put (see SimulatedObservationSource above)
        # - check they correspond to within 3 times the stated noise level
        self.assertAlmostEqual( 20.0, statevector[12], delta=0.3)
        self.assertAlmostEqual(-15.0, statevector[17], delta=0.3)
        self.assertAlmostEqual(  5.0, statevector[41], delta=0.3)

        # Also check entire state vector within outer bounds set by obs
        self.assertTrue(all(statevector < 20.0))
        self.assertTrue(all(statevector >-15.0))

        # And check output corresponds too
        # (evaluate result on output structure same as input)
        simulated_output_structure = SimulatedObservationStructure(0)
        result = analysis_system.evaluate_expected_value('MAP', simulated_output_structure, flag='POINTWISE')
        numpy.testing.assert_almost_equal(statevector[ [ 12, 17, 41 ] ], result)

	# test output gridding, pointwise limit
	outputstructure = OutputRectilinearGridStructure(2,  epoch_plus_days(2), latitudes=numpy.linspace(-89.875, 89.875, num=10), longitudes=numpy.linspace(-179.875, 179.875, num=20))
	pointwise_result = analysis_system.evaluate_expected_value('MAP', outputstructure, 'POINTWISE')
	pointwise_limit_result = analysis_system.evaluate_expected_value('MAP', outputstructure,  'GRID_CELL_AREA_AVERAGE', [1, 1], 3)
	numpy.testing.assert_array_almost_equal(pointwise_result, pointwise_limit_result)