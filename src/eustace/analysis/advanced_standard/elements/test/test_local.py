import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalDesign
from eustace.analysis.advanced_standard.elements.local import LocalPrior
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.mesh.geometry import cartesian_to_polar2d

class TestLocalHyperparameters(unittest.TestCase):
    
    def test_init(self):
        
        h = LocalHyperparameters(log_sigma=23.2, log_rho=87.9)
        self.assertEqual(87.9, h.log_rho)
        self.assertEqual(23.2, h.log_sigma)
        numpy.testing.assert_almost_equal(h.get_array(), [ 23.2, 87.9 ])
               
class TestLocalElement(unittest.TestCase):
    
    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 15.0, -7.0 ], [  5.0, 100.0 ] ])                 
    
    def test_init(self):
        
        local = LocalElement(n_triangulation_divisions=2)
        self.assertEqual(2, local.spde.triangulation.level)

    def test_element_design(self):

        local = LocalElement(n_triangulation_divisions=1)
        design = local.element_design(TestLocalElement.SimulatedObservationStructure())
        self.assertTrue(isinstance(design, LocalDesign))
        self.assertEqual(1, design.spde.triangulation.level)

    def test_design_prior(self):

        local = LocalElement(n_triangulation_divisions=1)
        prior = local.element_prior(LocalHyperparameters(0.0, numpy.log(10.0 * numpy.pi/180)))
        self.assertTrue(isinstance(prior, LocalPrior))
        numpy.testing.assert_almost_equal(prior.hyperparameters.get_array(), [ 0.0, numpy.log(10.0 * numpy.pi/180) ])
        self.assertEqual(1, prior.spde.triangulation.level)
        self.assertEqual(2, prior.alpha)

class TestLocalDesign(unittest.TestCase):
                
    def test_design_number_of_state_parameters(self):

        spde = SphereMeshSPDE(level=1)
        design = LocalDesign(TestLocalElement.SimulatedObservationStructure(), spde)
        self.assertEqual(42, design.design_number_of_state_parameters())

    def test_design_matrix(self):

        spde = SphereMeshSPDE(level=1)
        design = LocalDesign(TestLocalElement.SimulatedObservationStructure(), spde)

        A = design.design_matrix()

        # Check sparse format
        self.assertEqual(SPARSEFORMAT, A.getformat())

        # Shape should be number of obs x number of vertices (basis functions)
        self.assertEqual((2, 42), A.shape)
        
        # Two triangles means 6 points are non-zero
        self.assertEqual(6, len(A.nonzero()[0]))
        
        # Get and check indices of nonzero design elements
        observation_indices, vertex_indices = A.sorted_indices().nonzero()
        numpy.testing.assert_equal(observation_indices, [ 0, 0, 0, 1, 1, 1 ])
        self.assertEqual((6,), vertex_indices.shape)

        # Vertices of observations 0 and 1
        # (one row per vertex, one column per coordinate)
        vertices0 = spde.triangulation.points[vertex_indices[0:3],:]
        vertices1 = spde.triangulation.points[vertex_indices[3:6],:]

        # Multiply vertices by the weights from A and sum to get cartesian locations
        testpoint0 = A[0, vertex_indices[0:3] ] * vertices0
        testpoint1 = A[1, vertex_indices[3:6] ] * vertices1

        # Check results correspond to original polar coordinates
        numpy.testing.assert_almost_equal(cartesian_to_polar2d(testpoint0), [ [ 15.0,  -7.0 ] ])
        numpy.testing.assert_almost_equal(cartesian_to_polar2d(testpoint1), [ [  5.0, 100.0 ] ])

    def test_design_function(self):

        # Build design
        design = LocalDesign(TestLocalElement.SimulatedObservationStructure(), SphereMeshSPDE(level=1))

        # Indices from design matrix indicate vertex weights that produce observation locations
        observation_indices, vertex_indices = design.design_matrix().sorted_indices().nonzero()

        # A vector that is 176.0 
        # except in location 0 where it is 82.0
        # and in location 1 where it is -8.888
        statevector = numpy.tile(176.0, (42,1))
        statevector[vertex_indices[0:3]] = 272.0
        statevector[vertex_indices[3:6]] = -8.888
        
        # Evaluate this - should get the two numbers back
        numpy.testing.assert_almost_equal(design.design_function(statevector), [ [ 272.0 ], [ -8.888 ] ])

    def test_isnonlinear(self):
        
        design = LocalDesign(TestLocalElement.SimulatedObservationStructure(), SphereMeshSPDE(level=1))
        self.assertFalse(design.isnonlinear())

    def test_design_jacobian(self):

        design = LocalDesign(TestLocalElement.SimulatedObservationStructure(), SphereMeshSPDE(level=1))
        numpy.testing.assert_equal(design.design_matrix().todense(), design.design_jacobian(currentstate=None).todense())
        
class TestLocalPrior(unittest.TestCase):

    def setUp(self):

        spde = SphereMeshSPDE(level=1)
        self.prior = LocalPrior(LocalHyperparameters(0.0, numpy.log(10.0 * numpy.pi/180)), spde)

    def test_prior_number_of_state_parameters(self):

        self.assertEqual(42, self.prior.prior_number_of_state_parameters())

    def test_prior_precision(self):
        
        Q = self.prior.prior_precision()
        self.assertEqual(SPARSEFORMAT, Q.getformat())
        self.assertEqual((42, 42), Q.shape)

    def test_prior_precision_derivative(self):
        
        dQ = self.prior.prior_precision_derivative(1)
        # self.assertTrue(scipy.sparse.isspmatrix_csc(dQ))
        self.assertEqual((42, 42), dQ.shape)
