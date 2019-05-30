"""Tests for seasonal element."""

import unittest
import numpy
import scipy.sparse
from datetime import datetime

from eustace.analysis.advanced_standard.elements.seasonal import datetime_to_decimal_year
from eustace.analysis.advanced_standard.elements.seasonal import sine_seasonal_harmonic
from eustace.analysis.advanced_standard.elements.seasonal import cosine_seasonal_harmonic
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElementDesign
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElementPrior
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.mesh.geometry import cartesian_to_polar2d
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT

class TestTemporalBasisFunctions(unittest.TestCase):

    def test_datetime_to_decimal_year(self):

        # easy cases
        self.assertAlmostEqual(1789.0, datetime_to_decimal_year(datetime(1789, 1, 1, 0, 0, 0)))
        self.assertAlmostEqual(2017.0, datetime_to_decimal_year(datetime(2017, 1, 1)))

        # note 365 = 73 x 5 so in non-leap year 15th March is 0.2 of year
        self.assertAlmostEqual(1977.2, datetime_to_decimal_year(datetime(1977, 3, 15)))

        # 2nd July is start of day 183 (zero-based) in leap year, which is exactly half-way through 
        self.assertAlmostEqual(1904.5, datetime_to_decimal_year(datetime(1904, 7, 2)))

        # start of last day is 364/365 in non-leap year
        self.assertAlmostEqual(2015 + (364.0 / 365.0), datetime_to_decimal_year(datetime(2015, 12, 31)))

        # start of last day is 365/366 in leap year
        self.assertAlmostEqual(2016 + (365.0 / 366.0), datetime_to_decimal_year(datetime(2016, 12, 31)))

        # 12 noon should be 0.5/365 of year
        self.assertAlmostEqual(1850 + (0.5 / 365.0), datetime_to_decimal_year(datetime(1850, 1, 1, 12, 0, 0)))
        

    def test_cosine_seasonal_harmonic(self):

        self.assertAlmostEqual(1.0, cosine_seasonal_harmonic(1926.0, 1))
        self.assertAlmostEqual(1.0, cosine_seasonal_harmonic(1926.333333333, 3))
        self.assertAlmostEqual(1.0 / numpy.sqrt(2.0), cosine_seasonal_harmonic(1926.125 , 1))
        self.assertAlmostEqual(0.0, cosine_seasonal_harmonic(1926.125, 2))
        self.assertAlmostEqual(1.0 / numpy.sqrt(2.0), cosine_seasonal_harmonic(1926.0625, 2))

    def test_sine_seasonal_harmonic(self):

        self.assertAlmostEqual(0.0, sine_seasonal_harmonic(1926.0, 1))
        self.assertAlmostEqual(0.0, sine_seasonal_harmonic(1926.333333333, 3))
        self.assertAlmostEqual(1.0 / numpy.sqrt(2.0), sine_seasonal_harmonic(1926.125 , 1))
        self.assertAlmostEqual(1.0, sine_seasonal_harmonic(1926.125, 2))
        self.assertAlmostEqual(1.0 / numpy.sqrt(2.0), sine_seasonal_harmonic(1926.0625, 2))

class TestSeasonalHyperparameters(unittest.TestCase):

    def test_init(self):

        h = SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=1.1, common_log_rho=1.2)
        numpy.testing.assert_equal(h.harmonics, [ 1.1, 1.2, 1.1, 1.2, 1.1, 1.2 ])

    def test_get_array(self):

        h = SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0)
        numpy.testing.assert_equal(h.get_array(), [ 0.8, 1.0, 0.8, 1.0, 0.8, 1.0 ])

    def test_set_array(self):

        h = SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0)
        h.set_array([ 88.2, 103.4, 0.2, 0.3, 0.4, 0.5 ])
        numpy.testing.assert_equal(h.harmonics, [ 88.2, 103.4, 0.2, 0.3, 0.4, 0.5 ])

    def test_harmonic_parameters(self):

        h = SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0)
        h.set_array([ 88.2, 103.4, 0.2, 0.3, 0.4, 0.5 ])
        self.assertEqual(0.2, h.harmonic_parameters(1).log_sigma)
        self.assertEqual(103.4, h.harmonic_parameters(0).log_rho)
        self.assertEqual({ 'log_sigma': 0.4, 'log_rho': 0.5 }, h.harmonic_parameters(2).get_dictionary())


class TestSeasonalElement(unittest.TestCase):

    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 15.0, -7.0 ], [  5.0, 100.0 ] ])      

        def time_datetime(self):

            return datetime(1977, 3, 15)
          

    def test_init(self):

        s = SeasonalElement(n_triangulation_divisions=1, n_harmonics=3, include_local_mean=True)
        self.assertEqual(1, s.spde.triangulation.level)
        self.assertEqual(3, s.n_harmonics)
        self.assertEqual(True, s.include_local_mean)
        self.assertEqual(2, s.alpha)

    def test_element_design(self):

        element = SeasonalElement(n_triangulation_divisions=1, n_harmonics=3, include_local_mean=True)
        design = element.element_design(TestSeasonalElement.SimulatedObservationStructure())
        self.assertTrue(isinstance(design, SeasonalElementDesign))
        self.assertTrue(isinstance(design.observationstructure, TestSeasonalElement.SimulatedObservationStructure))
        self.assertEqual(1, design.spde.triangulation.level)
        self.assertEqual(3, design.n_harmonics)
        self.assertEqual(True, design.include_local_mean)

    def test_element_prior(self):

        element = SeasonalElement(n_triangulation_divisions=1, n_harmonics=3, include_local_mean=False)
        prior = element.element_prior(SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.5, common_log_rho=0.7))
        self.assertTrue(isinstance(prior, SeasonalElementPrior))
        numpy.testing.assert_equal([ 0.5, 0.7, 0.5, 0.7, 0.5, 0.7 ], prior.hyperparameters.get_array())
        self.assertEqual(3, prior.n_harmonics)
        self.assertEqual(False, prior.include_local_mean)
        self.assertEqual(2, prior.alpha)
        
class TestSeasonalElementDesign(unittest.TestCase):

    def test_init(self):

        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), SphereMeshSPDE(level=2), n_harmonics=5, include_local_mean=True)
        self.assertTrue(isinstance(design.observationstructure, TestSeasonalElement.SimulatedObservationStructure))
        self.assertEqual(2, design.spde.triangulation.level)
        self.assertEqual(5, design.n_harmonics)
        self.assertEqual(True, design.include_local_mean)

    def test_design_number_of_state_parameters(self):

        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), SphereMeshSPDE(level=1), n_harmonics=3, include_local_mean=True)
        self.assertEqual(294, design.design_number_of_state_parameters())
    
    def test_build_A_time(self):

        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), SphereMeshSPDE(level=2), n_harmonics=3, include_local_mean=True)
        A = design.build_A_time()

        # Simulated datetime(1977, 3, 15) is precisely 0.2 through the year (which is 72 degrees of 360 degree cycle)
        # Should give coefficients:
        #   1.0
        #   sin(72 degrees)
        #   cos(72 degrees)
        #   sin(144 degrees)
        #   cos(144 degrees)        
        #   sin(216 degrees)
        #   cos(216 degrees)        
        numpy.testing.assert_almost_equal(A, 
                                          [ [ 1.0,
                                              0.951056516295,
                                              0.309016994375,
                                              0.587785252292,
                                             -0.809016994375,
                                             -0.587785252292,
                                             -0.809016994375 ] ])
        
    def test_build_A_space(self):

        spde = SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT)
        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), spde, n_harmonics=3, include_local_mean=True)
        A = design.build_A_space()

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

    def test_design_matrix(self):

        spde = SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT)
        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), spde, n_harmonics=3, include_local_mean=True)
        A = design.design_matrix()

        # Check shape and sparse format
        self.assertEqual((2, 294), A.shape)
        self.assertEqual(SPARSEFORMAT, A.getformat())

        # Check that kronecker expansion worked as expected
        numpy.testing.assert_almost_equal(A[:,0:42].todense() *  0.951056516295, A[:,  42:84].todense())
        numpy.testing.assert_almost_equal(A[:,0:42].todense() * -0.809016994375, A[:,252:294].todense())
        
    def test_design_function(self):

        spde = SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT)
        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), spde, n_harmonics=3, include_local_mean=True)

        # Initialise zero state vector
        state_vector = numpy.zeros((294,1))

        # Get indices of nonzero design elements
        A_simulated = spde.build_A(TestSeasonalElement.SimulatedObservationStructure().location_polar_coordinates())
        observation_indices, vertex_indices = A_simulated.sorted_indices().nonzero()

        # This should select observation 1
        state_vector[vertex_indices[3:6],0] = 1.0

        # And this will add 0.309016994375 to it
        state_vector[vertex_indices[3:6]+84,0] = 1.0

        # Apply state vector
        y = design.design_function(state_vector)

        # Check observation 1 selected with appropriate sum
        numpy.testing.assert_almost_equal(y, numpy.array([ [ 0.0 ], [ 1.309016994375 ] ]))

    def test_isnonlinear(self):
        
        spde = SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT)
        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), spde, n_harmonics=3, include_local_mean=True)
        self.assertFalse(design.isnonlinear())

    def test_design_jacobian(self):

        spde = SphereMeshSPDE(level=2, sparse_format=SPARSEFORMAT)
        design = SeasonalElementDesign(TestSeasonalElement.SimulatedObservationStructure(), spde, n_harmonics=2, include_local_mean=False)
        A = design.design_matrix()
        J = design.design_jacobian(numpy.array([ ]))
        self.assertEqual(SPARSEFORMAT, J.getformat())
        numpy.testing.assert_almost_equal(J.todense(), A.todense())


class TestSeasonalElementPrior(unittest.TestCase):

    
    def test_init(self):

        prior = SeasonalElementPrior(
            SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0),
            spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
            n_harmonics=2, include_local_mean=True, alpha=3)

        numpy.testing.assert_equal(prior.hyperparameters.get_array(), [ 0.8, 1.0, 0.8, 1.0, 0.8, 1.0 ])
        self.assertEqual(2, prior.n_harmonics)
        self.assertTrue(prior.include_local_mean)
        self.assertEqual(3, prior.alpha)

    def test_n_spatial_components(self):

        # With local mean 2 harmonics plus mean = (1 + 2*2) = 5 parameters
        self.assertEqual(5, SeasonalElementPrior(
                SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0),
                spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
                n_harmonics=2, include_local_mean=True, alpha=2).n_spatial_components())

        # Without local mean 2 harmonics = 2*2 = 4 parameters
        self.assertEqual(4, SeasonalElementPrior(
                SeasonalHyperparameters(n_spatial_components=2, common_log_sigma=0.8, common_log_rho=1.0),
                spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
                n_harmonics=2, include_local_mean=False, alpha=2).n_spatial_components())

        # Also check state parameters
        self.assertEqual(210, SeasonalElementPrior(
                SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0),
                spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
                n_harmonics=2, include_local_mean=True, alpha=2).prior_number_of_state_parameters())

    def test_harmonic_index_for_spatial_component(self):
        
        prior_a = SeasonalElementPrior(
            SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0),
            spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
            n_harmonics=2, include_local_mean=True, alpha=2)       
        self.assertEqual(0, prior_a.harmonic_index_for_spatial_component(0))
        self.assertEqual(1, prior_a.harmonic_index_for_spatial_component(1))
        self.assertEqual(1, prior_a.harmonic_index_for_spatial_component(2))
        self.assertEqual(2, prior_a.harmonic_index_for_spatial_component(3))
        self.assertEqual(2, prior_a.harmonic_index_for_spatial_component(4))

        prior_b = SeasonalElementPrior(
            SeasonalHyperparameters(n_spatial_components=2, common_log_sigma=0.8, common_log_rho=1.0),
            spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
            n_harmonics=2, include_local_mean=False, alpha=2)
        self.assertEqual(0, prior_b.harmonic_index_for_spatial_component(0))
        self.assertEqual(0, prior_b.harmonic_index_for_spatial_component(1))
        self.assertEqual(1, prior_b.harmonic_index_for_spatial_component(2))
        self.assertEqual(1, prior_b.harmonic_index_for_spatial_component(3))
    

    def test_prior_precision(self):

        prior = SeasonalElementPrior(
            SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0),
            spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
            n_harmonics=2, include_local_mean=True, alpha=2)

        Q = prior.prior_precision()

        # For now just check size and format of result
        self.assertEqual(SPARSEFORMAT, Q.getformat())
        self.assertEqual((210,210), Q.shape)

    def test_prior_precision_derivative(self):

        prior = SeasonalElementPrior(
            SeasonalHyperparameters(n_spatial_components=3, common_log_sigma=0.8, common_log_rho=1.0),
            spde=SphereMeshSPDE(level=1, sparse_format=SPARSEFORMAT),
            n_harmonics=2, include_local_mean=True, alpha=2)

        dQ = prior.prior_precision_derivative(1)

        # For now just check size and format of result
        self.assertEqual(SPARSEFORMAT, dQ.getformat())
        self.assertEqual((210,210), dQ.shape)

    
