"""Tests for grand mean element."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanElement
from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanDesign
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT

class TestGrandMeanElement(unittest.TestCase):
    
    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 15.0, -7.0 ], [  5.0, 100.0 ] ])                 
    
    def test_element_design(self):

        grandmean = GrandMeanElement()
        design = grandmean.element_design(TestGrandMeanElement.SimulatedObservationStructure())
        self.assertTrue(isinstance(design, GrandMeanDesign))
        self.assertEqual(2, design.number_of_observations)

    def test_element_prior(self):

        prior = GrandMeanElement().element_prior(CovariateHyperparameters(3.3333))
        self.assertTrue(isinstance(prior, CovariatePrior))
        self.assertEqual(1, prior.number_of_state_parameters)
        self.assertEqual(3.3333, prior.hyperparameters.value)


class TestGrandMeanDesign(unittest.TestCase):
               
    def test_design_number_of_state_parameters(self):

        self.assertEqual(1, GrandMeanDesign(number_of_observations=3).design_number_of_state_parameters())

    def test_isnonlinear(self):
        
        self.assertFalse(GrandMeanDesign(number_of_observations=3).isnonlinear())

    def test_design_matrix(self):

        A = GrandMeanDesign(number_of_observations=7).design_matrix()
        self.assertEqual(SPARSEFORMAT, A.getformat())
        numpy.testing.assert_almost_equal(A.todense(), [ [ 1.0 ], [ 1.0 ], [ 1.0 ], [ 1.0 ], [ 1.0 ], [ 1.0 ], [ 1.0 ] ])

    def test_design_jacobian(self):

        J = GrandMeanDesign(number_of_observations=6).design_jacobian(currentstate=None)
        self.assertEqual(SPARSEFORMAT, J.getformat())
        numpy.testing.assert_almost_equal(J.todense(), [ [ 1.0 ], [ 1.0 ], [ 1.0 ], [ 1.0 ], [ 1.0 ], [ 1.0 ] ])

    def test_design_function(self):
        
        y = GrandMeanDesign(number_of_observations=3).design_function(currentstate=numpy.array([ [ 27.9 ] ]))
        numpy.testing.assert_almost_equal(y, [ [ 27.9 ], [ 27.9 ], [ 27.9 ] ])
    
