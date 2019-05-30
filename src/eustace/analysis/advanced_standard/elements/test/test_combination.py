"""Test combination of analysis elements."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT

from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationDesign
from eustace.analysis.advanced_standard.elements.combination import CombinationPrior
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters

from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalDesign
from eustace.analysis.advanced_standard.elements.local import LocalPrior

from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanElement
from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanDesign
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.element import ObservationStructure

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE

class TestCombinationHyperparameters(unittest.TestCase):

    def test_init(self):

        p = CombinationHyperparameters( [ CovariateHyperparameters(23.6), CovariateHyperparameters(22.9) ] )
        self.assertEqual(2, len(p.elementparameters))
        self.assertEqual(23.6, p.elementparameters[0].value)
        self.assertEqual(22.9, p.elementparameters[1].value)

    def test_get_array(self):

        p = CombinationHyperparameters( [ CovariateHyperparameters(23.6), LocalHyperparameters(log_sigma=0.1, log_rho=1.2) ] )
        numpy.testing.assert_equal([ 23.6, 0.1, 1.2 ], p.get_array())

    def test_set_array(self):

        p = CombinationHyperparameters( [ CovariateHyperparameters(23.6), LocalHyperparameters(log_sigma=0.1, log_rho=1.2), CovariateHyperparameters(24.7) ] )
        p.set_array(numpy.array([1.4, 1.5, 1.6, 1.7]))
        self.assertEqual(3, len(p.elementparameters))
        self.assertEqual(1.4, p.elementparameters[0].value)
        self.assertEqual(1.5, p.elementparameters[1].log_sigma)
        self.assertEqual(1.6, p.elementparameters[1].log_rho)
        self.assertEqual(1.7, p.elementparameters[2].value)
        numpy.testing.assert_equal([1.4, 1.5, 1.6, 1.7], p.get_array())

    def test_get_element_ranges(self):

        p = CombinationHyperparameters( [ CovariateHyperparameters(23.6), LocalHyperparameters(log_sigma=0.1, log_rho=1.2), CovariateHyperparameters(24.7) ] )
        self.assertEqual([ [ 0 ], [ 1, 2 ], [ 3 ] ], p.get_element_ranges())

class TestCombinationElement(unittest.TestCase):

    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 15.0, -7.0 ], [  5.0, 100.0 ] ])                 

    def test_init(self):

        c = CombinationElement([ GrandMeanElement(), LocalElement(0) ])
        self.assertEqual(2, len(c.combination))
        self.assertTrue(isinstance(c.combination[0], GrandMeanElement))
        self.assertTrue(isinstance(c.combination[1], LocalElement))

    def test_element_design(self):
        
        obs = TestCombinationElement.SimulatedObservationStructure()
        design = CombinationElement([ GrandMeanElement(), LocalElement(0) ]).element_design(obs)
        self.assertTrue(isinstance(design, CombinationDesign))
        self.assertEqual(2, len(design.designlist))
        self.assertTrue(isinstance(design.designlist[0], GrandMeanDesign))
        self.assertTrue(isinstance(design.designlist[1], LocalDesign))

    def test_element_prior(self):

        hyperparameters = CombinationHyperparameters([CovariateHyperparameters(23.6), LocalHyperparameters(log_sigma=0.1, log_rho=1.2) ])
        prior = CombinationElement([ GrandMeanElement(), LocalElement(0) ]).element_prior(hyperparameters)
        self.assertTrue(isinstance(prior, CombinationPrior))
        self.assertEqual(2, len(prior.priorlist))
        self.assertTrue(isinstance(prior.priorlist[0], CovariatePrior))
        self.assertTrue(isinstance(prior.priorlist[1], LocalPrior))

class TestCombinationDesign(unittest.TestCase):

    def test_init(self):

        obs = TestCombinationElement.SimulatedObservationStructure()
        spde = SphereMeshSPDE(0)
        design = CombinationDesign([ GrandMeanDesign(2), LocalDesign(obs, spde) ])
        self.assertEqual(2, len(design.designlist))
        self.assertTrue(isinstance(design.designlist[0], GrandMeanDesign))
        self.assertTrue(isinstance(design.designlist[1], LocalDesign))

    def test_element_states(self):

        obs = TestCombinationElement.SimulatedObservationStructure()
        spde = SphereMeshSPDE(0)
        design = CombinationDesign([ GrandMeanDesign(obs), LocalDesign(obs, spde) ])
        states = design.element_states(numpy.array( [ 270.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 ] ))
        self.assertTrue(isinstance(states, list))
        self.assertEqual(2, len(states))
        numpy.testing.assert_equal(states[0], [ 270.0 ])
        numpy.testing.assert_equal(states[1], [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0 ])
    
    def test_design_function(self):

        obs = TestCombinationElement.SimulatedObservationStructure()
        spde = SphereMeshSPDE(1)
        design = CombinationDesign([ GrandMeanDesign(2), LocalDesign(obs, spde) ])

        A = spde.build_A(obs.location_polar_coordinates())
        observation_indices, vertex_indices = A.sorted_indices().nonzero()

        # Make a state vector full of a test value 561 which should never actually get accessed
        state_vector = numpy.tile(561.0, (47,1))

        # Set the grand mean to 20.0
        state_vector[0, 0] = 20.0

        # Choose some numbers for observations of interest
        state_vector[(vertex_indices[0:3] + 1),0] = 43.5
        state_vector[(vertex_indices[3:6] + 1),0] = -27.0

        # Should end up with (20.0 + 43.5) and (20.0 - 27.0)
        y = design.design_function(state_vector)
        numpy.testing.assert_almost_equal(y, [ [ 63.5 ], [ -7.0 ] ])

    def test_design_jacobian(self):

        obs = TestCombinationElement.SimulatedObservationStructure()
        spde = SphereMeshSPDE(1)
        design = CombinationDesign([ GrandMeanDesign(2), LocalDesign(obs, spde) ])

        # Make a state vector full of a test value 561 which should never actually get accessed
        state_vector = numpy.tile(561.0, (47,1))

        J = design.design_jacobian(state_vector)

        # A-matrix for comparison
        A = spde.build_A(obs.location_polar_coordinates())

        # Should have ones at start for global mean
        J_expected = numpy.hstack([ numpy.array([ [ 1.0 ], [ 1.0 ] ]), A.todense() ] )

        # Check as exepcted
        self.assertEqual(SPARSEFORMAT, J.getformat())
        numpy.testing.assert_equal(J_expected, J.todense())

class TestCombinationPrior(unittest.TestCase):

    def test_prior_number_of_state_parameters(self):

        h = CombinationHyperparameters([CovariateHyperparameters(0.1), CovariateHyperparameters(0.2)])
        priorlist = [ 
            CovariatePrior(h.elementparameters[0], number_of_state_parameters=4),
            CovariatePrior(h.elementparameters[1], number_of_state_parameters=5) ]
        prior = CombinationPrior(h, priorlist)
        self.assertEqual(9, prior.prior_number_of_state_parameters())

    def test_prior_precision(self):

        h = CombinationHyperparameters([CovariateHyperparameters(-0.5 * numpy.log(23.6)), CovariateHyperparameters(-0.5*numpy.log(88.9))])
        priorlist = [ 
            CovariatePrior(h.elementparameters[0], number_of_state_parameters=2),
            CovariatePrior(h.elementparameters[1], number_of_state_parameters=3) ]
        prior = CombinationPrior(h, priorlist)
        Q = prior.prior_precision()
        self.assertEqual(SPARSEFORMAT, Q.getformat())
        numpy.testing.assert_almost_equal(
            Q.todense(),
            [ [ 23.6, 0.0, 0.0, 0.0, 0.0 ],
              [  0.0,23.6, 0.0, 0.0, 0.0 ],
              [  0.0, 0.0,88.9, 0.0, 0.0 ],
              [  0.0, 0.0, 0.0,88.9, 0.0 ],
              [  0.0, 0.0, 0.0, 0.0,88.9 ] ])
        
    def test_prior_precision_derivative(self):

        h = CombinationHyperparameters([CovariateHyperparameters(-0.5 * numpy.log(23.6)), CovariateHyperparameters(-0.5*numpy.log(88.9))])
        priorlist = [ 
            CovariatePrior(h.elementparameters[0], number_of_state_parameters=2),
            CovariatePrior(h.elementparameters[1], number_of_state_parameters=3) ]
        prior = CombinationPrior(h, priorlist)
        dQ = prior.prior_precision_derivative(1)
        self.assertEqual(SPARSEFORMAT, dQ.getformat())
        self.assertEqual((5,5), dQ.shape)
        self.assertEqual(3, dQ.nnz)
        self.assertAlmostEqual(-2.0 * 88.9, dQ[2,2])
        self.assertAlmostEqual(-2.0 * 88.9, dQ[3,3])
        self.assertAlmostEqual(-2.0 * 88.9, dQ[4,4])
