"""Tests for bias element."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.elements.bias import BiasDesign
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT

class TestBiasElement(unittest.TestCase):

    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 23.0, 3.2 ], [ 15.0, -7.0 ], [  5.0, 100.0 ] ])

        def covariate_effect(self, groupname):

            if groupname == 'Bob':

                # obs 1 --> bias 0
                # obs 2 --> bias 4
                return numpy.array( [ [ 1, 0 ],
                                      [ 2, 4 ] ], numpy.int64 )

            elif groupname == 'NotBob':

                # obs 0 --> bias 2
                return numpy.array( [ [ 0, 2 ] ], numpy.int64 )

            else:

                # (empty map)
                return numpy.array( [ [ ] ], numpy.int64)

        def number_of_observations(self):

            return self.location_polar_coordinates().shape[0]

    def test_init(self):

        b = BiasElement('Bob', 42)
        self.assertEqual('Bob', b.groupname)
        self.assertEqual(42, b.number_of_biases)


    def test_element_design(self):

        design = BiasElement('Bob', 42).element_design(TestBiasElement.SimulatedObservationStructure())

        # Check type
        self.assertTrue(isinstance(design, BiasDesign))

        # Covariate info should be copied from constructor
        self.assertEqual('Bob', design.groupname)
        self.assertEqual(42, design.number_of_biases)

        # And effect info comes from observation structure
        numpy.testing.assert_equal(design.effect, [ [ 1, 0 ], [ 2, 4 ] ])


    def test_element_prior(self):

        prior = BiasElement('Bob', 42).element_prior(CovariateHyperparameters(9.9))

        self.assertTrue(isinstance(prior, CovariatePrior))
        self.assertEqual(9.9, prior.hyperparameters.value)

        # The number of state parameters should be the total number of biases (irrespective of number observed)
        self.assertEqual(42, prior.number_of_state_parameters)



class TestBiasDesign(unittest.TestCase):

    def test_init(self):

        design_bob = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Bob', 5)
        self.assertEqual('Bob', design_bob.groupname)
        self.assertEqual(5, design_bob.number_of_biases)
        self.assertEqual(3, design_bob.number_of_observations)
        numpy.testing.assert_equal(design_bob.effect, [ [ 1, 0 ], [ 2, 4 ] ])

        design_notbob = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'NotBob', 6)
        self.assertEqual('NotBob', design_notbob.groupname)
        self.assertEqual(6, design_notbob.number_of_biases)
        self.assertEqual(3, design_notbob.number_of_observations)
        numpy.testing.assert_equal(design_notbob.effect, [ [ 0, 2 ] ])

        design_nobody = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Nobody', 7)
        self.assertEqual('Nobody', design_nobody.groupname)
        self.assertEqual(7, design_nobody.number_of_biases)
        self.assertEqual(3, design_nobody.number_of_observations)
        numpy.testing.assert_equal(design_nobody.effect, [ [ ] ])

    def test_design_number_of_state_parameters(self):

        design = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Bob', 5)
        self.assertEqual(5, design.design_number_of_state_parameters())

    def test_design_matrix(self):

        A = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Bob', 5).design_matrix()
        self.assertEqual(SPARSEFORMAT, A.getformat())
        self.assertEqual((3, 5), A.shape)
        self.assertEqual(2, A.nnz)

        # This mapping comes from observation definition
        # in TestBiasElement.SimulatedObservationStructure:
        # obs 1 --> bias 0, obs 2 --> bias 4
        numpy.testing.assert_equal(A.todense(), 
                                   [ [ 0.0, 0.0, 0.0, 0.0, 0.0 ],
                                     [ 1.0, 0.0, 0.0, 0.0, 0.0 ],
                                     [ 0.0, 0.0, 0.0, 0.0, 1.0 ] ])

    def test_design_function(self):

        design = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Bob', 5)
        y = design.design_function(numpy.array([ [ 22.2 ], [ 33.3 ], [ 44.4 ], [ 55.5 ], [ 66.6 ] ]))

        # Again this comes from observation definition
        # in TestBiasElement.SimulatedObservationStructure:
        # obs 0 has no associated bias (so is zero)
        # obs 1 --> bias 0 (which is 22.2)
        # obs 2 --> bias 4 (which is 66.6)
        numpy.testing.assert_almost_equal(y, numpy.array([ [ 0.0 ], [ 22.2 ], [ 66.6 ] ]))

    def test_isnonlinear(self):
        
        self.assertFalse(BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Bob', 5).isnonlinear())

    def test_design_jacobian(self):

        A = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Bob', 5).design_matrix()
        J = BiasDesign(TestBiasElement.SimulatedObservationStructure(), 'Bob', 5).design_jacobian(numpy.array([ ]))
        self.assertEqual(SPARSEFORMAT, J.getformat())
        numpy.testing.assert_almost_equal(J.todense(), A.todense())
