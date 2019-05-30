"""Test spatial component."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spatial import SpatialComponentSolution
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.elements.element import Element
from eustace.analysis.advanced_standard.elements.element import ElementPrior
from eustace.analysis.advanced_standard.elements.element import ElementLinearDesign
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.element import  ObservationStructure

class TestSpatialComponentSolution(unittest.TestCase):

    class TestDesign(ElementLinearDesign):

        def design_matrix(self):

            return scipy.sparse.csc_matrix( [ [ -1.5 , 2.2 ], [ 0.0, 3.3 ] ] )

    class TestElement(Element):

        def element_design(self, observationstructure):

            return TestSpatialComponentSolution.TestDesign()

        def element_prior(self, hyperparameters):

            return CovariatePrior(hyperparameters, number_of_state_parameters=2)
        
    class TestObservations(ObservationStructure):

        def time_index(self):

            return 21

        def number_of_observations(self):

            return 2

        def observation_vector(self):

            return numpy.array( [ 7.0,  9.0 ] )

        def observation_precision(self):

            return scipy.sparse.csc_matrix( [ [ 10.0, 3.0 ], [ 3.0, 10.0 ] ] )


    class TestSomeOtherTime(ObservationStructure):

        def time_index(self):

            return 5927

        def number_of_observations(self):

            return 7

    def test_process_observations(self):

        # Our test system is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5  0.0 ] [ 10.0  3.0 ] [ -1.5  2.2 ] ) x = [ -1.5  0.0 ] [ 10.0  3.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [  2.2  3.3 ] [  3.0 10.0 ] [  0.0  3.3 ] )     [  2.2  3.3 ] [  3.0 10.0 ] [ 9.0 - 3.0 ]
        #
        # [  24.5   -47.85 ] x = [ -102.0 ]
        # [ -47.85  202.86 ]     [  397.1 ]
        #
        # => x = [ -0.63067268 ]
        #        [  1.80874649 ]
        #

        # Input data
        test_offset = numpy.array([ 2.0 , 3.0 ])
        test_obs = TestSpatialComponentSolution.TestObservations()

        # Make component and check it's of the correct class
        c = SpatialComponent(
            ComponentStorage_InMemory(TestSpatialComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
            SpatialComponentSolutionStorage_InMemory())
        s = c.component_solution()
        self.assertIsInstance(s, SpatialComponentSolution)
        self.assertFalse(s.compute_uncertainties)
        self.assertFalse(s.compute_sample)

        # Do the processing
        s.process_observations(test_obs, test_offset)
        s.update_time_step()

        # Should be the value of x (see matrices in comment above)
        numpy.testing.assert_almost_equal(s.solutionstorage.partial_state_read(21), numpy.array([ [ -0.63067268 ], [ 1.80874649 ] ]))

        # Should be the value of x multiplied by the design matrix
        numpy.testing.assert_almost_equal(
            s.solution_observation_expected_value(TestSpatialComponentSolution.TestObservations()),
            numpy.array([ [ -1.5*-0.63067268 + 2.2*1.80874649 ],  [ 3.3*1.80874649 ] ]))

        # This should be zero as we have no information for other times
        numpy.testing.assert_almost_equal(
            s.solution_observation_expected_value(TestSpatialComponentSolution.TestSomeOtherTime()),
            numpy.zeros((7,)))

	# In this case we did not allowed the computation of uncertainties or samples, hence the state of marginal variances should be None, while the expected uncertainties zero
	self.assertEqual(s.solutionstorage.partial_state_marginal_std_read(21), None)
	self.assertEqual(s.solutionstorage.state_marginal_std_at_time, {})
	numpy.testing.assert_almost_equal( s.solution_observation_expected_uncertainties(TestSpatialComponentSolution.TestObservations()),0.)

	self.assertEqual(s.solutionstorage.partial_state_sample_read(21), None)
	self.assertEqual(s.solutionstorage.state_sample_at_time, {})
	numpy.testing.assert_almost_equal( s.solution_observation_projected_sample(TestSpatialComponentSolution.TestObservations()),0.)


    def test_process_observations_compute_uncertainties(self):

        # Our test system is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5  0.0 ] [ 10.0  3.0 ] [ -1.5  2.2 ] ) x = [ -1.5  0.0 ] [ 10.0  3.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [  2.2  3.3 ] [  3.0 10.0 ] [  0.0  3.3 ] )     [  2.2  3.3 ] [  3.0 10.0 ] [ 9.0 - 3.0 ]
        #
        # [  24.5   -47.85 ] x = [ -102.0 ]
        # [ -47.85  202.86 ]     [  397.1 ]
        #
        # => x = [ -0.63067268 ]
        #        [  1.80874649 ]
        #

        # Input data
        test_offset = numpy.array([ 2.0 , 3.0 ])
        test_obs = TestSpatialComponentSolution.TestObservations()

        # Make component and check it's of the correct class
        c = SpatialComponent(
            ComponentStorage_InMemory(TestSpatialComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
            SpatialComponentSolutionStorage_InMemory(), compute_uncertainties = True, sample_size=666)
        s = c.component_solution()
        self.assertIsInstance(s, SpatialComponentSolution)
        self.assertTrue(s.compute_uncertainties)
        self.assertFalse(s.compute_sample)
        
        # Do the processing
        s.process_observations(test_obs, test_offset)
        s.update_time_step()

	
	# In this case we are considering the last iteration of model solving, hence marginal variances should have been stored
	expected_marginal_std = numpy.sqrt(numpy.diag(numpy.linalg.inv(numpy.array([[ 24.5, -47.85 ], [ -47.85,  202.86 ]]))))
	numpy.testing.assert_array_almost_equal(s.solutionstorage.partial_state_marginal_std_read(21), expected_marginal_std)
	self.assertEqual(len(s.solutionstorage.state_marginal_std_at_time),1)
	self.assertListEqual(s.solutionstorage.state_marginal_std_at_time.keys(), [21])
	
	# We also test the computation of prior marginal variances, and their projection onto observations
	
	expected_prior_std = numpy.sqrt(numpy.diag(numpy.linalg.inv(numpy.array([[ 2., 0. ], [ 0.,  2. ]]))))
	numpy.testing.assert_array_almost_equal(s.solution_prior_std(100), expected_prior_std, decimal=1)
	numpy.testing.assert_array_almost_equal(s.solution_prior_std(10000), expected_prior_std, decimal=2)
	
	design_matrix = numpy.array([[-1.5, 2.2], [0.0, 3.3]])
	expected_prior_std_projection = numpy.dot(design_matrix ,expected_prior_std)
	numpy.testing.assert_array_almost_equal(s.solution_observation_prior_uncertainties(None), expected_prior_std_projection, decimal=1)
	
	# In this case we did not allowed the computation of uncertainties samples
	self.assertEqual(s.solutionstorage.partial_state_sample_read(21), None)
	self.assertEqual(s.solutionstorage.state_sample_at_time, {})
	numpy.testing.assert_almost_equal( s.solution_observation_projected_sample(TestSpatialComponentSolution.TestObservations()),0.)
	self.assertEqual(666, s.solution_observation_projected_sample(TestSpatialComponentSolution.TestObservations()).shape[1])
	
    def test_process_observations_compute_sample(self):

        # Our test system is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5  0.0 ] [ 10.0  3.0 ] [ -1.5  2.2 ] ) x = [ -1.5  0.0 ] [ 10.0  3.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [  2.2  3.3 ] [  3.0 10.0 ] [  0.0  3.3 ] )     [  2.2  3.3 ] [  3.0 10.0 ] [ 9.0 - 3.0 ]
        #
        # [  24.5   -47.85 ] x = [ -102.0 ]
        # [ -47.85  202.86 ]     [  397.1 ]
        #
        # => x = [ -0.63067268 ]
        #        [  1.80874649 ]
        #

        # Input data
        test_offset = numpy.array([ 2.0 , 3.0 ])
        test_obs = TestSpatialComponentSolution.TestObservations()
	sample_size=300

        # Make component and check it's of the correct class
        c = SpatialComponent(
            ComponentStorage_InMemory(TestSpatialComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
            SpatialComponentSolutionStorage_InMemory(), compute_sample = True, sample_size=sample_size)
        s = c.component_solution()
        self.assertIsInstance(s, SpatialComponentSolution)
        self.assertFalse(s.compute_uncertainties)
        self.assertTrue(s.compute_sample)
 
	# Do the processing
        s.process_observations(test_obs, test_offset)
        numpy.random.seed(0)
        s.update_time_step()
	computed_sample = s.solutionstorage.partial_state_sample_read(21)
	self.assertEqual((2, sample_size), computed_sample.shape)        	

        numpy.random.seed(0)
	variate = scipy.random.normal(0.0, 1.0, (2, sample_size))

	expected_posterior_precision = numpy.array([[ 24.5, -47.85 ], [ -47.85,  202.86 ]])

	self.assertEqual(len(s.solutionstorage.state_sample_at_time),1)
	self.assertListEqual(s.solutionstorage.state_sample_at_time.keys(), [21])
        # check result
        numpy.testing.assert_almost_equal(numpy.dot(variate.T,variate), numpy.dot(computed_sample.T,numpy.dot(expected_posterior_precision, computed_sample)))

	# In this case we did not allowed the computation of uncertainties
	self.assertEqual(s.solutionstorage.partial_state_marginal_std_read(21), None)
	self.assertEqual(s.solutionstorage.state_marginal_std_at_time, {})
	numpy.testing.assert_almost_equal( s.solution_observation_expected_uncertainties(TestSpatialComponentSolution.TestObservations()),0.)
