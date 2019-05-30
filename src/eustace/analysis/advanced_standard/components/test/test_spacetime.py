"""Test space-time component."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponentSolution
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.elements.element import Element
from eustace.analysis.advanced_standard.elements.element import ElementPrior
from eustace.analysis.advanced_standard.elements.element import ElementLinearDesign
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.element import  ObservationStructure

class TestSpaceTimeComponentSolution(unittest.TestCase):

    class TestDesign(ElementLinearDesign):

        def __init__(self, t):

            self.t = t

        def design_matrix(self):

	    return scipy.sparse.csc_matrix( [[0.0, 3.3 ]])
	    
        def design_jacobian(self, currentstate):

            if self.t == 21:
                
                return scipy.sparse.csc_matrix( [ [-1.5, 2.2 ] ] )

            elif self.t == 532:

                return scipy.sparse.csc_matrix( [ [ 0.0, 3.3 ] ] )

        def design_initial_state(self):

            return numpy.zeros((2,))

        def design_function(self, currentstate, rectified=False):
	    if rectified:
	      return numpy.abs(self.design_jacobian(currentstate)).dot(currentstate)
	    else:
	      return self.design_jacobian(currentstate).dot(currentstate)

    class TestElement(Element):

        def isnonlinear(self):

            return False

        def element_design(self, observationstructure):

            return TestSpaceTimeComponentSolution.TestDesign(observationstructure.time_index())

        def element_prior(self, hyperparameters):

            return CovariatePrior(hyperparameters, number_of_state_parameters=2)
        
    class TestObservations(ObservationStructure):

        def __init__(self, t):

            self.t = t

        def time_index(self):

            return self.t

        def observation_vector(self):

            if self.t == 21:

                return numpy.array( [ 7.0 ] )

            elif self.t == 532:

                return numpy.array( [ 9.0 ] )

	def number_of_observations(self):
	    return 1

        def observation_precision(self):

            return scipy.sparse.csc_matrix( [ [ 5.0 ] ] )

    def test_process_observations_no_uncertainties(self):

        # Our test system is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5  0.0 ] [ 5.0  0.0 ] [ -1.5  2.2 ] ) x = [ -1.5  0.0 ] [ 10.0  3.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [  2.2  3.3 ] [ 0.0  5.0 ] [  0.0  3.3 ] )     [  2.2  3.3 ] [  3.0 10.0 ] [ 9.0 - 3.0 ]
        #
        # [  13.25 -16.5  ] x = [ -37.5 ]
        # [-16.5   80.65 ]     [ 154.0 ]
        #
        # => x = [ -0.60697861 ]
        #        [  1.78530506 ]
        #

        c = SpaceTimeComponent(
            ComponentStorage_InMemory(TestSpaceTimeComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
            SpaceTimeComponentSolutionStorage_InMemory(), sample_size=666)
        s = c.component_solution()
        self.assertIsInstance(s, SpaceTimeComponentSolution)
        self.assertFalse(s.compute_uncertainties)
        test_offset = numpy.array([ 2.0 , 3.0 ])
        s.process_observations(TestSpaceTimeComponentSolution.TestObservations(t=21), test_offset[0:1])
        s.update_time_step()
        s.process_observations(TestSpaceTimeComponentSolution.TestObservations(t=532), test_offset[1:2])
        s.update_time_step()
        s.update()

        numpy.testing.assert_almost_equal(s.solutionstorage.state, numpy.array([ -0.60697861 ,  1.78530506 ]))

        # No marginal variances should have been computed at all, same for the sample
        self.assertEqual(None, s.solutionstorage.state_marginal_std)
        self.assertEqual(None, s.solutionstorage.state_sample)

        for time, expected_array in zip([21, 532], [numpy.array([ -1.5*-0.60697861 + 2.2*1.78530506 ]), numpy.array([ 3.3*1.78530506 ])]):
            # Observation at time t=t* should be design matrix for that time multiplied by expected state
            numpy.testing.assert_almost_equal(s.solution_observation_expected_value( TestSpaceTimeComponentSolution.TestObservations(t=time)), expected_array)
            
            # In this case we are considering a generical model solving iteration for the model, no marginal variances stored, hence we expect 0. as observation uncertainties            
            numpy.testing.assert_array_equal(s.solution_observation_expected_uncertainties(TestSpaceTimeComponentSolution.TestObservations(t=time)), 0.)
            # check samples are zero
            numpy.testing.assert_array_equal(0., s.solution_observation_projected_sample(TestSpaceTimeComponentSolution.TestObservations(t=time)))
            # check default number of samples
            self.assertEqual(s.solution_observation_projected_sample(TestSpaceTimeComponentSolution.TestObservations(t=time)).shape[1], 666)


    def test_process_observations_compute_uncertainties(self):

        # Our test system is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5  0.0 ] [ 5.0  0.0 ] [ -1.5  2.2 ] ) x = [ -1.5  0.0 ] [ 10.0  3.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [  2.2  3.3 ] [ 0.0  5.0 ] [  0.0  3.3 ] )     [  2.2  3.3 ] [  3.0 10.0 ] [ 9.0 - 3.0 ]
        #
        # [  13.25 -16.5  ] x = [ -37.5 ]
        # [-16.5   80.65 ]     [ 154.0 ]
        #
        # => x = [ -0.60697861 ]
        #        [  1.78530506 ]
        #

        c = SpaceTimeComponent(
            ComponentStorage_InMemory(TestSpaceTimeComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
            SpaceTimeComponentSolutionStorage_InMemory(), compute_uncertainties=True)
        s = c.component_solution()
        self.assertIsInstance(s, SpaceTimeComponentSolution)
        self.assertTrue(s.compute_uncertainties)
        test_offset = numpy.array([ 2.0 , 3.0 ])
        s.process_observations(TestSpaceTimeComponentSolution.TestObservations(t=21), test_offset[0:1])
        s.update_time_step()
        s.process_observations(TestSpaceTimeComponentSolution.TestObservations(t=532), test_offset[1:2])
        s.update_time_step()
        s.update()

        # No sample should have been computed at all
        self.assertEqual(None, s.solutionstorage.state_sample)


        # In this case we are considering the last iteration of model solving, hence marginal variances should have been stored
        expected_marginal_std = numpy.sqrt(numpy.diag(numpy.linalg.inv(numpy.array([[  13.25, -16.5  ], [-16.5,  80.65 ]]))))
        numpy.testing.assert_array_almost_equal(s.solutionstorage.state_marginal_std, expected_marginal_std)
        
        # Now we compute the projection of marginal variances onto the given observations
        for time in [532]:
            # Observation at time t=t* should be design matrix for that time multiplied by expected state
            expected_projection = TestSpaceTimeComponentSolution.TestDesign(t=time).design_function(expected_marginal_std)
            numpy.testing.assert_almost_equal(s.solution_observation_expected_uncertainties( TestSpaceTimeComponentSolution.TestObservations(t=time)), expected_projection)
            # check samples are zero
            numpy.testing.assert_array_equal(0., s.solution_observation_projected_sample(TestSpaceTimeComponentSolution.TestObservations(t=time)))
            # check default number of samples
            self.assertEqual(s.solution_observation_projected_sample(TestSpaceTimeComponentSolution.TestObservations(t=time)).shape[1], 1)

    def test_process_observations_compute_sample(self):

        # Our test system is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5  0.0 ] [ 5.0  0.0 ] [ -1.5  2.2 ] ) x = [ -1.5  0.0 ] [ 10.0  3.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [  2.2  3.3 ] [ 0.0  5.0 ] [  0.0  3.3 ] )     [  2.2  3.3 ] [  3.0 10.0 ] [ 9.0 - 3.0 ]
        #
        # [  13.25 -16.5  ] x = [ -37.5 ]
        # [-16.5   80.65 ]     [ 154.0 ]
        #
        # => x = [ -0.60697861 ]
        #        [  1.78530506 ]
        #
        
        number_of_samples = 200
        c = SpaceTimeComponent(
            ComponentStorage_InMemory(TestSpaceTimeComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
            SpaceTimeComponentSolutionStorage_InMemory(), compute_sample=True, sample_size=number_of_samples)
        s = c.component_solution()
        self.assertIsInstance(s, SpaceTimeComponentSolution)
        self.assertTrue(s.compute_sample)
        test_offset = numpy.array([ 2.0 , 3.0 ])
        s.process_observations(TestSpaceTimeComponentSolution.TestObservations(t=21), test_offset[0:1])
        s.update_time_step()
        s.process_observations(TestSpaceTimeComponentSolution.TestObservations(t=532), test_offset[1:2])
        s.update_time_step()
        numpy.random.seed(0)
        s.update()
        # In this case we are considering the last iteration of model solving, hence sample should have been stored
        computed_sample = s.solutionstorage.state_sample
        self.assertEqual((2, number_of_samples), computed_sample.shape)        	
        
        numpy.random.seed(0)
        variate = scipy.random.normal(0.0, 1.0, (2, number_of_samples))
        expected_posterior_precision = numpy.array([[  13.25, -16.5  ], [-16.5,   80.65 ]])
        # check result
        numpy.testing.assert_almost_equal(numpy.dot(variate.T,variate), numpy.dot(computed_sample.T,numpy.dot(expected_posterior_precision, computed_sample)))