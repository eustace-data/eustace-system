"""Test delayed solve spatial component."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.components.spatialdelayed import DelayedSpatialComponent
from eustace.analysis.advanced_standard.components.spatialdelayed import DelayedSpatialComponentSolution
#from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
#from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_files import DelayedSpatialComponentSolutionStorage_Files
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory

from eustace.analysis.advanced_standard.elements.element import Element
from eustace.analysis.advanced_standard.elements.element import ElementPrior
from eustace.analysis.advanced_standard.elements.element import ElementLinearDesign
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.element import  ObservationStructure

class TestDelayedSpatialComponentSolution(unittest.TestCase):

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

            return TestDelayedSpatialComponentSolution.TestDesign(observationstructure.time_index())

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

        # Our test system for the first time step (key 21) is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5 ] [ 5.0 ] [ -1.5 2.2 ] ) x = [ -1.5 ] [ 5.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [ 2.2  ]                      )     [  2.2 ]         
        #
        # [ 13.25 -16.5  ] x = [ -37.5 ]
        # [-16.5   26.2 ]      [ 55.0  ]
        #
        # => x = [-1.00133511 ]
        #        [ 1.46862483 ]

        # Our test system for the first time step (key 21) is:
        #
        # ( [ 2.0  0.0 ]  +  [ 0.0 ] [ 5.0 ] [ 0.0 3.3 ] ) x = [  0.0 ] [ 5.0 ] [ 9.0 - 3.0 ]
        # ( [ 0.0  2.0 ]     [ 3.3  ]                    )     [  3.3 ]         
        #
        # [ 2.0   0.0  ] x = [ 0.0  ]
        # [ 0.0   56.45 ]     [ 99.0 ]
        #
        # => x = [ 0.         ]
        #        [ 1.75376439 ]
        
        for component_storage_class in DelayedSpatialComponentSolutionStorage_Files, SpatialComponentSolutionStorage_InMemory:

            c = DelayedSpatialComponent(
                ComponentStorage_InMemory(TestDelayedSpatialComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
                component_storage_class())
            c.solutionstorage.statefiledictionary_read = None
            c.solutionstorage.statefiledictionary_write = {21:'state_test.A.pickle', 532:'state_test.B.pickle'}
            c.solutionstorage.measurementfiledictionary_write = c.solutionstorage.measurementfiledictionary_read = {21: 'measurement_test.A.pickle', 532: 'measurement_test.B.pickle'}
            
            s = c.component_solution()
            self.assertIsInstance(s, DelayedSpatialComponentSolution)
            self.assertFalse(s.compute_uncertainties)
            test_offset = numpy.array([ 2.0 , 3.0 ])
            c.solutionstorage.state_time_index = 21
            c.solutionstorage.measurement_time_index_write = 21
            s.process_observations(TestDelayedSpatialComponentSolution.TestObservations(t=21), test_offset[0:1])
            s.update_time_step()
            c.solutionstorage.state_time_index = 532
            c.solutionstorage.measurement_time_index_write = 532
            s.process_observations(TestDelayedSpatialComponentSolution.TestObservations(t=532), test_offset[1:2])
            s.update_time_step()
            s.update()

            c.solutionstorage.statefiledictionary_read = c.solutionstorage.statefiledictionary_write    # now enable reading from the previously written files
            
            numpy.testing.assert_almost_equal(s.solutionstorage.partial_state_read( 21  ), numpy.array([ -1.00133511 ,  1.46862483 ]))
            numpy.testing.assert_almost_equal(s.solutionstorage.partial_state_read( 532 ), numpy.array([ 0.0         ,  1.75376439 ]))
            
            # No marginal variances should have been computed at all
            self.assertEqual(None, s.solutionstorage.partial_state_marginal_std_read(21))
            
            for time, expected_array in zip([21, 532], [numpy.array([ -1.5*-1.00133511 + 2.2*1.46862483 ]), numpy.array([ 3.3*1.75376439 ])]):
                # Observation at time t=t* should be design matrix for that time multiplied by expected state
                numpy.testing.assert_almost_equal(s.solution_observation_expected_value( TestDelayedSpatialComponentSolution.TestObservations(t=time)), expected_array)
                
                # In this case we are considering a generical model solving iteration for the model, no marginal variances stored, hence we expect 0. as observation uncertainties            
                numpy.testing.assert_array_equal(s.solution_observation_expected_uncertainties(TestDelayedSpatialComponentSolution.TestObservations(t=time)), 0.)
    
    def atest_process_observations_compute_uncertainties(self):

        # Our test system is:
        #
        # ( [ 2.0  0.0 ]  +  [ -1.5  0.0 ] [ 5.0  0.0 ] [ -1.5  2.2 ] ) x = [ -1.5  0.0 ] [ 10.0  3.0 ] [ 7.0 - 2.0 ]
        # ( [ 0.0  2.0 ]     [  2.2  3.3 ] [ 0.0  5.0 ] [  0.0  3.3 ] )     [  2.2  3.3 ] [  3.0 10.0 ] [ 9.0 - 3.0 ]
        #
        # [  3.25 -16.5  ] x = [ -37.5 ]
        # [-16.5   80.65 ]     [ 154.0 ]
        #
        # => x = [ -0.60697861 ]
        #        [  1.78530506 ]
        #

        c = DelayedSpatialComponent(
            ComponentStorage_InMemory(TestDelayedSpatialComponentSolution.TestElement(), CovariateHyperparameters(-0.5*numpy.log(2.0))),
            DelayedSpatialComponentSolutionStorage_Files(), compute_uncertainties=True)
        s = c.component_solution()
        self.assertIsInstance(s, DelayedSpatialComponentSolution)
        self.assertTrue(s.compute_uncertainties)
        test_offset = numpy.array([ 2.0 , 3.0 ])
        s.process_observations(TestDelayedSpatialComponentSolution.TestObservations(t=21), test_offset[0:1])
        s.update_time_step()
        s.process_observations(TestDelayedSpatialComponentSolution.TestObservations(t=532), test_offset[1:2])
        s.update_time_step()
        s.update()

        # In this case we are considering the last iteration of model solving, hence marginal variances should have been stored
        expected_marginal_std = numpy.sqrt(numpy.diag(numpy.linalg.inv(numpy.array([[  13.25, -16.5  ], [-16.5,  80.65 ]]))))
        numpy.testing.assert_array_almost_equal(s.solutionstorage.state_marginal_std, expected_marginal_std)
        
        # Now we compute the projection of marginal variances onto the given observations
        for time in [532]:
            # Observation at time t=t* should be design matrix for that time multiplied by expected state
            expected_projection = TestDelayedSpatialComponentSolution.TestDesign(t=time).design_function(expected_marginal_std)
            numpy.testing.assert_almost_equal(s.solution_observation_expected_uncertainties( TestDelayedSpatialComponentSolution.TestObservations(t=time)), expected_projection)
