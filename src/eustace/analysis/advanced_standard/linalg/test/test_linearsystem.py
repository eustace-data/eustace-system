import unittest
import numpy
import scipy.sparse
from eustace.analysis.advanced_standard.linalg.linearsystem import MeasurementUpdate
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystem
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_Cholesky
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_CG
from eustace.analysis.advanced_standard.elements.element import ElementPrior
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE

class TestMeasurementUpdate(unittest.TestCase):

    def test_create_blank(self):

        m = MeasurementUpdate.create_blank(dimension=7)
        self.assertIsInstance(m, MeasurementUpdate)
        numpy.testing.assert_equal(numpy.zeros((7,1)), m.information_vector.toarray())
        numpy.testing.assert_equal(numpy.zeros((7,7)), m.information_precision.toarray())

    def test_create_from_system_observations(self):
        """Test observations of an example system where:
               [ y0 ] = [ -1.5  2.2 ] [ x0 ]   +  noise(precision=diag(2.0, 4.0, 3.0))
               [ y1 ] = [  0.0  3.3 ] [ x1 ]
               [ y2 ] = [  1.0  0.0 ]
        """

        design_jacobian = scipy.sparse.csc_matrix( [ [ -1.5 , 2.2 ], [ 0.0, 3.3 ], [ 1.0, 0.0 ] ] )
        observation_precision = scipy.sparse.csc_matrix( [ 
                [ 2.0 , 0.0, 0.0 ],
                [ 0.0 , 4.0, 0.0 ],
                [ 0.0 , 0.0, 3.0 ] ] )
        relative_observations = numpy.array([ 10.0, 20.0, 50.0 ])
        result = MeasurementUpdate.create_from_system_observations(design_jacobian, observation_precision, relative_observations)
        self.assertIsInstance(result, MeasurementUpdate)
        numpy.testing.assert_almost_equal([ [ 120.0 ], [ 308.0 ] ], result.information_vector)
        numpy.testing.assert_almost_equal([ [ 7.5 , -6.6 ], [ -6.6 , 53.24 ] ], result.information_precision.toarray())
        

    def test_append_measurement_update(self):

        base = MeasurementUpdate.create_blank(dimension=5)
        
        update_a = MeasurementUpdate(
            scipy.sparse.csc_matrix( [ [ 4.0, 0.5, 0.0, 0.0, 0.0 ],
                                       [ 0.5, 4.0, 0.0, 0.0, 0.0 ],
                                       [ 0.0, 0.0, 3.0, 0.0, 0.0 ],
                                       [ 0.0, 0.0, 0.0, 3.0, 0.0 ],
                                       [ 0.0, 0.0, 0.0, 0.0, 5.0 ] ] ),
            scipy.sparse.csc_matrix( [ [ 0.0 ], [ 3.5 ], [ 0.0 ], [ 1.2 ], [ 0.0 ] ] ))

        update_b = MeasurementUpdate(
            scipy.sparse.csc_matrix( [ [ 1.0, 0.0, 0.0, 0.0, 0.0 ],
                                       [ 0.0, 2.0, 0.0, 0.0, 0.0 ],
                                       [ 0.0, 0.0, 1.0, 0.0, 0.0 ],
                                       [ 0.0, 0.0, 0.0, 2.0, 0.0 ],
                                       [ 0.0, 0.0, 0.0, 0.0, 1.0 ] ] ),
            scipy.sparse.csc_matrix( [ [ 0.0 ], [ 0.2 ], [ 0.0 ], [ 0.0 ], [ 0.1 ] ] ))

        base.append_measurement_update(update_a)
        base.append_measurement_update(update_b)

        numpy.testing.assert_almost_equal(            
            base.information_precision.toarray(),
            [ [ 5.0, 0.5, 0.0, 0.0, 0.0 ],
              [ 0.5, 6.0, 0.0, 0.0, 0.0 ],
              [ 0.0, 0.0, 4.0, 0.0, 0.0 ],
              [ 0.0, 0.0, 0.0, 5.0, 0.0 ],
              [ 0.0, 0.0, 0.0, 0.0, 6.0 ] ])

        numpy.testing.assert_almost_equal(
            base.information_vector.toarray(),
            [ [ 0.0 ], [ 3.7 ], [ 0.0 ], [ 1.2 ], [ 0.1 ] ])

        # Check that sparsity structure was recognised and preserved
        # (should know the number of non-zero entries if so)
        self.assertEqual(7, base.information_precision.nnz)
        self.assertEqual(3, base.information_vector.nnz)



class ExamplePrior(ElementPrior):
    """Mock prior to use as example."""

    def prior_number_of_state_parameters(self):

        return 3

    def prior_precision(self):

        return scipy.sparse.diags( [ 2.2, 3.3, 4.4 ], format=SPARSEFORMAT)

    def prior_precision_derivative(self, parameter_index):

        return scipy.sparse.diags( [ 1.0 if parameter_index == index else 0.0 for index in range(3) ] )


class TestLinearSystem(unittest.TestCase):

    def test_init(self):

        s = LinearSystem(ExamplePrior())
        self.assertTrue(isinstance(s.prior, ExamplePrior))
        self.assertTrue(isinstance(s.measurement, MeasurementUpdate))
        self.assertEqual((3,3), s.measurement.information_precision.shape)
        self.assertEqual((3,1), s.measurement.information_vector.shape)

    
    def test_append_measurement_update(self):

        s = LinearSystem(ExamplePrior())
        
        inputupdate = MeasurementUpdate(
            scipy.sparse.csc_matrix( [ [ 4.0, 0.5, 0.0 ],
                                       [ 0.5, 4.0, 0.0 ],
                                       [ 0.0, 0.0, 3.0 ] ] ),
            scipy.sparse.csc_matrix( [ [ 0.0 ], [ 3.5 ], [ 0.0 ] ] ))

        s.append_measurement_update(inputupdate)

        numpy.testing.assert_almost_equal(            
            s.measurement.information_precision.toarray(),
            [ [ 4.0, 0.5, 0.0 ],
              [ 0.5, 4.0, 0.0 ],
              [ 0.0, 0.0, 3.0 ] ])

        numpy.testing.assert_almost_equal(
            s.measurement.information_vector.toarray(),
            [ [ 0.0 ], [ 3.5 ], [ 0.0 ] ])

        # Check that sparsity structure was recognised and preserved
        # (should know the number of non-zero entries if so)
        self.assertEqual(5, s.measurement.information_precision.nnz)
        self.assertEqual(1, s.measurement.information_vector.nnz)
        
class TestLinearSystemSolution_Cholesky(unittest.TestCase):

    def setUp(self):

        # empty example linear system (will have dimension 3 due to dimension of prior)
        self.simulated_system = LinearSystem(ExamplePrior())
        
        # example measurement precision
        simulated_measurement_precision = scipy.sparse.csc_matrix( [ [ 4.0, 1.0, 0.0 ],
                                                                     [ 1.0, 4.0, 0.5 ],
                                                                     [ 0.0, 0.5, 3.0 ] ] )

        # the expected posterior result is the above plus ExamplePrior.prior_precision
        self.expected_posterior = numpy.array([ [ 6.2, 1.0, 0.0 ],
                                                [ 1.0, 7.3, 0.5 ],
                                                [ 0.0, 0.5, 7.4 ] ])
                                                
	# the expected marginal variance vector is given by the diagonal elements of posterior precision
	self.expected_marginal_variance = numpy.diagonal(numpy.linalg.inv(self.expected_posterior))

        # simulate solution of mean
        self.expected_mean = numpy.array([ [ 3.0 ] , [ 4.0 ] , [ 5.0 ] ])
        
        # and the observation vector that would produce it
        simulated_observations = numpy.dot(self.expected_posterior, self.expected_mean)

        # make the update object with this info...
        inputupdate = MeasurementUpdate(
            scipy.sparse.csc_matrix(simulated_measurement_precision),
            scipy.sparse.csc_matrix(simulated_observations))

        # ...and append to system
        self.simulated_system.append_measurement_update(inputupdate)

    def test_init(self):

        # Just check that initialisation runs (numerical tests are later)
        solution = LinearSystemSolution_Cholesky(self.simulated_system)
        self.assertEqual(3, solution.number_of_state_parameters)

    def test_maximum_a_posteriori(self):

        solution = LinearSystemSolution_Cholesky(self.simulated_system)
        numpy.testing.assert_almost_equal(solution.maximum_a_posteriori().squeeze(), [ 3.0 , 4.0 , 5.0 ])

    def test_marginal_std(self):
        solution = LinearSystemSolution_Cholesky(self.simulated_system)
        numpy.testing.assert_almost_equal(solution.marginal_std(), numpy.sqrt(self.expected_marginal_variance))

    def test_approximated_marginal_std(self):
	solution = LinearSystemSolution_Cholesky(self.simulated_system)
        numpy.testing.assert_almost_equal(solution.marginal_std(), solution.approximated_marginal_std(100), decimal=1)
        numpy.testing.assert_almost_equal(solution.marginal_std(), solution.approximated_marginal_std(10000), decimal=2)
        
    def test_white_noise_to_posterior_distribution(self):

        sample = numpy.array([ [ 88.0, 99.0, 34.1 ] ]).T
        solution = LinearSystemSolution_Cholesky(self.simulated_system)
        result = solution.white_noise_to_posterior_distribution(sample)

        # check that dot product is satisfied
        self.assertAlmostEqual(numpy.dot(sample.T, sample), numpy.dot(result.T, numpy.dot(self.expected_posterior, result)))

    def test_random_sample(self):

        # Make the solution
        solution = LinearSystemSolution_Cholesky(self.simulated_system)

        # Initialise random seed
        numpy.random.seed(987654)
        sample = scipy.random.normal(0.0, 1.0, (3, 200))

        # Re-initialise random seed so that random numbers
        # used by solution are same as above
        numpy.random.seed(987654)
        
        # Do random sample
        result = solution.random_sample(200)

        # check shape
        self.assertEqual((3, 200), result.shape)

        # check result
        numpy.testing.assert_almost_equal(numpy.dot(sample.T, sample), numpy.dot(result.T, numpy.dot(self.expected_posterior, result)))
        

class SimulatedLinearSystemOnMesh(LinearSystem):

    def __init__(self, level_number, number_of_observations):

        # Mesh
        self.spheremesh = SphereMeshSPDE(level=level_number)
        
        # Some random set of points in polar coordinates
        numpy.random.seed(12345)
        locations_polar = numpy.vstack( [ 180.0*numpy.random.rand(1,number_of_observations) -  90.0,
                                          360.0*numpy.random.rand(1,number_of_observations) - 180.0 ] )

        # Get design matrix
        A = self.spheremesh.build_A(locations_polar.T)

        # Assume noise level 1 and independent ==> Qe = identity
        self.measurement = MeasurementUpdate(
            information_precision=A.T.dot(A),
            information_vector=A.T.dot(10.0*numpy.random.randn(number_of_observations)))

        # This means self will get consulted to define the prior
        self.prior = self

    def prior_precision(self):

        # return scipy.sparse.eye(self.prior_number_of_state_parameters(), format='csc')
        return self.spheremesh.build_Q_stationary(log_sigma=0.0, log_rho=numpy.log(1.0), alpha=2)

    def prior_number_of_state_parameters(self):

        return self.measurement.information_precision.shape[0]

class TestLinearSystemSolutionsOnMesh(unittest.TestCase):
    """Checks using mesh as example."""

    def test_cholesky_solution(self):

        # Make simulated system
        simsystem = SimulatedLinearSystemOnMesh(5, 10000)

        # Solve it and check solution makes sense
        solution = LinearSystemSolution_Cholesky(simsystem)

        # This should be the matrix of the equation that was solved
        M = simsystem.prior.prior_precision() + simsystem.measurement.information_precision

        # Multiplying by state vector should give the information vector back again
        statevector = solution.maximum_a_posteriori().ravel()
        numpy.testing.assert_almost_equal(M.dot(statevector), simsystem.measurement.information_vector)

    def test_cholesky_inverse(self):
      
	# Make simulated system
        simsystem = SimulatedLinearSystemOnMesh(3, 1000)

        # Solve it and check solution makes sense
        solution = LinearSystemSolution_Cholesky(simsystem)

        # This should be the matrix of the equation that was solved
        M = simsystem.prior.prior_precision() + simsystem.measurement.information_precision

        # This should be the vector containing just the marginal viariances
        diag_inv_M = numpy.diag(numpy.linalg.inv(M.todense()))
        numpy.testing.assert_almost_equal(numpy.sqrt(diag_inv_M), solution.marginal_std())
        
    def test_cg(self):

        # Make simulated system
        simsystem = SimulatedLinearSystemOnMesh(5, 10000)

        # Solve it and check solution makes sense
        solution = LinearSystemSolution_CG(simsystem)

        # This should be the matrix of the equation that was solved
        M = simsystem.prior.prior_precision() + simsystem.measurement.information_precision

        # Multiplying by state vector should give the information vector back again
        statevector = solution.maximum_a_posteriori().ravel()
        numpy.testing.assert_almost_equal(M.dot(statevector), simsystem.measurement.information_vector, decimal=3)
