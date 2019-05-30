"""Define linear (Gaussian) system and solutions."""

import scipy.sparse
import numpy
#from extendedcholmodwrapper import ExtendedCholmodWrapper as directsolver
from pardisowrapper import PardisoWrapper as directsolver

class MeasurementUpdate(object):

    @staticmethod
    def create_blank(dimension):
        """Make an empty measurement update object with the given state vector dimension."""

        return MeasurementUpdate(
            scipy.sparse.dok_matrix((dimension, dimension), dtype=numpy.float64),
            scipy.sparse.dok_matrix((dimension, 1), dtype=numpy.float64))

    @staticmethod
    def create_from_system_observations(design_jacobian, observation_precision, relative_observations):
        """Use Jacobian to compute observation increment."""

        # Calculate linear system elements
        design_obs_precision = design_jacobian.T.dot(observation_precision)
        information_precision = design_obs_precision.dot(design_jacobian)
        information_vector = design_obs_precision.dot(relative_observations.reshape(relative_observations.shape[0],1))

        # Build observation slice
        return MeasurementUpdate(information_precision, information_vector)

    def __init__(self, information_precision, information_vector):
        """Initialise with given precision increment and observation information."""

        self.information_precision = information_precision
        self.information_vector = information_vector

    def append_measurement_update(self, newupdate):
        """Extend with new evidence."""
        
        self.information_precision += newupdate.information_precision
        self.information_vector += newupdate.information_vector

class LinearSystem(object):
    """Linear (Gaussian) system comprises measurements with Gaussian error together with a prior covariance."""

    def __init__(self, prior):

        # Store prior
        self.prior = prior

        # Initialise with no evidence for the given size of system
        self.measurement = MeasurementUpdate.create_blank(prior.prior_number_of_state_parameters())

    def append_measurement_update(self, newupdate):

        self.measurement.append_measurement_update(newupdate)

class LinearSystemSolution(object):
    """Solution to linear (Gaussian) system and capabilities to generate samples from posterior."""

    def __init__(self, linearsystem):
        """Initialise and compute cholesky factorisation."""

        pass

    def maximum_a_posteriori(self):
        """Retrieve MAP estimate (which is also mean of distribution)."""

        raise NotImplementedError

    def white_noise_to_posterior_distribution(self, variate):
        """Apply the transform x = B.y where B is chosen such that xT.A.x = yT.y for all y."""

        raise NotImplementedError
    
    def random_sample(self, number_of_samples):
        """Generate samples from posterior distribution."""

        # Generate white noise
        variate = scipy.random.normal(0.0, 1.0, (self.number_of_state_parameters, number_of_samples))

        # Convert to sample from posterior
        return self.white_noise_to_posterior_distribution(variate)
    

class LinearSystemSolution_Cholesky(LinearSystemSolution):
    """Solution to linear (Gaussian) system and capabilities to generate samples from posterior."""

    def __init__(self, linearsystem, printstats=False):
        """Initialise and compute cholesky factorisation."""

        # Call base class
        super(LinearSystemSolution_Cholesky, self).__init__(linearsystem)

        # Evaluate cholesky decomposition
        self.posterior_precision = linearsystem.prior.prior_precision() + linearsystem.measurement.information_precision

        self.number_of_state_parameters = self.posterior_precision.shape[0]
        self.posterior_cholesky_factor = directsolver.cholesky(self.posterior_precision, printstats=printstats)
        self.observation_information = linearsystem.measurement.information_vector
        
    def maximum_a_posteriori(self):
        """Solve to get system mean."""

        return self.posterior_cholesky_factor.solve_A(self.observation_information)
        
    def marginal_std(self):
      """Compute marginal std from posterior precision for error propagation"""
      
      return numpy.sqrt(self.posterior_cholesky_factor.diag_inv_A())
      
    def approximated_marginal_std(self, number_of_samples=100):
      
      sample = self.random_sample(number_of_samples)

      return sample.std(axis=1)
      
    def white_noise_to_posterior_distribution(self, variate):
        """Apply the transform x = B.y where B is chosen such that xT.A.x = yT.y for all y."""

        return self.posterior_cholesky_factor.solve_backward_substitution(variate)
    
    def __del__(self):
        
        """If the object returned by the solve has clean up methods then call them"""
        free_memory_op = getattr(self, "free_memory", None)
        
        if callable(free_memory_op):
            free_memory_op(self.path.parent_op)
        
    
class LinearSystemSolution_CG(LinearSystemSolution):
    """Solution to linear (Gaussian) system and capabilities to generate samples from posterior."""

    counter = 0

    @staticmethod
    def callback(xk):

        # print 'LinearSystemSolution_CG iteration: ', LinearSystemSolution_CG.counter
        LinearSystemSolution_CG.counter += 1

    def __init__(self, linearsystem):
        """Initialise and compute cholesky factorisation."""

        # Call base class
        super(LinearSystemSolution_CG, self).__init__(linearsystem)

        # Evaluate posterior precision
        self.A = linearsystem.prior.prior_precision() + linearsystem.measurement.information_precision

        # And get information vector
        self.b = linearsystem.measurement.information_vector

        # Initial estimate just based on diagonal
        d = self.A.diagonal()
        
        # Precondition
        precondition = scipy.sparse.diags(1.0 / d, format='csc')

        # Initial estimate
        x0 = precondition.dot(self.b)

        # Solve to get maximum-a-posteriori
        x, info = scipy.sparse.linalg.cgs(self.A, self.b, x0, M=precondition, callback=LinearSystemSolution_CG.callback, tol=1.0E-6)

        # print 'LinearSystemSolution_CG iterations completed: ', self.counter

        # Solve to get MAP
        self.map = x

    def maximum_a_posteriori(self):
        """Solve to get system mean."""

        return self.map

    def white_noise_to_posterior_distribution(self, variate):
        """Apply the transform x = B.y where B is chosen such that xT.A.x = yT.y for all y."""

        # Relative to mean
        x = variate - self.map

        # Haven't done the maths on how to solve this yet(!)

        raise NotImplementedError
        
