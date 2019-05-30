"""Cost functions for optimisation of SPDE parameters"""

import numpy
import sksparse.cholmod
import sparse
import qinv

class CostFunction(object):
    """
    
    Abstract base class for cost functions
    
    """
    
    def cost_function(self, parameters, *args):
        """Evaluate cost function"""
        
        raise NotImplementedError
    
    def cost_derivative(self, parameters, *args):
        """Evaluate derivative of cost function with respect to its parameters"""
        
        raise NotImplementedError

class MarginalGaussianCF(CostFunction):
    """
    
    Marginal Log Likelihood cost function.
    
    This is computed using the conditional form where the log likelihood is 
    computed in it's conditional form at a fixed value of the state variable.
    
    """
    
    def __init__(self, set_point, Q_function, dQ_function):
        
        self.set_point = set_point
        self.Q_function = Q_function
        self.dQ_function = dQ_function
    
    def cost_function(self, parameters, prior_mean, i_increment, Q_increment, Qargs):
        """Compute cost function for given parameters"""
        
        # pre-compute prior and posterior precision matrices and setup linear system solvers
        prior_precision     = self.Q_function( parameters, *Qargs )
        posterior_precision = prior_precision + Q_increment

        prior_precision_solver     = sksparse.cholmod.cholesky( prior_precision )
        posterior_precision_solver = sksparse.cholmod.cholesky( posterior_precision )
        
        # compute posterior mean using information vector form of the calculation
        prior_i = prior_precision.dot( prior_mean )
        posterior_mean = posterior_precision_solver.solve_A( prior_i + i_increment )

        # compute terms of the marginal log likelihood function
        likelihood_term = 0.0   # Not actually zero but is independent of the parameters 
                                # and requires both the raw observations and design matrices 
                                # to compute. We instead neglect it as an additive constant. 
        
        prior_term  = self.gaussian_loglikelihood(parameters, prior_mean - self.set_point, prior_precision, prior_precision_solver)
        posterior_term  = self.gaussian_loglikelihood(parameters, posterior_mean - self.set_point, posterior_precision, posterior_precision_solver)

        cf_value = - likelihood_term - prior_term + posterior_term
        
        return numpy.asscalar(cf_value)

    def cost_derivative(self, parameters, prior_mean, i_increment, Q_increment, Qargs):
        """Derivative of cost function with respect to each of its parameters"""
        
        # initialise output
        number_of_parameters = len( parameters )
        cf_derivative = numpy.zeros( number_of_parameters )
        
        # pre-compute prior and posterior precision matrices and setup linear system solvers
        prior_precision     = self.Q_function( parameters, *Qargs )
        posterior_precision = prior_precision + Q_increment

        prior_precision_solver     = sksparse.cholmod.cholesky( prior_precision )
        posterior_precision_solver = sksparse.cholmod.cholesky( posterior_precision )
        
        # compute posterior mean using information vector form of the calculation
        prior_i = prior_precision.dot( prior_mean )
        posterior_mean = posterior_precision_solver.solve_A( prior_i + i_increment )
        
        # pre-compute partial inverses that will be used in cost function gradient calculations
        prior_Q_pinv = qinv.Q_inv( prior_precision, prior_precision_solver )
        posterior_Q_pinv = qinv.Q_inv( posterior_precision, posterior_precision_solver )
        
        # compute terms of the marginal log likelihood function derivatives
        likelihood_term = 0.0   # derivative of likelihood term is zero when using a fixed set point for conditioning
        
        # compute cost function derivative with respect to each parameter
        for parameter_index in range( number_of_parameters ):
        
            prior_precision_derivative = self.dQ_function( parameters, *Qargs + (parameter_index,) )
            posterior_precision_derivative = prior_precision_derivative
        
            posterior_mean_derivative = -posterior_precision_solver.solve_A( prior_precision_derivative * (posterior_mean - prior_mean) )
            
            prior_term     = self.gaussian_loglikelihood_derivative(parameters, 
                                                                    prior_mean - self.set_point, 
                                                                    prior_precision, 
                                                                    numpy.zeros(prior_mean.shape), 
                                                                    prior_precision_derivative, 
                                                                    prior_precision_solver, prior_Q_pinv)
                                                                    
            posterior_term = self.gaussian_loglikelihood_derivative(parameters, 
                                                                    posterior_mean - self.set_point, 
                                                                    posterior_precision, 
                                                                    posterior_mean_derivative, 
                                                                    posterior_precision_derivative, 
                                                                    posterior_precision_solver, 
                                                                    posterior_Q_pinv)
            
            cf_derivative[parameter_index] = - likelihood_term - prior_term + posterior_term
            
        return cf_derivative
    
    @staticmethod
    def gaussian_loglikelihood(parameters, mean_departure, precision, precision_solver):
        """Evaluate a gaussian log likelihood function"""
        
        n_variables = precision.shape[0]
        
        loglikelihood_value = - 0.5 * n_variables * numpy.log( 2.0 * numpy.pi ) + 0.5 * precision_solver.logdet() \
                              - 0.5 * mean_departure.T.dot( precision.dot( mean_departure ) )
        
        return loglikelihood_value
    
    
    @staticmethod
    def gaussian_loglikelihood_derivative(parameters, mean_departure, precision, mean_departure_derivative, precision_derivative, precision_solver, iQp = None):
        """Evaluate the derivative of a gaussian log likelihood function"""
        
        derivative_value = 0.5 * sparse.trace_iQ_dQdp(precision, precision_derivative, fQ = precision_solver, iQp = iQp) \
                            - 0.5 * mean_departure.T.dot( precision_derivative.dot( mean_departure ) ) \
                            - mean_departure_derivative.T.dot( precision.dot( mean_departure ) )
        
        return derivative_value



class MarginalGaussianReplicateCF_InMemory(MarginalGaussianCF):
    """
    
    Marginal Log Likelihood cost function for replicates of the same process.
    
    This is computed using the conditional form where the log likelihood is 
    computed in it's conditional form at a fixed value of the state variable.
    
    Here a set_point is specified for each replicate.
    
    """
    
    def __init__(self, Q_function, dQ_function):
        
        self.Q_function = Q_function
        self.dQ_function = dQ_function

    def cost_function(self, parameters, set_points, prior_means, i_increments, Q_increments, Qargs):
        """Compute cost function for given parameters"""
        
        print "parameters:", parameters
        
        # pre-compute prior and posterior precision matrices and setup linear system solvers
        prior_precision     = self.Q_function( parameters, *Qargs )
        prior_precision_solver     = sksparse.cholmod.cholesky( prior_precision )
        
        # initialise the cost function value
        cf_value = 0

        # loop over replicates evaluating the cost function
        n_replicates = len(set_points)
        
        for n in range(n_replicates):
        
            posterior_precision = prior_precision + Q_increments[n]
            posterior_precision_solver = sksparse.cholmod.cholesky( posterior_precision )
        
            # compute posterior mean using information vector form of the calculation
            prior_i = prior_precision.dot( prior_means[n] )
            posterior_mean = posterior_precision_solver.solve_A( prior_i + i_increments[n] )
            print posterior_mean
            # compute terms of the marginal log likelihood function
            likelihood_term = 0.0   # Not actually zero but is independent of the parameters 
                                    # and requires both the raw observations and design matrices 
                                    # to compute. We instead neglect it as an additive constant. 
        
            prior_term  = self.gaussian_loglikelihood(parameters, prior_means[n] - set_points[n], prior_precision, prior_precision_solver)
            posterior_term  = self.gaussian_loglikelihood(parameters, posterior_mean - set_points[n], posterior_precision, posterior_precision_solver)

            cf_value += - likelihood_term - prior_term + posterior_term
        
        
        print "cf_value:", cf_value
        
        return numpy.asscalar(cf_value)

    def cost_derivative(self, parameters, set_points, prior_means, i_increments, Q_increments, Qargs):
        """Derivative of cost function with respect to each of its parameters"""
        
        print "parameters:", parameters
        
        # initialise output
        number_of_parameters = len( parameters )
        
        # pre-compute prior precision matrix and setup linear system solver
        prior_precision     = self.Q_function( parameters, *Qargs )
        prior_precision_solver     = sksparse.cholmod.cholesky( prior_precision )
        
        # pre-compute prior partial inverses that will be used in cost function gradient calculations
        prior_Q_pinv = qinv.Q_inv( prior_precision, prior_precision_solver )
        
        # initialise the cost function value
        cf_derivative = numpy.zeros( number_of_parameters )

        # loop over replicates evaluating the cost function
        n_replicates = len(set_points)
        
        for n in range(n_replicates):
        
            # pre-compute posterior precision matrix for this replicate and setup linear system solver
            prior_precision     = self.Q_function( parameters, *Qargs )
            posterior_precision = prior_precision + Q_increments[n]
            posterior_precision_solver = sksparse.cholmod.cholesky( posterior_precision )
        
            # compute posterior mean using information vector form of the calculation
            prior_i = prior_precision.dot( prior_means[n] )
            posterior_mean = posterior_precision_solver.solve_A( prior_i + i_increments[n] )
        
            # pre-compute each prior precision derivative matrix
            prior_precision_derivatives = []
            for parameter_index in range( number_of_parameters ):
                prior_precision_derivatives.append( self.dQ_function( parameters, *Qargs + (parameter_index,) ) )
        
            # pre-compute posterior partial inverse that will be used in cost function gradient calculations
            posterior_Q_pinv = qinv.Q_inv( posterior_precision, posterior_precision_solver )
        
            # compute terms of the marginal log likelihood function derivatives
            likelihood_term = 0.0   # derivative of likelihood term is zero when using a fixed set point for conditioning
        
            # compute cost function derivative with respect to each parameter
            for parameter_index in range( number_of_parameters ):
            
                prior_precision_derivative = prior_precision_derivatives[parameter_index]
                posterior_precision_derivative = prior_precision_derivative
            
                posterior_mean_derivative = -posterior_precision_solver.solve_A( prior_precision_derivative * (posterior_mean - prior_means[n]) )
                
                prior_term     = self.gaussian_loglikelihood_derivative(parameters, 
                                                                        prior_means[n] - set_points[n], 
                                                                        prior_precision, 
                                                                        numpy.zeros(prior_means[n].shape), 
                                                                        prior_precision_derivative, 
                                                                        prior_precision_solver, prior_Q_pinv)
                                                                        
                posterior_term = self.gaussian_loglikelihood_derivative(parameters, 
                                                                        posterior_mean - set_points[n], 
                                                                        posterior_precision, 
                                                                        posterior_mean_derivative, 
                                                                        posterior_precision_derivative, 
                                                                        posterior_precision_solver, 
                                                                        posterior_Q_pinv)
                
                cf_derivative[parameter_index] += - likelihood_term - prior_term + posterior_term
        
        
        print "cf_derivative:", cf_derivative
        
        return cf_derivative

class MarginalGaussianReplicatesCF_Concept(object):
    """
    
    Marginal Log Likelihood cost function for replicates of the same process.
    
    This is computed using the conditional form where the log likelihood is 
    computed in it's conditional form at a fixed value of the state variable.
    
    Here a set_point is specified for each replicate.
    
    """

    def __init__(self, set_point, Q_function, dQ_function, Qargs):
        
        self.set_point = set_point
        self.Q_function = Q_function
        self.dQ_function = dQ_function
        
        self.cost_function_method = MarginalGaussianCF
        
        self.initalise_storage()
    
    def initialise_prior(self, parameters):
        
        self.precision = self.Q_function(parameters, Qargs)
    
    def Q_function():
        
        return self.prior_precision
    
    def dQ_function(parameter_index):
        
        return self.prior_precision_derivative[parameter_index]
    
    def initalise_storage(self):
        
        self.value = None
        self.value_derivative = None
        
        self.precision = None
        self.precision_derivative_list = []
        
        self.precision_solver = None
        self.precision_derivative_solver_list = []
        
        
    def cost_function_update(self, parameters, set_point, prior_mean, i_increment, Q_increment, Qargs):
        
        cost_function = MarginalGaussianCF(set_point, self.Q_function, self.dQ_function)
        
        self.value += self.cost_function.cost_function(parameters, prior_mean, i_increment, Q_increment, Qargs)
    
    def cost_derivative_update(self, set_point, parameters, prior_mean, i_increment, Q_increment, Qargs):
        
        cost_function = MarginalGaussianCF(set_point, self.Q_function, self.dQ_function)
        
        self.value_derivative += self.cost_function.cost_function(parameters, prior_mean, i_increment, Q_increment, Qargs)