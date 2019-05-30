import unittest
import numpy
import scipy.sparse
import eustace.analysis.advanced_standard.stats.spde.lattice
import eustace.analysis.advanced_standard.stats.spde.spherical
import eustace.analysis.advanced_standard.linalg.costfunction
import sksparse.cholmod

class TestMarginalGaussianCF(unittest.TestCase):
    
    def test_cost_function_derivative(self):
        """Check that the analytical cost function derivative is close to a numerical approximation of the derivative"""
        
        print "TestMarginalGaussianCF"
        
        # initialise either a SphereMeshSPDE or LatticeSPDE
        for mode in ['lattice']:
            
            if mode == 'lattice':
                dimension_specification = [[0., 10., 11],
                                           [0., 10., 11],]
                basis_function = eustace.analysis.advanced_standard.stats.spde.lattice.WendlandC4Basis
                overlap_factor = 2.5
                
                my_spde = eustace.analysis.advanced_standard.stats.spde.lattice.LatticeSPDE.construct(dimension_specification, basis_function, overlap_factor)
                Qargs = (None, None, 2, 1)
                
                # define wrapper functions to standardise inputs to the Q_function and its derivative with parameter_index as the last
                def Q_wrapper(Q_parameters, sigma_design_matrix, rho_design_matrix, alpha, H):
                    log_sigma = Q_parameters[0]
                    log_rho   = Q_parameters[1]
                    
                    return my_spde.build_Q_stationary(log_sigma, log_rho, alpha, H, sparse_format = 'csc')
                def dQdp_wrapper(Q_parameters, sigma_design_matrix, rho_design_matrix, alpha, H, parameter_index):
		    log_sigma = Q_parameters[0]
                    log_rho   = Q_parameters[1]
                    
                    return my_spde.build_dQdp_stationary(log_sigma, log_rho, alpha, H, parameter_index, sparse_format = 'csc')
                
                Q_function = Q_wrapper
                dQdp_function = dQdp_wrapper
                
            elif mode == 'sphere':
                my_spde = eustace.analysis.advanced_standard.stats.spde.spherical.SphereMeshSPDE(level = 3)
                Qargs = (None, None, 2)
            
            # initialise the precision matrix parameters
            parameters = numpy.ones(2) 
            
            # initialise the prior mean
            number_of_variables = my_spde.n_latent_variables()        
            prior_mean = numpy.matrix( numpy.zeros( (number_of_variables,1) ) )
            
            # generate observation increments
            i_increment = numpy.matrix( numpy.ones( (number_of_variables,1) ) )
            Q_increment = 2.0*scipy.sparse.eye(number_of_variables)
            
            # pick a "set_point" value for the latent variables at which the cost function is computed 
            set_point_type = 'initial_guess'
            if set_point_type == 'initial_guess':
                
                Q_prior = Q_function( parameters, *Qargs )
                
                Q_posterior = Q_prior + Q_increment
                Q_posterior_solver = sksparse.cholmod.cholesky(Q_posterior)
                mu_posterior = Q_posterior_solver.solve_A( Q_prior*prior_mean + i_increment )
            
                set_point = numpy.asmatrix( mu_posterior )
                
            elif set_point_type == 'zeros':
                set_point = numpy.matrix( numpy.zeros( (number_of_variables,1) ) )
            
            # initialise cost function object
            my_cost_function = eustace.analysis.advanced_standard.linalg.costfunction.MarginalGaussianCF( set_point, Q_function, dQdp_function ) 
    
            # compute numerical approximation of cost function gradient
            epsilon = 0.000001
            
            numerical_gradient = numpy.zeros( len(parameters) )
            for parameter_index in range( len(parameters) ):
                
                pshift = numpy.zeros(len(parameters))
                pshift[parameter_index] = epsilon
                
                numerical_gradient[parameter_index] = (my_cost_function.cost_function(parameters+pshift , prior_mean, i_increment, Q_increment, Qargs) \
                                                       - my_cost_function.cost_function(parameters-pshift , prior_mean, i_increment, Q_increment, Qargs) ) \
                                                        / (2.0*epsilon)
            # compute the analytical gradient
            analytical_gradient =  my_cost_function.cost_derivative(parameters, prior_mean, i_increment, Q_increment, Qargs)
            
            #print analytical_gradient
            #print numerical_gradient
            
            # check that the error is small
            numpy.testing.assert_allclose(analytical_gradient, numerical_gradient)
            
            
class TestMarginalGaussianReplicatesCF_InMemory(unittest.TestCase):
    
    def test_cost_function_derivative(self):
        """Check that the analytical cost function derivative is close to a numerical approximation of the derivative"""
        
        print "TestMarginalGaussianReplicatesCF_InMemory"
        
        # initialise either a SphereMeshSPDE or LatticeSPDE
        for mode in ['lattice']:
            
            if mode == 'lattice':
                dimension_specification = [[0., 10., 11],
                                           [0., 10., 11],]
                basis_function = eustace.analysis.advanced_standard.stats.spde.lattice.WendlandC4Basis
                overlap_factor = 2.5
                
                my_spde = eustace.analysis.advanced_standard.stats.spde.lattice.LatticeSPDE.construct(dimension_specification, basis_function, overlap_factor)
                Qargs = (None, None, 2, 1)
                
                # define wrapper functions to standardise inputs to the Q_function and its derivative with parameter_index as the last
                def Q_wrapper(Q_parameters, sigma_design_matrix, rho_design_matrix, alpha, H):
                    log_sigma = Q_parameters[0]
                    log_rho   = Q_parameters[1]
                    
                    return my_spde.build_Q_stationary(log_sigma, log_rho, alpha, H, sparse_format = 'csc')
                def dQdp_wrapper(Q_parameters, sigma_design_matrix, rho_design_matrix, alpha, H, parameter_index):
                    log_sigma = Q_parameters[0]
                    log_rho   = Q_parameters[1]
                    
                    return my_spde.build_dQdp_stationary(log_sigma, log_rho, alpha, H, parameter_index, sparse_format = 'csc')
                
                Q_function = Q_wrapper
                dQdp_function = dQdp_wrapper
                
            elif mode == 'sphere':
                my_spde = eustace.analysis.advanced_standard.stats.spde.spherical.SphereMeshSPDE(level = 3)
                Qargs = (None, None, 2)
            
            # initialise the precision matrix parameters
            parameters = numpy.ones(2) 
            
            # setup replicates
            
            n_replicates = 1
            
            set_points   = []
            prior_means  = []
            i_increments = []
            Q_increments = []
            
            for n in range(n_replicates):
                        
                # initialise the prior mean
                number_of_variables = my_spde.n_latent_variables()        
                prior_mean = numpy.matrix( numpy.zeros( (number_of_variables,1) ) )
                
                # generate observation increments
                i_increment = numpy.matrix( numpy.ones( (number_of_variables,1) ) )
                Q_increment = 2.0*scipy.sparse.eye(number_of_variables)
                
                # pick a "set_point" value for the latent variables at which the cost function is computed 
                set_point_type = 'initial_guess'
                if set_point_type == 'initial_guess':
                    
                    Q_prior = Q_function( parameters, *Qargs )
                    
                    Q_posterior = Q_prior + Q_increment
                    Q_posterior_solver = sksparse.cholmod.cholesky(Q_posterior)
                    mu_posterior = Q_posterior_solver.solve_A( Q_prior*prior_mean + i_increment )
                
                    set_point = numpy.asmatrix( mu_posterior )
                    
                elif set_point_type == 'zeros':
                    set_point = numpy.matrix( numpy.zeros( (number_of_variables,1) ) )
                    
                set_points.append(set_point)
                prior_means.append(prior_mean)
                i_increments.append(i_increment)
                Q_increments.append(Q_increment)
            
            # initialise cost function object
            my_cost_function = eustace.analysis.advanced_standard.linalg.costfunction.MarginalGaussianReplicateCF_InMemory( Q_function, dQdp_function ) 
    
            # compute numerical approximation of cost function gradient
            epsilon = 0.000001
            
            numerical_gradient = numpy.zeros( len(parameters) )
            for parameter_index in range( len(parameters) ):
                
                pshift = numpy.zeros(len(parameters))
                pshift[parameter_index] = epsilon
                
                numerical_gradient[parameter_index] = (my_cost_function.cost_function(parameters+pshift , set_points, prior_means, i_increments, Q_increments, Qargs) \
                                                       - my_cost_function.cost_function(parameters-pshift , set_points, prior_means, i_increments, Q_increments, Qargs) ) \
                                                        / (2.0*epsilon)
            # compute the analytical gradient
            analytical_gradient =  my_cost_function.cost_derivative(parameters, set_points, prior_means, i_increments, Q_increments, Qargs)
            
            #print analytical_gradient
            #print numerical_gradient
            
            # check that the error is small
            numpy.testing.assert_allclose(analytical_gradient, numerical_gradient)