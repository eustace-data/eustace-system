'''
Created on Jun 13, 2018

@author: hadcc
'''

import scipy.optimize

DEFAULT_METHOD = 'BFGS'
DEFAULT_TOL = 1e-7
DEFAULT_OPTIONS = {'disp': False, 'gtol': 1e-05, 'eps': 1e-08, 'return_all': False, 'maxiter': 20}

import time
def timer(function):
    """Simple decorator for printing run times of functions"""
    
    def timed_function(*args, **kwargs):
        
        tic = time.time()
        output = function(*args, **kwargs)
        toc = time.time()
        print function
        print toc - tic
        
        return output
    
    return timed_function

class Optimizer(object):
    """
    
    Setup and run optimization of a CostFunction
    
    """
    
    def __init__(self, initial_parameters, cost_function, method = DEFAULT_METHOD, tol = DEFAULT_TOL, options = DEFAULT_OPTIONS, cost_function_args = (), verbose = False):
        
        self.initial_parameters = initial_parameters
        self.cost_function = cost_function
        self.method   = method
        self.tol = tol  
        self.options  = options        
        self.cost_function_args = cost_function_args
        self.verbose  = verbose        
        self.callback = self.reporter
        
        self.optim_path = []
        self.result   = None
        
    def fit(self):
        """Run the cost function optimization"""
        
        self.optim_path = [self.initial_parameters,]
        
        fit_object = scipy.optimize.minimize( self.cost_function.cost_function,
                                              self.initial_parameters,
                                              self.cost_function_args,
                                              method=self.method,
                                              jac=self.cost_function.cost_derivative,
                                              tol=self.tol,
                                              callback = self.callback,
                                              options=self.options )
        
        self.result = fit_object
        
        if self.verbose:
            print self.result

    # setup callback function for recording path of optimisation 

    
    def reporter(self, current_state):
        """Reporter function to capture intermediate states of optimization."""
        
        if self.verbose:
            print "Current state:", current_state
            
        self.optim_path.append( current_state )
        

def demonstrate_optimisation():
    
    import numpy
    import scipy.sparse
    import eustace.advanced_standard.stats.spde.spherical
    import eustace.advanced_standard.linalg.costfunction
    
    # initialise SphereMeshSPDE
    my_spde = eustace.advanced_standard.stats.spde.spherical.SphereMeshSPDE(level = 6)
    dummy_spde = eustace.advanced_standard.stats.spde.spherical.SphereMeshSPDE(level = 3)    
    my_spde = eustace.advanced_standard.stats.spde.spherical.SphereMeshViewLocal(my_spde, 3, 10)
    
    
    
    import matplotlib.pyplot as plt
    dummy_spde.plot_triangles()
    my_spde.plot_triangles(ax = plt.gca())
    plt.show()
#    
    Qargs = (None, None, 2)
        
    
    # initialise prior mean and precision matrices
    Q_parameters = numpy.log( [1.0, 10.0 * numpy.pi / 180.0] )
    Q_prior = my_spde.build_Q( Q_parameters, *Qargs )
    
    mu_c = numpy.matrix(numpy.zeros((Q_prior.shape[0],1)))
    
    # generate random observations
    number_of_observations = 5000
    observation_sigma = 0.2
    
    locations = numpy.random.normal(0.0, 1.0, (number_of_observations,3))
    locations = ( locations.T * ( 1. / numpy.linalg.norm(locations, axis = 1) ) ).T
    
    design_matrix = my_spde.build_A(locations)
    
    Q_y = scipy.sparse.eye(number_of_observations) / observation_sigma**2
    observations = numpy.matrix( design_matrix * eustace.advanced_standard.linalg.sparse.sample_mvn(Q_prior, 1) + eustace.advanced_standard.linalg.sparse.sample_mvn(Q_y, 1) )

    # compute obervation increments
    i_increment = design_matrix.T * Q_y * observations
    Q_increment = design_matrix.T * Q_y * design_matrix
    
    # compute set_point as a first guess of the posterior using initial parameter values
    import sksparse.cholmod
    
    Q_prior = my_spde.build_Q( Q_parameters, *Qargs )
    
    Q_posterior = Q_prior + Q_increment
    Q_posterior_solver = sksparse.cholmod.cholesky(Q_posterior)
    mu_posterior = Q_posterior_solver.solve_A( Q_prior*mu_c + i_increment )
    
    #set_point = numpy.matrix( numpy.zeros( (number_of_variables,1) ) )
    set_point = numpy.asmatrix( mu_posterior )

    my_cost_function = eustace.advanced_standard.linalg.costfunction.MarginalGaussianCF( set_point, my_spde.build_Q, my_spde.build_dQdp ) 

    initial_hyperparameters = numpy.log( [2.0, 35.0 * numpy.pi / 180.0] )
    cf_args = (mu_c, i_increment, Q_increment, Qargs)
    
    # setup and run the optimizer
    my_optimizer = Optimizer(initial_hyperparameters, my_cost_function, method = 'BFGS', cost_function_args = cf_args, verbose = True)
    my_optimizer.fit()
    
    fit_object = my_optimizer.result
    optim_path = numpy.array( my_optimizer.optim_path )
    optim_path[:,0] = numpy.exp(optim_path[:,0])
    optim_path[:,1] = numpy.exp(optim_path[:,1]) * 180.0 / numpy.pi
    
    print optim_path
    
    print "number of local regions:", dummy_spde.n_triangles()
    print "number of variables per region:",  my_spde.n_latent_variables()
    
    plot_output = True
    if plot_output:
        # generate grid of parameter values over which to evaluate the cost function
        h0 = numpy.log( numpy.linspace(0.5, 3.0, 21) )
        h1 =  numpy.log( numpy.linspace(5, 40, 41) * numpy.pi / 180.0 )
        H0, H1 = numpy.meshgrid(h0, h1)
        
        # GP sparse
        ll = numpy.zeros(H0.shape).ravel()
    
        for n, qp in enumerate(zip(H0.ravel(),H1.ravel())):
            #print numpy.exp( qp )
            
            ll[n] = my_cost_function.cost_function(qp, *cf_args)
        
        
        import matplotlib.pyplot as plt
        
        plt.figure()
        # plot contour map of cost function
        ll = ll.reshape( H0.shape )
        plt.contour(numpy.exp( H0 ), numpy.exp(H1) * 180.0 / numpy.pi, ll, numpy.linspace(numpy.min(ll), numpy.max(ll), 100), cmap='jet')
        
        # plot parameter values at each step of the optimisation
        ps = numpy.array(optim_path)
        plt.plot(ps[:,0], ps[:,1], '-ro')
        plt.scatter(ps[-1,0], ps[-1,1], s = 80, marker='x')
        
        plt.show()
        
if __name__ == '__main__':
    demonstrate_optimisation()