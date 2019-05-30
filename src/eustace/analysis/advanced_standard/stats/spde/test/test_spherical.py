"""Tests for sphere mesh SPDE and original experimental display code."""

import unittest

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
import numpy
import scipy
from eustace.analysis.mesh.geometry import cartesian_to_polar2d
import matplotlib.pyplot

class TestSphereMeshSPDE(unittest.TestCase):

    @staticmethod
    def matern(distance, sigma, nu, rho):
        """Matern function which should be approximated by SPDE."""

        s = numpy.sqrt(2.0*nu) * distance / rho
        result = (sigma ** 2) * ((2 ** (1.0 - nu)) / scipy.special.gamma(nu)) * (s ** nu) * scipy.special.kv(nu, s)
        result[distance < 0.000001] = (sigma ** 2)
        return result

    def test_covariance_function(self):
        """Check that inverse precision follows matern covariance patterns."""
    
        # Generate
        spde = SphereMeshSPDE(level=5)

        # Index 12 in the mesh is 0 latitude, 0 longitude
        # - The covariance evaluated will relate to distance from this chosen point
        #   but the same shape of graph should be obtained for any point -
        #   there is nothing special about this choice.
        reference_index = 12

        # nu and rho as used in matern definition
        nu = 1.0
        rho = 0.3

        # Make precision
        precision = spde.build_Q_stationary(log_sigma=numpy.log(numpy.sqrt(1.0)), log_rho=numpy.log(2.0*rho), alpha=2)

        # Sample covariance at reference index

        # This is the vector with 1 in the reference location and zeros everywhere else
        # like [ 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 ... ]
        reference_unit = scipy.sparse.csc_matrix(([1.0], ([reference_index], [0])), shape=(precision.shape[0],1))

        # Solve Q.x = reference unit
        # This is like picking the column of inverse(Q) corresponding to reference index
        # and so is the covariance with respect to that reference location
        covariance_wrt_reference = scipy.sparse.linalg.spsolve(precision, reference_unit)

        # Great arc distances from reference point
        points = spde.triangulation.points
        refpoint = points[reference_index,:]
        dotprod = points.dot(refpoint)
        distance_from_reference = numpy.arccos(dotprod)

        # There is a scale factor difference between estimated covariance and Matern model
        sigma0 = numpy.sqrt(covariance_wrt_reference[reference_index])

        # Compute disparity (taking into account scale factor difference)
        disparity = max(numpy.abs(covariance_wrt_reference - TestSphereMeshSPDE.matern(distance_from_reference, sigma0, nu, rho)))

        # Convert to percentage and check
        disparity_percent = 100.0 * disparity / (sigma0 ** 2)
        # print 'Disparity (%): ', disparity_percent
        self.assertTrue(disparity_percent < 2.0)

    def test_derivative(self):
	"""Check that numerical and analytical derivatives of Q match."""
	
        spde = SphereMeshSPDE(level=3)
        Q = spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0), 2)

        epsilon = 0.000001

        expected_dQ0 = \
            ((spde.build_Q_stationary(numpy.log(1.0)+epsilon, numpy.log(1.0), 2) - 
              spde.build_Q_stationary(numpy.log(1.0)-epsilon, numpy.log(1.0), 2)) / (2.0 * epsilon)).todense()

        expected_dQ1 = \
            ((spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0)+epsilon, 2) - 
              spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0)-epsilon, 2)) / (2.0 * epsilon)).todense()

        dQ0 = spde.build_dQdp_stationary(numpy.log(1.0), numpy.log(1.0), 2, 0).todense()
        dQ1 = spde.build_dQdp_stationary(numpy.log(1.0), numpy.log(1.0), 2, 1).todense()

        # print numpy.abs(dQ0 - expected_dQ0).ravel().max() / numpy.abs(expected_dQ0).ravel().max()
        # print numpy.abs(dQ1 - expected_dQ1).ravel().max() / numpy.abs(expected_dQ1).ravel().max()
        
        numpy.testing.assert_almost_equal(dQ0, expected_dQ0, decimal=7)
        numpy.testing.assert_almost_equal(dQ1, expected_dQ1, decimal=7)

#------- Original experimental display code is below

def spde_demo():
    """
    Demonstrate construction of an SPDE model on a tringular mesh on a sphere.
    
    The mesh is constructed as a subdivided icosahedron. Observations are simulated
    and a latent variable model is learnt.
    
    """
    
    import matplotlib.pyplot as plt
    from advanced_standard.linalg import sparse
    from advanced_standard.utilities import plot_tools
    
    print "############"
    print "# Setup    #"
    print "############"
    
    my_spde =  SphereMeshSPDE(level = 4, project_onto_sphere = True)
    print "number of latent variables:"
    print my_spde.n_latent_variables()
    
    print "############"
    print "# Q matrix #"
    print "############"
    
    Q = my_spde.build_Q(Q_parameters = np.log(np.array([2., 0.5])), alpha = 2)
    #dQdp0 = my_spde.build_dQdp(Q_parameters = np.log(np.array([1., 0.5])), alpha = 2, parameter_index = 0)
    #dQdp1 = my_spde.build_dQdp(Q_parameters = np.log(np.array([1., 0.5])), alpha = 2, parameter_index = 1)
    print "############"
    print "# Q sample #"
    print "############"

    c = sparse.sample_mvn(Q, 1)
    
    print "############"
    print "# Obs gen  #"
    print "############"
    
    number_of_observations = 2592
    test_points = np.random.uniform(-1,1,(number_of_observations, 3))
    
    standard_error = 0.5
    Q_e = scipy.sparse.eye(number_of_observations) * (1.0 / standard_error**2)
    
    observation_error = sparse.sample_mvn(Q_e, 1)
    
    print "############"
    print "# A matrix #"
    print "############"
    
    A = my_spde.build_A(test_points)
    truth = A*c
    
    y = truth + observation_error
    
    print "############"
    print "# solving #"
    print "############"
    
    Q_estimate = sparse.compute_Q_u(Q_e, A, Q)
    i_estimate = sparse.compute_i_u(y, Q_e, A, np.zeros((my_spde.n_latent_variables())))
    c_estimate = sparse.compute_mu_u(i_estimate, Q_estimate)
    
    reconstruction = A * c_estimate
      
    print "############"
    print "# plotting #"
    print "############"
    
    vmin = -3
    vmax = 3
    
    plt.figure()
    vertices_coords_2d =cartesian_to_polar(my_spde.triangulation.vertices)
    plot_tools.scatter2d(vertices_coords_2d[:,1]* 180/np.pi,
                             vertices_coords_2d[:,0]* 180/np.pi, cs=c, vmax=vmax, vmin=vmin)
    
    plt.figure()
    test_point_coords_2d = cartesian_to_polar(test_points)
    plot_tools.scatter2d(test_point_coords_2d[:,1]* 180/np.pi,
                             test_point_coords_2d[:,0]* 180/np.pi, cs=truth, vmax=vmax, vmin=vmin)
    
    plt.figure()
    test_point_coords_2d = cartesian_to_polar(test_points)
    plot_tools.scatter2d(test_point_coords_2d[:,1]* 180/np.pi,
                             test_point_coords_2d[:,0]* 180/np.pi, cs=y, vmax=vmax, vmin=vmin)
    
    plt.figure()
    vertices_coords_2d =cartesian_to_polar(my_spde.triangulation.vertices)
    plot_tools.scatter2d(vertices_coords_2d[:,1]* 180/np.pi,
                             vertices_coords_2d[:,0]* 180/np.pi, cs=c_estimate.ravel(), vmax=vmax, vmin=vmin)
    
    plt.figure()
    test_point_coords_2d = cartesian_to_polar(test_points)
    plot_tools.scatter2d(test_point_coords_2d[:,1]* 180/np.pi,
                             test_point_coords_2d[:,0]* 180/np.pi, cs=reconstruction.ravel(), vmax=vmax, vmin=vmin)
    
    plt.show()
    
def optimisation_demo():
    """
    Demonstrate construction of an SPDE model on a tringular mesh on a sphere.
    
    The mesh is constructed as a subdivided icosahedron. Observations are simulated
    and a latent variable model is learnt.
    
    """
    
    import matplotlib.pyplot as plt
    from advanced_standard.linalg import sparse
    from advanced_standard.utilities import plot_tools
    
    print "############"
    print "# Setup    #"
    print "############"
    
    my_spde =  SphereMeshSPDE(level = 4, project_onto_sphere = True)
    print "number of latent variables:"
    print my_spde.n_latent_variables()
    
    print "############"
    print "# Q matrix #"
    print "############"
    
    Q_parameters = np.log(np.array([2., 1.0]))
    print "Q_parameters:", Q_parameters
    Q = my_spde.build_Q(Q_parameters = Q_parameters, alpha = 2)
    #dQdp0 = my_spde.build_dQdp(Q_parameters = np.log(np.array([1., 0.5])), alpha = 2, parameter_index = 0)
    #dQdp1 = my_spde.build_dQdp(Q_parameters = np.log(np.array([1., 0.5])), alpha = 2, parameter_index = 1)
    print "############"
    print "# Q sample #"
    print "############"

    c = sparse.sample_mvn(Q, 1)
    
    print "############"
    print "# Obs gen  #"
    print "############"
    
    number_of_observations = 25000#2592
    test_points = np.random.uniform(-1,1,(number_of_observations, 3))
    
    standard_error = 1.0
    Q_e = scipy.sparse.eye(number_of_observations) * (1.0 / standard_error**2)
    
    observation_error = sparse.sample_mvn(Q_e, 1)
    
    print "############"
    print "# A matrix #"
    print "############"
    
    A = my_spde.build_A(test_points)
    truth = A*c
    
    #y = np.atleast_2d(truth + observation_error).T
    y = truth + observation_error
    print "############"
    print "# solving #"
    print "############"
    
    Q_increment = A.T * Q_e * A
    i_increment = A.T * Q_e * y
    
    import advanced_standard.linalg.optimisation as optimisation
    
    Qkwargs = {'alpha': 2, 'sparse_format': 'csc'}
    
#    interval = 1e-5
#    
#    cost_function = optimisation.cost_function(Q_parameters=Q_parameters,
#                                               i_c = np.zeros( (Q.shape[0],1) ),
#                                               i_increment= i_increment,
#                                               Q_increment=Q_increment,
#                                               Q_function = my_spde.build_Q,
#                                               dQ_function = my_spde.build_dQdp, set_point_mode = 'fixed',
#                                               fixed_set_point = np.zeros( (Q.shape[0],1) ),
#                                               observation_loglik=None,
#                                               **Qkwargs)
#    
#    vost_function = optimisation.cost_function(Q_parameters=Q_parameters +np.array([0.0, interval]),
#                                               i_c = np.zeros( (Q.shape[0],1) ),
#                                               i_increment= i_increment,
#                                               Q_increment=Q_increment,
#                                               Q_function = my_spde.build_Q,
#                                               dQ_function = my_spde.build_dQdp, set_point_mode = 'fixed',
#                                               fixed_set_point = np.zeros( (Q.shape[0],1) ),
#                                               observation_loglik=None,
#                                               **Qkwargs)
#    
#    dcost_function = optimisation.cost_function_derivative(Q_parameters=Q_parameters,
#                                               i_c = np.zeros( (Q.shape[0],1) ),
#                                               i_increment= i_increment,
#                                               Q_increment=Q_increment,
#                                               Q_function = my_spde.build_Q,
#                                               dQ_function = my_spde.build_dQdp, set_point_mode = 'fixed',
#                                               fixed_set_point = np.zeros( (Q.shape[0],1) ),
#                                               **Qkwargs)
#    
#    
#    print (vost_function - cost_function) / interval
#    
#    print dcost_function
#    print cost_function
    
    Q_init = np.log( np.array([2.0, np.pi/4]) )
    
    mu_c = np.zeros( (Q.shape[0],1) )
    i_increment= i_increment
    Q_increment=Q_increment
    Q_function = my_spde.build_Q
    dQ_function = my_spde.build_dQdp
    set_point_mode = 'fixed'
    fixed_set_point = np.zeros( (Q.shape[0],1) )
    observation_loglik = None
    
    opt_args = (mu_c, i_increment, Q_increment, Q_function, dQ_function, set_point_mode, fixed_set_point, observation_loglik, Qkwargs)
    
    
    print scipy.optimize.check_grad(optimisation.cost_function, optimisation.cost_function_derivative, Q_init, *opt_args)
    
    def opt_printer(x):
        print x
    
    fit = scipy.optimize.minimize(optimisation.cost_function,
                                                    Q_init,
                                                    args=opt_args,
                                                    method='BFGS',
                                                    jac=optimisation.cost_function_derivative,
                                                    tol=None,
                                                    callback = opt_printer,
                                                    options={'disp': False, 'gtol': 1e-04, 'eps': 1e-08, 'return_all': False, 'maxiter': 20})
    
    print fit
    
    fitted_hyperparamters = fit['x']
    
    new_Q = Q_function(Q_parameters = fitted_hyperparamters, **Qkwargs)
    
    Q_estimate = sparse.compute_Q_u(Q_e, A, new_Q)
    i_estimate = sparse.compute_i_u(y, Q_e, A, np.zeros((my_spde.n_latent_variables(),1)))
    c_estimate = sparse.compute_mu_u(i_estimate, Q_estimate)
    
    reconstruction = A * c_estimate
      
    print "############"
    print "# plotting #"
    print "############"
    
    vmin = -6
    vmax = 6
    
    plt.figure()
    vertices_coords_2d =cartesian_to_polar(my_spde.triangulation.vertices)
    plot_tools.scatter2d(vertices_coords_2d[:,1]* 180/np.pi,
                             vertices_coords_2d[:,0]* 180/np.pi, cs=c.ravel(), vmax=vmax, vmin=vmin)
    
    plt.figure()
    test_point_coords_2d = cartesian_to_polar(test_points)
    plot_tools.scatter2d(test_point_coords_2d[:,1]* 180/np.pi,
                             test_point_coords_2d[:,0]* 180/np.pi, cs=truth.ravel(), vmax=vmax, vmin=vmin)
    
    plt.figure()
    test_point_coords_2d = cartesian_to_polar(test_points)
    plot_tools.scatter2d(test_point_coords_2d[:,1]* 180/np.pi,
                             test_point_coords_2d[:,0]* 180/np.pi, cs=np.asarray(y).ravel(), vmax=vmax, vmin=vmin)
    
    plt.figure()
    vertices_coords_2d =cartesian_to_polar(my_spde.triangulation.vertices)
    plot_tools.scatter2d(vertices_coords_2d[:,1]* 180/np.pi,
                             vertices_coords_2d[:,0]* 180/np.pi, cs=c_estimate.ravel(), vmax=vmax, vmin=vmin)
    
    plt.figure()
    test_point_coords_2d = cartesian_to_polar(test_points)
    plot_tools.scatter2d(test_point_coords_2d[:,1]* 180/np.pi,
                             test_point_coords_2d[:,0]* 180/np.pi, cs=reconstruction.ravel(), vmax=vmax, vmin=vmin)
    
    plt.show()
    
def mesh_neighbour_demo():
    
    import advanced_standard.stats.spde.convex
    import matplotlib.pyplot as plt
    
    spde_level = 3
    neighbour_level = 1
    region_centre_index = 0
    
    my_spde =  SphereMeshSPDE(level = spde_level, project_onto_sphere = True)
    in_indices = my_spde.neighbours_at_level(neighbour_level, region_centre_index)
    vertex_indices = np.unique( my_spde.triangulation.triangles[in_indices,1:].ravel() )
    
    convex_mesh_triangulation = advanced_standard.stats.spde.convex.ConvexMeshSPDE(my_spde.triangulation.vertices[vertex_indices,:])
    convex_mesh_triangulation.plot_triangles(linewidths=0.5)
    
    #my_spde.plot_triangles()
    
    my_spde =  SphereMeshSPDE(level = neighbour_level, project_onto_sphere = True)
    in_indices = my_spde.neighbours_at_level(neighbour_level, region_centre_index)
    vertex_indices = np.unique( my_spde.triangulation.triangles[in_indices,1:].ravel() )
    
    ax = plt.gca()
    convex_mesh_triangulation = advanced_standard.stats.spde.convex.ConvexMeshSPDE(my_spde.triangulation.vertices[vertex_indices,:])
    convex_mesh_triangulation.plot_triangles(colors = 'r', ax=ax, linewidths=1.0)
    
    
    centre_vertex_indices = np.unique( my_spde.triangulation.triangles[region_centre_index,1:].ravel() )
    centre_vertex_triangulation = advanced_standard.stats.spde.convex.ConvexMeshSPDE(my_spde.triangulation.vertices[centre_vertex_indices,:])
    centre_vertex_triangulation.plot_triangles(colors = 'b', ax=ax, linewidths=1.5)
    plt.show()
    #print low_resolution_mesh_indices
    
if __name__ == '__main__':
    
    #spde_demo()
    #optimisation_demo()  
    mesh_neighbour_demo()  
