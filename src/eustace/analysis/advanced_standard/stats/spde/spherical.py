'''
Created on May 10, 2017

@author: hadcc

Updated 20170812 JRM

'''

import numpy as np
import scipy.sparse
from eustace.analysis.mesh.mesh import MeshIcosahedron
from eustace.analysis.mesh.mesh import MeshIcosahedronSubdivision
from eustace.analysis.mesh.locator import TriangleLocatorHierarchical
from eustace.analysis.mesh.geometry import polar2d_to_cartesian

class SphereMeshSPDE(object):
    """
    
    Compute precision matrices resulting from SPDE on spherical mesh.
    
    """
    
    def __init__(self, level, sparse_format='csr'):
        """Build for specified level of subdivision of icosahedron."""
        
        self.triangulation = MeshIcosahedronSubdivision.build(level)
        self.triangle_locator = TriangleLocatorHierarchical(MeshIcosahedronSubdivision.build(0))
        self.sparse_format = sparse_format
    
    @staticmethod
    def n_triangles_at_level(level):
        """Utility routine for obtaining the number of triangulations at a given level.
        
        Intended to be used to infer the number of views that make up a full sphere.
        
        """
        
        n_triangles = 20
        
        for this_level in range(1,level+1):        
            n_triangles = n_triangles * 4
            
        return n_triangles
        
    def n_latent_variables(self):
        """Return number of latent variables in the SPDE model.
        
        Equal to the number of basis functions/number of nodes in the
        triangular mesh.
        
        """
        
        return self.triangulation.points.shape[0]

    def build_A(self, polar_coordinates):
        """
        Takes polar coordinates at which to evaluate basis functions (at triangulation level if specified) and returns the A matrix of basis weights.
        Polar coordinates should be in NumPy 2D array with 2 rows (for latitude, longitude) and one column per point.
        
        """
        
        # Convert to cartesian coordinates on unit sphere
        points = polar2d_to_cartesian(polar_coordinates)
        
        # Evaluate weights by computing barycentres
        triangle_id, basis_weights = self.triangle_locator.locate_vector_hierarchical(points, levels=int(self.triangulation.level))

        # Triangle numbers are for maximum subdivisions - revise to use only the subdivisions required
        triangle_index = triangle_id >> (MeshIcosahedron.BASE_HASH_SHIFT - np.uint64(2*self.triangulation.level))
        
        # Put into matrix
        n_bases = self.n_latent_variables()
        n_locations = polar_coordinates.shape[0]
        point_indices  = np.repeat(np.arange(n_locations), 3).ravel() 
        vertex_indices = self.triangulation.triangles[triangle_index,1:].ravel()
        A = scipy.sparse.dok_matrix((n_locations, n_bases))
        A[point_indices, vertex_indices] = basis_weights.ravel()
        return A.asformat(self.sparse_format)


    @staticmethod
    def triangle_area(vertices, triangle_vertex_indices):
        """
        Returns the area of a triangle with corners with coordinates in rows of triangle_vertex_indices.
        
        Args:
        
            * vertices:
                Locations of triangle vertices for all triangles with dimensions [n_vertices, 3]
                where each row contains the 3D coordinates of a vertex.
                
            * triangle_vertex_indices:
                Indices to the vertices of a set of triangles with dimensions [n_triangles, 3]
                where each row contains the indices of the three corners of a triangle.
        
        """
        
        all_vertices = vertices[triangle_vertex_indices.ravel(),:] # ordered as [triangle vertex number, 3d vertex coordinates] with each consecutive 3 rows corresponding to a single triangle
    
        edges_01  = all_vertices[1::3,:] - all_vertices[0::3,:]
        edges_02 = all_vertices[2::3,:] - all_vertices[0::3,:]
        
        triangle_areas = 0.5 * np.linalg.norm(np.cross(edges_01, edges_02), axis = 1) 
        
        return triangle_areas
    
    @staticmethod
    def triangle_edges(vertices, triangle_vertex_indices):
        """
        Returns the edge vectors of triangles with corners with coordinates in rows of triangle_vertex_indices.
        
        Returned edges are ordered such that the zeroth edge is opposite the zeroth vertex, first edge is 
        opposite the first vertex and second edge is opposite the second vertex. 
        
        Args:
        
            * vertices:
                Locations of triangle vertices for all triangles with dimensions [n_vertices, 3]
                where each row contains the 3D coordinates of a vertex.
                
            * triangle_vertex_indices:
                Indices to the vertices of a set of triangles with dimensions [n_triangles, 3]
                where each row contains the indices of the three corners of a triangle.
        
        Returns:
            
            * edges:
                Returns a 3D array of triangle edge vectors formatted as [3*n_triangles,3] where
                [0::3,:] contains the edge vectors for the zeroth edge of each triangle,
                [1::3,:] contains the edge vectors for the first edge of each triangle and
                [2::3,:] contains the edge vectors for the last edge of each triangle.
        
        """
        
        n_triangles = triangle_vertex_indices.shape[0]
        all_vertices = vertices[triangle_vertex_indices.ravel(),:] # ordered as [triangle vertex number, 3d vertex coordinates] with each consecutive 3 rows corresponding to a single triangle
        
        edges = np.zeros((3*n_triangles,3))
        
        edges[0::3] = all_vertices[2::3] - all_vertices[1::3]
        edges[1::3] = all_vertices[0::3] - all_vertices[2::3]
        edges[2::3] = all_vertices[1::3] - all_vertices[0::3]
        
        return edges
        
    def sigma_design_matrix_stationary(self):
        """Get sigma matrix for stationary process."""

        number_of_latent_variables = self.n_latent_variables()        
        
        sigma_design_matrix = scipy.sparse.dok_matrix((number_of_latent_variables, 2))
        sigma_design_matrix[:,0] = 1.0
        
        return sigma_design_matrix.asformat( self.sparse_format )

    def rho_design_matrix_stationary(self):
        """Get rho matrix for stationary process."""
        
        number_of_latent_variables = self.n_latent_variables()        
        
        rho_design_matrix = scipy.sparse.dok_matrix((number_of_latent_variables, 2))
        rho_design_matrix[:,1] = 1.0
        
        return rho_design_matrix.asformat( self.sparse_format )
    
    def build_Q_stationary(self, log_sigma, log_rho, alpha):
        """Call build_Q for the stationary case."""
        
        # Vector of parameters
        Q_parameters = np.array([ log_sigma, log_rho ])
        
        return self.build_Q(Q_parameters, alpha, self.sigma_design_matrix_stationary(), self.rho_design_matrix_stationary())

    def build_dQdp_stationary(self, log_sigma, log_rho, alpha, parameter_index):
        """Call build_dQdp for the stationary case."""        
        
        # Vector of parameters
        Q_parameters = np.array([ log_sigma, log_rho ])
        
        return self.build_dQdp(Q_parameters, alpha, self.sigma_design_matrix_stationary(), self.rho_design_matrix_stationary(), parameter_index)

    def sigma_design_matrix_expanded(self):
        """Get sigma matrix for non-stationary process with expanded parameterisation."""

        number_of_latent_variables = self.n_latent_variables()        
        
        sigma_design_matrix = scipy.sparse.hstack( [scipy.sparse.eye(number_of_latent_variables),
                                                    scipy.sparse.dok_matrix((number_of_latent_variables, number_of_latent_variables))] )
        
        return sigma_design_matrix.asformat( self.sparse_format )

    def rho_design_matrix_expanded(self):
        """Get rho matrix for non-stationary process with expanded parameterisation."""
        
        number_of_latent_variables = self.n_latent_variables()        
        
        rho_design_matrix = scipy.sparse.hstack( [scipy.sparse.dok_matrix((number_of_latent_variables, number_of_latent_variables)),
                                                  scipy.sparse.eye(number_of_latent_variables)] )
        
        return rho_design_matrix.asformat( self.sparse_format )

    def build_Q_expanded(self, log_sigma, log_rho, alpha):
        """Call build_Q for the non-stationary case for parameter expanded representation with parameters specified at each node."""
        
        # Vector of parameters
        Q_parameters = np.concatenate([ log_sigma, log_rho ])
        
        return self.build_Q(Q_parameters, alpha, self.sigma_design_matrix_expanded(), self.rho_design_matrix_expanded())
    
    def build_dQdp_expanded(self, log_sigma, log_rho, alpha, parameter_index):
        """Call build_dQdp for the non-stationary case for parameter expanded representation with parameters specified at each node.
        
        The expanded representation is intended for use for non-stationary models with pre-computed parameters.
        
        It is strongly recommended that the expanded representation is not used for optimisation except for
        trivially small models as it would optimise a (sigma, rho) pair at each individual node!!!
        
        """       
        
        # Vector of parameters
        Q_parameters = np.concatenate([ log_sigma, log_rho ])
        
        return self.build_dQdp(Q_parameters, alpha, self.sigma_design_matrix_expanded(), self.rho_design_matrix_expanded(), parameter_index)
    
    
    def build_Q(self, Q_parameters, alpha, sigma_design_matrix, rho_design_matrix):
        """Build precision matrix for regression type non-stationary precision model.
        
        log( \sigma ) = \sigma_0 + \sum_{k}^p b_k theta_k
        
        log( \rho ) = \rho_0 + \sum_{k}^p c_k theta_k
        
        In terms of input kwargs this may be computed as:
        
        log( \sigma ) = \sigma_0 + sigma_design_matrix * Q_parameters
        log( \rho )   = \rho_0   + rho_design_matrix * Q_parameters
        
        b_k values are stored in the kth column of sigma_design_matrix,
        c_k values are stored in the kth column of rho_design_matrix, and 
        theta_k is the value in the kth row of Q_parameters. \sigma_0, \rho_0 are 
        fixed internally to each equal zero.
        
        By default a stationary precision is constructed such that:
        
        log( \sigma ) = \sigma_0 + b_0 theta_0        
        log( \rho )   = \rho_0   + b_1 theta_1
        
        For default setup, with sigma_design_matrix=None and rho_design_matrix = None, Q_parameters 
        are equivalent to log(sigma), log(rho). This setup permits optimisation of log(sigma), log(rho)
        and has a reasonably intuitive interpretation of Q_parameters. In the more general case, the
        interpretation of the parameters is less intuitive.
        
        This default setup of Q_parameters = np.array([[0.0, 0.0]]).T results in sigma = 1, rho = 1,
        i.e. theta_1 = exp(1), sigma = exp(1).
        
        Note: the default stationary model has been removed from the implementation. Convenience
        methods are instead provided to construct the design matrices for a stationary model. Methods
        are also provided for an expanded parameter representation when sigma and rho are specified
        at each individual node.
        
        Internally, the (sigma, rho) parameterisation is converted into the native (kappa, tau) parameterisation,
        by computation of \kappa_0, kappa_design_matrix, \tau_0, and tau_design_matrix. This forms a similar
        pair of log linear models that use the same Q_parameters matrix as the (sigma, rho) parameterisation.
        
        log( \kappa ) = \kappa_0 + kappa_design_matrix * Q_parameters
        log( \tau )   = \tau_0   + tau_design_matrix * Q_parameters
        
        """
        
        from scipy.special import gamma

        # convert sigma, rho parameterisation into kappa, tau
        d = 2   # domain dimension for manifold
        nu = alpha - d / 2.0
        
        if 2 * alpha <= d:
            raise ValueError('alpha must be greater than the manifold dimension / 2')
        
        # Convert provided design matrices into appropriate design matrices for the internal (tau, kappa) representation.        
        b0_sigma = 0.0 # Equivalent to a unity base scale factor on sigma
        b0_rho   = 0.0 # Equivalent to a unity base scale factor on rho
        
        b0_kappa = np.log( 8.0 * nu ) / 2.0 - b0_rho
        b0_tau   = 0.5 * np.log( gamma(nu) / (gamma(alpha)* (4*np.pi)**(0.5*d) ) ) - b0_sigma - nu * b0_kappa
        
        kappa_design_matrix = -rho_design_matrix
        tau_design_matrix   = -sigma_design_matrix - nu * kappa_design_matrix
        
        # compute non-stationary kappa and tau vectors
        kappa = np.exp( b0_kappa + kappa_design_matrix.dot(Q_parameters) ).ravel()
        tau   = np.exp( b0_tau + tau_design_matrix.dot(Q_parameters) ).ravel()
    
        kappa = scipy.sparse.diags(kappa, format=self.sparse_format)
        tau = scipy.sparse.diags(tau, format=self.sparse_format)
    
        Q = self.build_Q_native(kappa=kappa, alpha= alpha, tau = tau)        

        return Q
    
    def build_dQdp(self, Q_parameters, alpha, sigma_design_matrix, rho_design_matrix, parameter_index):
        """Build precision matrix for regression type non-stationary precision model.
        
        log( \sigma ) = \sigma_0 + \sum_{k}^p b_k theta_k
        
        log( \rho ) = \rho_0 + \sum_{k}^p c_k theta_k
        
        In terms of input kwargs this may be computed as:
        
        log( \sigma ) = \sigma_0 + sigma_design_matrix * Q_parameters
        log( \rho )   = \rho_0   + rho_design_matrix * Q_parameters
        
        b_k values are stored in the kth column of sigma_design_matrix,
        c_k values are stored in the kth column of rho_design_matrix, and 
        theta_k is the value in the kth row of Q_parameters. \sigma_0, \rho_0 are 
        fixed internally to each equal zero.
        
        By default a stationary precision is constructed such that:
        
        log( \sigma ) = \sigma_0 + b_0 theta_0        
        log( \rho )   = \rho_0   + b_1 theta_1
        
        For default setup, with sigma_design_matrix=None and rho_design_matrix = None, Q_parameters 
        are equivalent to log(sigma), log(rho). This setup permits optimisation of log(sigma), log(rho)
        and has a reasonably intuitive interpretation of Q_parameters. In the more general case, the
        interpretation of the parameters is less intuitive.
        
        This default setup of Q_parameters = np.array([[0.0, 0.0]]).T results in sigma = 1, rho = 1,
        i.e. theta_1 = exp(1), sigma = exp(1).
        
        Internally, the (sigma, rho) parameterisation is converted into the native (kappa, tau) parameterisation,
        by computation of \kappa_0, kappa_design_matrix, \tau_0, and tau_design_matrix. This forms a similar
        pair of log linear models that use the same Q_parameters matrix as the (sigma, rho) parameterisation.
        
        log( \kappa ) = \kappa_0 + kappa_design_matrix * Q_parameters
        log( \tau )   = \tau_0   + tau_design_matrix * Q_parameters
        
        """
        
        from scipy.special import gamma

        # convert sigma, rho parameterisation into kappa, tau
        d = 2   # domain dimension for manifold
        nu = alpha - d / 2.0
        
        if 2 * alpha <= d:
            raise ValueError('alpha must be greater than the manifold dimension / 2')
        
        # Convert provided design matrices into appropriate design matrices for the internal (tau, kappa) representation.        
        b0_sigma = 0.0 # Equivalent to a unity base scale factor on sigma
        b0_rho   = 0.0 # Equivalent to a unity base scale factor on rho
        
        b0_kappa = np.log( 8.0 * nu ) / 2.0 - b0_rho
        b0_tau   = 0.5 * np.log( gamma(nu) / (gamma(alpha)* (4*np.pi)**(0.5*d) ) ) - b0_sigma - nu * b0_kappa
        
        kappa_design_matrix = -rho_design_matrix
        tau_design_matrix   = -sigma_design_matrix - nu * kappa_design_matrix
        
        # compute non-stationary kappa and tau vectors evaluated at each vertex
        kappa = np.exp( b0_kappa + kappa_design_matrix.dot( Q_parameters) ).ravel()
        tau   = np.exp( b0_tau + tau_design_matrix.dot( Q_parameters) ).ravel()
    
        # compute derivatives of kappa and tau w.r.t. the specified hyperparameter
        dQ_parameters = np.zeros(Q_parameters.shape)
        dQ_parameters[parameter_index] = 1.0    
        dkappa = (kappa_design_matrix.dot(dQ_parameters) * kappa).ravel()
        dtau   = (tau_design_matrix.dot(dQ_parameters) * tau).ravel()
    
        # remap kappa and tau and their derivatives into diagonal sparse matrices
        kappa = scipy.sparse.diags(kappa, format=self.sparse_format)
        tau = scipy.sparse.diags(tau, format=self.sparse_format)
        
        dkappa = scipy.sparse.diags(dkappa, format=self.sparse_format)
        dtau = scipy.sparse.diags(dtau, format=self.sparse_format)
        
        dQdp = self.build_dQdp_native(alpha= alpha, kappa=kappa, tau = tau, dkappa=dkappa, dtau=dtau)

        return dQdp
    
    
    
    def build_Q_native(self, kappa, tau, alpha):
        """
        Build the precision matrix for the latent process vector following method of Lindgren and Rue.
        
        This routine uses the native kappa, tau parameterisation.
        
        """

        C_tilde = self.build_Q_C_tilde()
        
        # C_tilde is diagonal so invert as follows
        inv_C_tilde = scipy.sparse.diags( 1.0 / C_tilde.diagonal(), format=self.sparse_format)
        
        G = self.build_Q_G()
        
        K = kappa * C_tilde * kappa + G
        
        if alpha % 2 == 0:
            
            Q = K * inv_C_tilde * K
            for i in range(4, alpha, 2):
                Q = K * inv_C_tilde * Q * inv_C_tilde * K
        else:
            Q = K
            for i in range(3, alpha, 2):
                Q = K * inv_C_tilde * Q * inv_C_tilde * K
        
        Q = tau * Q * tau
        
        return Q

    def build_dQdp_native(self, kappa= 1.0, tau = 1.0, alpha = 2, dkappa=None, dtau=None, level = None):
        """
        Build the precision matrix for the latent process vector following method of Lindgren and Rue.
        
        This routine uses the native kappa, tau parameterisation.
        
        """
        
        C_tilde = self.build_Q_C_tilde()
        
        # C_tilde is diagonal so invert as follows
        inv_C_tilde = scipy.sparse.diags( 1.0 / C_tilde.diagonal(), format=self.sparse_format)
        
        G = self.build_Q_G()
        
        K = kappa * C_tilde * kappa + G        
        dKdp = dkappa * C_tilde * kappa + kappa * C_tilde * dkappa
        
        if alpha % 2 == 0:
            
            Q = K * inv_C_tilde * K
            dQdp = dKdp * inv_C_tilde * K + K * inv_C_tilde * dKdp
            
            for i in range(4, alpha, 2):
                
                dQdp = dKdp * inv_C_tilde * Q * inv_C_tilde * K + \
                       K * inv_C_tilde * dQdp * inv_C_tilde * K + \
                       K * inv_C_tilde * Q * inv_C_tilde * dKdp
                Q = K * inv_C_tilde * Q * inv_C_tilde * K
        else:
            Q = K
            dQdp = dKdp
            for i in range(3, alpha, 2):
                
                dQdp = dKdp * inv_C_tilde * Q * inv_C_tilde * K + \
                       K * inv_C_tilde * dQdp * inv_C_tilde * K + \
                       K * inv_C_tilde * Q * inv_C_tilde * dKdp
                Q = K * inv_C_tilde * Q * inv_C_tilde * K
        
        dQdp = dtau * Q * tau + \
                tau * dQdp * tau + \
                tau * Q * dtau
        
        return dQdp

    def build_Q_C_tilde(self):
        """
        Build the Q_C matrix.
        
        """
                
        triangle_vertex_indices = self.triangulation.triangles[:,1:]
        
        dT = self.triangle_area(self.triangulation.points, triangle_vertex_indices)
        
        indices = triangle_vertex_indices.ravel()
        data = np.repeat(dT,3) / 3.0
        
        Q_C_tilde = scipy.sparse.coo_matrix( (data, (indices,indices) ) )        
        Q_C_tilde = Q_C_tilde.asformat(self.sparse_format)
        
        return Q_C_tilde
    
    def build_Q_G(self):
        """
        Build the Q_G matrix.
        
        """
        
        number_of_triangles = self.triangulation.triangles.shape[0]
        
        triangle_vertex_indices = self.triangulation.triangles[:,1:]
        dT = self.triangle_area(self.triangulation.points, triangle_vertex_indices)
        
        edges = self.triangle_edges(self.triangulation.points, triangle_vertex_indices).reshape( (number_of_triangles, 3, 3) )
        
        # Following einsum for edges_dot is equivalent to np.dot(edges[i,:,:], edges[i,:,:].T)
        # for all triangles index by i. Can shown to be with:
        # np.testing.assert_allclose(edges_dot[0,:,:], np.dot(edges[0,:,:], edges[0,:,:].T))       
        edges_dot = np.einsum('ijk,ilk->ijl', edges, edges)
        
        g = 1.0 / (4.0 * np.repeat(dT,9)) * edges_dot.ravel()
        
        data = g.ravel()
        
        indices_i = np.repeat( triangle_vertex_indices.ravel(), 3 )
        indices_j = np.tile( triangle_vertex_indices, (1,3) ).ravel()
        
        Q_G = scipy.sparse.coo_matrix( (data, (indices_i,indices_j) ) )
        Q_G = Q_G.asformat(self.sparse_format)
        
        return Q_G
    
    def build_Q_B(self):
        """
        See Lindgren and Rue. Not implemented.
        
        """
        
        raise NotImplementedError

    def mesh_at_level(self, level = 0, project_onto_sphere = True):
        """Returns a subdivided icosohedron mesh at a given level"""
        return MeshIcosahedronSubdivision.build(level, project_onto_sphere = project_onto_sphere)
    
    def neighbours_at_level(self, neighbourhood_level = 0, centre_index_at_level = None):
        """Return indices of vertices the enclosing triangle and surrounding triangles at a given level"""
        
        # extract triangles at the subdivision level of the desired set of neighbours defining the subregion of interest
        mesh_at_neighbourhood_level = self.mesh_at_level(neighbourhood_level)
        
        full_resolution_triangle_ids = self.triangulation.triangles[:,0]
        neighbour_resolution_triangles = mesh_at_neighbourhood_level.triangles
        
        # get indices to triangles in neighbourhood that include at least one vertex of the the centre triangle
        centre_triangle_vertex_indices = neighbour_resolution_triangles[centre_index_at_level,1:]
        in_neighbourhood = np.any( neighbour_resolution_triangles[:,1:] == centre_triangle_vertex_indices[0], axis = 1 ) | \
                           np.any( neighbour_resolution_triangles[:,1:] == centre_triangle_vertex_indices[1], axis = 1 ) | \
                           np.any( neighbour_resolution_triangles[:,1:] == centre_triangle_vertex_indices[2], axis = 1 )
        
        neighbour_triangle_ids = neighbour_resolution_triangles[in_neighbourhood,0]
        
        neighbourhood_triangle_indices = neighbour_triangle_ids >> MeshIcosahedron.BASE_HASH_SHIFT - np.uint64(2*neighbourhood_level)
        self_triangle_indices_at_neighbourhood_level = full_resolution_triangle_ids >> MeshIcosahedron.BASE_HASH_SHIFT - np.uint64(2*neighbourhood_level)
        
        # accumulate flags for full resolution triangle that are in the desired neighbourhood
        n_neighbours = len(neighbourhood_triangle_indices)

        in_neighbourhood_full_resolution = np.zeros(self.triangulation.triangles.shape[0], dtype = np.int16)
        for n in range(n_neighbours):
            in_neighbourhood_full_resolution += self_triangle_indices_at_neighbourhood_level == neighbourhood_triangle_indices[n]

        # find indices to full resolution triangles that are in the neighbourhood
        indices_to_neighbourhood = np.where(in_neighbourhood_full_resolution > 0)[0]
        
        return indices_to_neighbourhood
    
    def super_triangle_at_level(self, super_triangle_level = 0, super_triangle_index = None):
        """Return indices of vertices within a triangle at a given level"""
        
        # bit shift triangle indices to get indices of containing super triangles at the desired level
        full_resolution_triangle_ids = self.triangulation.triangles[:,0]        
        indices_at_super_level = full_resolution_triangle_ids >> MeshIcosahedron.BASE_HASH_SHIFT - np.uint64(2*super_triangle_level)
        
        # get indices to those triangles in the super triangle with the specified triangle index
        in_super_triangle = indices_at_super_level == super_triangle_index
        indices_in_super_triangle = np.where(in_super_triangle)[0]
        
        return indices_in_super_triangle
