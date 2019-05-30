"""

Sub-classes of SphereMeshSPDE for construction of the SPDE on subregions of the
nested triangulation. In future it may be desirable to update the base
SphereMeshSPDE class to include this functionality.

"""

import numpy as np
import scipy.sparse
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.mesh.mesh import MeshIcosahedron
from eustace.analysis.mesh.geometry import polar2d_to_cartesian

class SphereMeshView(SphereMeshSPDE):
    """
    
    Defines a "view" to a subregion of the subdivided icosohedral triangulation
    
    Sub-classes SphereMeshSPDE and manages indexing into the triangles and 
    vertices within the subregion. Note that precision matrices produced
    by the class are not the same as would be obtained by extraction of 
    block of the precision matrices produced by SphereMeshSPDE because the 
    boundary conditions on the SPDE implied by slicing are different. This is 
    apparent in the marginal variances at locations near to the boundary.
    
    """
    
    def __init__(self, level, active_triangle_indices, sparse_format='csr'):
        
        super(SphereMeshView, self).__init__(level = level, sparse_format = sparse_format)
        
        self.active_triangle_indices = np.sort(active_triangle_indices)
        self.active_vertex_indices = self.get_active_vertex_indices()
    
        #self.vertex_mapping_matrix = self.get_vertex_mappings( full_sphere.triangulation.vertices )
        #self.triangulation_mapping_matrix = self.get_triangulation_mappings( full_sphere.triangulation.triangles )
    
    @staticmethod
    def n_triangles_at_level(level):
        """Utility routine for obtaining the number of triangulations at a given level.
        
        Intended to be used to infer the number of views that make up a full sphere.
        
        """
        
        n_triangles = 20
        
        for this_level in range(1,level+1):        
            n_triangles = n_triangles * 4
            
        return n_triangles
    
    def get_active_vertex_indices(self):
        """Return indices to subset of nodes that are vertices of the active triangles"""
        
        active_vertex_indices = np.unique( self.triangulation.triangles[self.active_triangle_indices,1:].ravel() )
        
        return np.sort( active_vertex_indices )
        
    def get_vertex_mappings(self):
        """Mapping matrix between full mesh node indices to active mesh node indices"""
        
        number_of_vertices_total  = self.triangulation.points.shape[0]
        number_of_vertices_active = len(self.active_vertex_indices)
        
        input_indices  = self.active_vertex_indices
        output_indices = np.arange(number_of_vertices_active)
        
        
        mapping_matrix = scipy.sparse.dok_matrix( (number_of_vertices_active, number_of_vertices_total), dtype=np.int64 )
        mapping_matrix[output_indices,input_indices] = np.ones(number_of_vertices_active, dtype=np.int64)
                                                  
        return mapping_matrix.tocsr()
    
    def get_triangulation_mappings(self):
        """Mapping matrix from full mesh triangle indices to active mesh triangle indices"""
        
        number_of_triangles_total = self.triangulation.triangles.shape[0]
        number_of_triangles_active = len(self.active_triangle_indices)
        
        input_indices = self.active_triangle_indices
        output_indices = np.arange(number_of_triangles_active)
        
        mapping_matrix = scipy.sparse.dok_matrix((number_of_triangles_active, number_of_triangles_total), dtype=np.int64)        
        mapping_matrix[output_indices,input_indices] = np.ones(number_of_triangles_active, dtype=np.int64)
                                                  
        return mapping_matrix.tocsr()
        
    def n_latent_variables(self):
        """Return number of nodes in the mesh"""
        
        total_number = len( self.active_vertex_indices )
        
        return total_number
    
    def build_Q_C_tilde(self):
        """Build the Q_C matrix for a subregion"""
        
        # need indices to vertices/triangles in the view ordering, not the full domain ordering of vertices/triangles
        view_triangles = self.triangulation.triangles[self.active_triangle_indices,:]

        dT = self.triangle_area(self.triangulation.points, view_triangles[:,1:])
        
        invalid_index = -self.triangulation.points.shape[0] # initialise to invalid index outside of index range
        vertex_index_map = np.zeros( self.triangulation.points.shape[0], dtype = np.int64 ) + invalid_index
        vertex_index_map[self.active_vertex_indices] = np.arange( len(self.active_vertex_indices) )
        
        indices = vertex_index_map[view_triangles[:,1:].ravel()]
        
        data = np.repeat(dT,3) / 3.0
        
        Q_C_tilde = scipy.sparse.coo_matrix( (data, (indices,indices) ) )
        Q_C_tilde = Q_C_tilde.asformat(self.sparse_format)
        
        return Q_C_tilde
    
    def build_Q_G(self):
        """Build the Q_G matrix for a subregion"""
        
        view_triangles = self.triangulation.triangles[self.active_triangle_indices,:]
        dT = self.triangle_area(self.triangulation.points, view_triangles[:,1:])
        
        number_of_triangles = view_triangles.shape[0]
        
        edges = self.triangle_edges(self.triangulation.points, view_triangles[:,1:]).reshape( (number_of_triangles, 3, 3) )

        invalid_index = -self.triangulation.points.shape[0] # initialise to invalid index outside of index range
        vertex_index_map = np.zeros( self.triangulation.points.shape[0], dtype = np.int64 ) + invalid_index
        vertex_index_map[self.active_vertex_indices] = np.arange( len(self.active_vertex_indices) )
        
        variable_indices0 = vertex_index_map[view_triangles[:,1]]
        variable_indices1 = vertex_index_map[view_triangles[:,2]]
        variable_indices2 = vertex_index_map[view_triangles[:,3]]
        
        triangle_vertex_indices = np.vstack( (variable_indices0, variable_indices1, variable_indices2) ).T
        
        
        # Following einsum for edges_dot is equivalent to np.dot(edges[i,:,:], edges[i,:,:].T)
        # for all triangles index by i. Can be shown to be with:
        # np.testing.assert_allclose(edges_dot[0,:,:], np.dot(edges[0,:,:], edges[0,:,:].T))       
        edges_dot = np.einsum('ijk,ilk->ijl', edges, edges)
        
        g = 1.0 / (4.0 * np.repeat(dT,9)) * edges_dot.ravel()
        
        data = g.ravel()
        
        indices_i = np.repeat( triangle_vertex_indices.ravel(), 3 )
        indices_j = np.tile( triangle_vertex_indices, (1,3) ).ravel()
        
        Q_G = scipy.sparse.coo_matrix( (data, (indices_i,indices_j) ) )
        Q_G = Q_G.asformat(self.sparse_format)
        
        return Q_G
    
    def build_A(self, polar_coordinates):
        """
        Takes a set of points at which to evaluate basis functions (at triangulation level if specified) and returns the A matrix of basis weights.
        
        """
        
        # Convert to cartesian coordinates on unit sphere
        points = polar2d_to_cartesian(polar_coordinates)
        
        # Evaluate weights by computing barycentres
        triangle_id, basis_weights = self.triangle_locator.locate_vector_hierarchical(points, levels = self.triangulation.level)
        
        triangle_index = triangle_id >> MeshIcosahedron.BASE_HASH_SHIFT - np.uint64(2*self.triangulation.level)
        
        # filter locations that are within the domain of the active triangles
        in_domain = np.in1d(triangle_index, self.active_triangle_indices)

        # Put into matrix
        n_bases_full = self.triangulation.points.shape[0]
        n_locations = polar_coordinates.shape[0]
        point_indices  = np.repeat( np.arange(n_locations)[in_domain], 3).ravel() 
        vertex_indices = self.triangulation.triangles[triangle_index[in_domain],1:].ravel()
        A = scipy.sparse.dok_matrix((n_locations, n_bases_full))
        A[point_indices, vertex_indices] = basis_weights[in_domain,:].ravel()
        
        # extract the columns corresponding to the active vertices in the local view
        A = A.tocsc()[:,self.active_vertex_indices]
        
        return A.asformat(self.sparse_format)
        
    def in_domain(self, polar_coordinates):
        """Flag whether input locations lie within the view's domain"""
        
        # Convert to cartesian coordinates on unit sphere
        points = polar2d_to_cartesian(polar_coordinates)
        
        # Evaluate weights by computing barycentres
        triangle_id, basis_weights = self.triangle_locator.locate_vector_hierarchical(points, levels = self.triangulation.level)
        
        triangle_index = triangle_id >> MeshIcosahedron.BASE_HASH_SHIFT - np.uint64(2*self.triangulation.level)
        
        # filter locations that are within the domain of the active triangles
        in_domain = np.in1d(triangle_index, self.active_triangle_indices)
        
        in_flags = np.zeros(points.shape[0], dtype = np.bool)
        in_flags[in_domain] = True
        
        return in_flags
    
class SphereMeshViewGlobal(SphereMeshView):
    """
    
    A "view" of the full domain of a SphereMeshSPDE
    
    """
    
    def __init__(self, level, sparse_format='csr'):
        
        full_sphere = SphereMeshSPDE(level)
        
        active_triangle_indices = np.arange(full_sphere.triangulation.triangles.shape[0])
        
        super(SphereMeshViewGlobal, self).__init__(level, active_triangle_indices, sparse_format=sparse_format)

    def merge_local_parameterisations(self, local_spde_view_list, local_hyperparameter_list, merge_method = 'direct_average'):
        """Generate merged parameterisation for non-stationary model from local parameterisations.
        
        Merges hyperparameter design matrices that map hyperparameters onto each latent variable
        and merges the hyperparameters into a single hyperparameter vector. Where multiple
        overlapping local models contribute to the same latent variable the projected 
        hyperparameters are effectively averaged with equal weight across the contributing
        local models.
        
        TODO: Revise this method to incrementally compute the non-expanded merged parameterisation
        
        """
        
        # intialise accumulators for local rho, sigma design matrices and local hyperparameter representations
        sigma_design_accumulator = []
        rho_design_accumulator = []
        theta_accumulator = []
        n_latent_variables = self.n_latent_variables()
        
        # initialise counter for number of local contributors for each latent variable
        sigma_contribution_counter = np.zeros( n_latent_variables )
        rho_contribution_counter   = np.zeros( n_latent_variables )
        
        for (local_spde, hyperparameters) in zip(local_spde_view_list, local_hyperparameter_list):
            
            # map the local design matrices for local representation to the global set of vertices
            vertex_map = local_spde.get_vertex_mappings().T
            
            local_sigma_design = vertex_map.dot( local_spde.sigma_design_matrix_stationary() )
            sigma_contribution_counter += np.int64( local_sigma_design.getnnz(axis = 1)  > 0 )
            
            local_rho_design = vertex_map.dot( local_spde.rho_design_matrix_stationary() )
            rho_contribution_counter += np.int64( local_rho_design.getnnz(axis = 1)  > 0 )
            
            # accumulate the design matrices
            sigma_design_accumulator.append( local_sigma_design )
            rho_design_accumulator.append( local_rho_design )
            theta_accumulator.append( hyperparameters )
        
        # concetenate the local hyperparameters into one hyperparameter vector
        merged_hyperparameter_values = np.concatenate( theta_accumulator )
        
        sigma_design = scipy.sparse.hstack(sigma_design_accumulator)
        rho_design = scipy.sparse.hstack(rho_design_accumulator) 
        if merge_method == 'direct_average':
            # averages the hyperparameters in their provided form, e.g. averages the log_sigma values.
            
            # merge the sigma local design matrices and normalise contributions at vertices with more than one local contributor
            sigma_normaliser = scipy.sparse.diags( 1.0 / sigma_contribution_counter, 0 )
            sigma_design = sigma_normaliser * sigma_design
        
            # merge the rho local design matrices and normalise contributions at vertices with more than one local contributor
            rho_normaliser = scipy.sparse.diags( 1.0 / rho_contribution_counter, 0 )
            rho_design = sigma_normaliser * rho_design
            
            # concetenate the local hyperparameters into one hyperparameter vector
            merged_hyperparameter_values = np.concatenate( theta_accumulator )
            
            return merged_hyperparameter_values, sigma_design, rho_design
            
        elif merge_method == 'exp_average':
            # returns an expanded parameterisation which averages the exponential of the local values, e.g. averages the sigma values rather than log_sigma.
            
            # merge the sigma local design matrices and normalise contributions at vertices with more than one local contributor
            sigma_normaliser = scipy.sparse.diags( 1.0 / sigma_contribution_counter, 0 )
            sigma_design = sigma_normaliser * sigma_design
        
            log_sigma = np.log( sigma_design.dot( np.exp(merged_hyperparameter_values) ) )
        
            # merge the rho local design matrices and normalise contributions at vertices with more than one local contributor
            rho_normaliser = scipy.sparse.diags( 1.0 / rho_contribution_counter, 0 )
            rho_design = sigma_normaliser * rho_design
            
            log_rho = np.log( rho_design.dot( np.exp(merged_hyperparameter_values) ) )
            
            expanded_sigma_design = scipy.sparse.hstack( [scipy.sparse.eye(n_latent_variables),
                                                          scipy.sparse.dok_matrix((n_latent_variables, n_latent_variables))] )
                                                          
            expanded_rho_design = scipy.sparse.hstack( [scipy.sparse.dok_matrix((n_latent_variables, n_latent_variables)),
                                                        scipy.sparse.eye(n_latent_variables)] )
            
            return np.concatenate( [log_sigma, log_rho] ), expanded_sigma_design, expanded_rho_design
            
        else:
            
            raise ValueError('Unrecognised merge method: {}'.format(merge_method)) 
    
    @staticmethod
    def accumulate_local_parameterisations(sigma_accumulator, rho_accumulator, contribution_counter, local_spde, local_hyperparameters):
        """Accumulate local parameterisations towards building a global nonstationary parameterisation.
        
        Output should be passed to finalise_local_parameterisations for final merged non-stationary parameterisation.
        
        """

        # map the local design matrices for local representation to the global set of vertices
        global_vertex_map = local_spde.get_vertex_mappings().T
        contribution_indicator = np.int64( global_vertex_map.getnnz(axis = 1)  > 0 )
        
        local_sigma_design = local_spde.sigma_design_matrix_stationary()
        local_rho_design = local_spde.rho_design_matrix_stationary()


        #sigma_contribution = np.ones(len(sigma_contribution))
        #rho_contribution = np.ones(len(rho_contribution))

        sigma = global_vertex_map.dot( np.exp( local_sigma_design.dot( local_hyperparameters ) ) )
        rho = global_vertex_map.dot( np.exp( local_rho_design.dot( local_hyperparameters ) ) )
        
        
        if sigma_accumulator is None or rho_accumulator is None or contribution_counter is None:
            sigma_accumulator = sigma
            rho_accumulator = rho            
            contribution_counter = contribution_indicator
        else:
            sigma_accumulator += sigma
            rho_accumulator += rho            
            contribution_counter += contribution_indicator

        return sigma_accumulator, rho_accumulator, contribution_counter

    @staticmethod
    def finalise_local_parameterisation_sigma_rho(sigma_accumulator, rho_accumulator, contribution_counter):
        """Apply normaliser global nonstationary parameterisation with (log_sigma, log_rho) parameterisation.
        
        """
        
        log_sigma = np.log( sigma_accumulator / contribution_counter )
        log_rho   = np.log( rho_accumulator / contribution_counter )
    
        return log_sigma, log_rho
        
    @staticmethod
    def finalise_local_parameterisation_q(sigma_accumulator, rho_accumulator, contribution_counter):
        """Apply normaliser global nonstationary parameterisation with (q_parameter, log_sigma_design, log_rho_design) parameterisation."""
        
        n_sigma = len(sigma_accumulator)
        n_rho = len(rho_accumulator)
        
        q_parameters = np.concatenate( [ np.log( sigma_accumulator / contribution_counter ),
                                            np.log( rho_accumulator / contribution_counter )] )
        
        log_sigma_design = scipy.sparse.hstack( [scipy.sparse.eye( n_sigma ),
                                                 scipy.sparse.dok_matrix((n_rho, n_rho))] )
                                                          
        log_rho_design = scipy.sparse.hstack( [scipy.sparse.dok_matrix((n_sigma, n_sigma)),
                                               scipy.sparse.eye(n_rho)] )
    
        return q_parameters, log_sigma_design, log_rho_design
        
        

class SphereMeshViewLocal(SphereMeshView):
    """
    
    A "view" of a local domain of a SphereMeshSPDE
    
    """
    
    def __init__(self, level, neighbourhood_level, centre_index_at_level, sparse_format='csr'):
        
        full_sphere = SphereMeshSPDE(level)
        
        active_triangle_indices = full_sphere.neighbours_at_level( neighbourhood_level = neighbourhood_level, centre_index_at_level = centre_index_at_level)
        
        super(SphereMeshViewLocal, self).__init__(level, active_triangle_indices, sparse_format=sparse_format)    

class SphereMeshViewSuperTriangle(SphereMeshView):
    """
    
    A "view" of a SphereMeshSPDE defined by a low resolution super triangle
    
    """
    
    def __init__(self, level, super_triangle_level, super_triangle_index, sparse_format='csr'):
        
        full_sphere = SphereMeshSPDE(level)
        
        active_triangle_indices = full_sphere.super_triangle_at_level( super_triangle_level = super_triangle_level, super_triangle_index = super_triangle_index)
        
        super(SphereMeshViewSuperTriangle, self).__init__(level, active_triangle_indices, sparse_format=sparse_format)
        
   