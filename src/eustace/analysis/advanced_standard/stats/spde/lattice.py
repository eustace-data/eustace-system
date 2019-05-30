'''
Created on Apr 11, 2017

@author: hadcc
'''

import numpy as np
import scipy.sparse

class BasisFunction(object):
    
    """
    
    Class for representation of a compactly supported ND basis function centered at the origin.
    
    """
    
    def __init__(self, basis_span=np.array([1.0, 1.0])):
        """
        
        Kwargs:
        
        * basis_span (ndarray):
            A 1D array of distances beyond which the basis function is zero valued.
            with elements ordered by dimension.
        
        """
        
        self.basis_span = basis_span

    def evaluate(self, eval_distances):
        """
        Evaluate the basis function at the specified distances from the basis function origin.
        
        Args:
        
        * eval_distances (ndarray):
            An ND array of distances of shape [number_of_observations, number_of_dimensions].
            Each row takes the value of the distances in each dimension of the location for 
            evaluation away from the basis function origin.
        
        """
        
        scaled_distances = np.linalg.norm(eval_distances / self.basis_span, axis = 1)
        
        scaled_distances[scaled_distances>1.0] = 1.0
        
        return 1.0 - scaled_distances

class WendlandC4Basis(BasisFunction):
    """A class for representing Wendland C4 basis functions."""
    
    def __init__(self, **kwargs):
        
        super(WendlandC4Basis, self).__init__(**kwargs)
        
    def evaluate(self, eval_distances):
        
        scaled_distances = np.linalg.norm(eval_distances / self.basis_span, axis = 1)
    
        scaled_distances[scaled_distances>1.0] = 1.0
        
        function_values = (1.0 - scaled_distances)**6 * (35.0*scaled_distances**2 +18.0*scaled_distances +3.0) / 3

        return function_values

class PiecewiseLinearBasis(BasisFunction):
    """A class for representing Piecewise Linear basis functions."""
    
    def __init__(self, *args, **kwargs):
        
        super(PiecewiseLinearBasis, self).__init__(*args, **kwargs)
        
    def evaluate(self, eval_distances):

        normalied_distances = np.abs( eval_distances / self.basis_span )
    
        normalied_distances[normalied_distances>1.0]  = 1.0
        #normalied_distances[normalied_distances<-1.0] = -1.0
        
        function_values = np.prod(1.0 - normalied_distances, axis = 1)

        return function_values        

class LatticeMesh(object):
    
    """
    
    Class for representation of a meshed 1D or 2D domain on regularly spaced lattice.
    
    """
    
    def __init__(self, dimension_specification = [(0., 10., 11),
                                                  (20., 30., 11)], basis_function = BasisFunction(), overlap_factor = None, wrap_dimensions = None):
        """
        
        A class for a domain represented as an ND lattice.
        
        Args:
        
        * dimension_specification (list):
            A list of tuples formated as (min_value, max_value, n_nodes) where
            each tuple describes the spacing of node locations in one dimension. 
            
        * basis_function (BasisFunction):
            A BasisFunction object used to compute values of basis functions
            centered on the nodes of the lattice. One BasisFunction object is used 
            for all basis functions on the lattice.
        
        * wrap_dimensions:
            List of booleans. Speficies cyclical dimensions where the min_value and 
            max_value in dimension_specification  represent the same location in that 
            dimension but a full period apart. The domain spanned by the process is
            [min_value, max_value) with nodes evenly spaced between min_value and
            max_value - (max_value - min_value) / (n_nodes).
            
        TODO (CPM) Add capability to wrap around dimensions for locate_nearby_basis_indices 
            and build_A.
        
        """
        
        self.basis_function = basis_function
        
        self.n_dims = len(dimension_specification)
        self.wrap_dimensions = wrap_dimensions
        
        self.domain_bounds = np.array([np.float64(this_dim[0:2]) for this_dim in dimension_specification]).T

        self.number_of_nodes_per_dimension = np.array([this_dim[2] for this_dim in dimension_specification])
        self.node_bounds = self.domain_bounds

        if self.wrap_dimensions is None:
            # Nodes evenly spaced between start and end points in dimension_specification. Nodes at 
            # the start and end points.
            self.node_spacing = (self.domain_bounds[1,:] - self.domain_bounds[0,:]) / (np.float64(self.number_of_nodes_per_dimension) - 1.0)
        else:
            # Same for dimensions that do not wrap. For wrapped dimensions the nodes are space between 
            # the start point in dimension_specification and the end_point - node_separation. The
            # boundaries in dimension_specification effectively represent the same location in that 
            # dimension but a full period apart.            
            self.node_spacing = (self.domain_bounds[1,:] - self.domain_bounds[0,:]) / (np.float64(self.number_of_nodes_per_dimension + self.wrap_dimensions) - 1.0)
            self.node_bounds[1,:] = self.node_bounds[1,:] - self.wrap_dimensions * self.node_spacing
        
        
        axis_coordinates = [np.linspace(this_dim[0], this_dim[1], n_nodes) for this_dim, n_nodes in zip(self.node_bounds.T, self.number_of_nodes_per_dimension )  ]
        self.axis_coordinates = np.array(axis_coordinates).T
        
        if len(axis_coordinates) > 1:
            meshgrid_nodes = np.meshgrid(*axis_coordinates, indexing='ij')
            self.nodes = np.array(meshgrid_nodes).T.reshape(-1,self.n_dims)
        else:
            self.nodes = self.axis_coordinates

        self.n_nodes = self.nodes.shape[0]
        
        if overlap_factor is not None:
            self.basis_function.basis_span = self.node_spacing * overlap_factor


    def locate_nearby_basis_indices(self, locations):
        """Compute the indices for non-zero valued basis functions and the specified locations.
        
        Output is formated as an ndarray with rows valued as:
        
        [obs_indices, [dim0_indices, ...]]
        
        """ 
        
        n_locations = locations.shape[0]
        
        node_origin     = self.domain_bounds[0,:]
        node_separation = self.node_spacing
        
        # find closest index to node
        nearest_nodes   = np.floor( (locations - np.tile(node_origin-0.5*node_separation, (n_locations, 1))) / np.tile(node_separation, (n_locations, 1)) ).astype(np.int64)

        # find number of indices to search in each dimension
        node_search_range = self.basis_function.basis_span / node_separation + 0.5 # node search range is basis span in decimal number of nodes + 0.5 to account for observations being up to half of the node separation away from the closest node
        node_index_search_range = np.atleast_1d( np.int64( np.floor( node_search_range ) ) )
        
        # compute index shifts in each dimension
        slices = [slice(-nind, nind+1) for nind in node_index_search_range] # Is the +1 needed here? Check this.
        search_cube_inds = np.mgrid[slices]
        search_cube_inds = search_cube_inds.reshape(search_cube_inds.shape[0],-1)
        
        n_shifts = search_cube_inds.shape[1]
        
        # apply shifts to compute indices of nearby nodes for each observation location and generate corresponding indices for the repeated observation locations
        node_indices = (nearest_nodes[:,None,:] + search_cube_inds.T).reshape(n_locations*n_shifts,self.n_dims, order='C')
        obs_indices  = np.arange(n_locations).repeat(n_shifts)
        
        # wrap indexing for circular dimensions
        if self.wrap_dimensions is not None:
            for dim_index in range(self.n_dims):
                if self.wrap_dimensions[dim_index]:
                    node_indices[:,dim_index] = node_indices[:,dim_index] % (self.number_of_nodes_per_dimension[dim_index])
            
            
        # remove entries with indices less than zero or greater than self.number_of_nodes_per_dimension
        lower_bound_test = np.all(node_indices >= 0, axis=1)
        upper_bound_test = np.all(node_indices / np.tile( (self.number_of_nodes_per_dimension-1.0), (node_indices.shape[0], 1)) <= 1, axis=1 )
        indices_in_bounds = np.logical_and(lower_bound_test, upper_bound_test)
        
        obs_indices  = obs_indices[indices_in_bounds]
        node_indices = node_indices[indices_in_bounds,:]
        
        # convert the N dimensional node indices into 1D indices 
        #node_indices_1d =  np.ravel_multi_index(node_indices[:,[2,0,1]].T, self.number_of_nodes_per_dimension)
        node_indices_1d =  np.ravel_multi_index(node_indices.T, self.number_of_nodes_per_dimension, order = 'F')

        return obs_indices, node_indices_1d
        
    def build_A(self, locations, sparse_format = 'csc'):    
        """Construct sparse design matrix by evaluating basis functions for the specified locations.

        Args:
        
        * locations (ndarray):
            An ND array of locations of shape [number_of_observations, number_of_dimensions].
            Each row takes the value of the locations in each dimension at which the design
            matrix is to be evaluated.
        
        Kwargs:
        
        * sparse_format (str):
            String description of sparse matrix format for output design
            matrix, e.g. 'csc', 'csr', 'dok', 'coo'.
            
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.

        """
        
        if locations.ndim == 1:
            # catches 1d vector pass in case of observations of a 1D process
            locations = np.atleast_2d(locations).T
        
        n_locations = locations.shape[0]
        n_bases = self.nodes.shape[0]

        # find all required non-zero basis functions for each observations
        obs_indices, node_indices = self.locate_nearby_basis_indices(locations)
        
        # form a vector of all required basis evaluations - multi observation, multi basis
        distances = locations[obs_indices,...] - self.nodes[node_indices,...]
        
        # account for basis functions overlapping wrapped dimensions
        if self.wrap_dimensions is not None:
            for dimension in range(self.n_dims):
                if self.wrap_dimensions[dimension]:
                    #print "distances pre wrap"
                    #print distances
                    # find the full range of circular dimension
                    dimension_span = self.number_of_nodes_per_dimension[dimension] * self.node_spacing[dimension]
                    
                    distances[:,dimension] = distances[:,dimension] % dimension_span
                    #print distances
                    distant_mappings = distances[:,dimension] > dimension_span / 2.0
                    distances[distant_mappings,dimension] = distances[distant_mappings,dimension] % -dimension_span
                    #print distances
                    ## compute distances as distnaces to nodes within the closest half circle of the circular dimension
                    #distant_mappings = distances[:,dimension] > dimension_span / 2.0
                    
                    #distances[distant_mappings,dimension] = distances[distant_mappings,dimension] % -dimension_span
                    
                    #print "distances post wrap + "
                    #print distances
                    
                    #distant_mappings = distances[:,dimension] < -dimension_span / 2.0
                    
                    #distances[distant_mappings,dimension] = distances[distant_mappings,dimension] % dimension_span
                    #print "distances post wrap - "
                    #print distances
        
        a_values = self.basis_function.evaluate(distances)
        
        # construct the A matrix
        A_matrix = scipy.sparse.dok_matrix((n_locations, n_bases))
        A_matrix[obs_indices, node_indices] = a_values
        
        return A_matrix.asformat(sparse_format)

# import advanced_standard.linalg.sparse
    
class LatticeSPDE():
    """
    Class for representation of random fields on triangulated spheres.
        
    """
    
    def __init__(self, lattice = None):
        
        self.lattice = lattice
    
    @staticmethod
    def construct(dimension_specification = None, basis_function = None, overlap_factor = None, wrap_dimensions = None):
    
        lattice = LatticeMesh(dimension_specification = dimension_specification, basis_function = basis_function, overlap_factor = overlap_factor, wrap_dimensions = wrap_dimensions)
        spde_model = LatticeSPDE(lattice)
    
        return spde_model
    
    ###########################################################################
    # Tools for building Q
    ###########################################################################

    def n_latent_variables(self):
        """Return the number of latent variables in the SPDE model.
        
        This is equal to the number of nodes in the lattice.
        
        """
        return self.lattice.n_nodes
    
    def n_bases(self):
        
        return self.n_latent_variables()
    
    def build_A(self, locations, sparse_format = 'csc', Q_matrix = None, fQ_matrix = None):
        """Construct the design matrix.
        
        The design matrix may be normalised to ensure unit marginal variance each each location,
        compensating for changes in marginal variance due to artefacts of the basis function 
        representation. Note that this will result A * Q * A.T having unit marginal variance,
        compensating for any intended non-unit variance behaviour. Use of this feature while
        optimising variance adjusting parameters of Q is not advised!  
        
        Args:
        
            * locations:
                Locations at which to evaluate the design matrix.
            
            * sparse_format:
                Sparse matrix format for output.
                
            * Q_matrix:
                A sparse precision matrix.  If provided then the design matrix
                will be scaled for unit marginal variance at all locations.
                
            * fQ_matrix:
                A pre-computed cholesky decomposition of the Q_matrix. If provided
                then this is used rather than Q_matrix.
                
        """
        
        A_matrix =  self.lattice.build_A(locations, sparse_format = sparse_format)
        
        if fQ_matrix is not None:

            raise NotImplementedError
            # Taken from original prototype code but not yet implemented here due to reliance on cholmod:
            # normaliser = 1.0 / np.sqrt( advanced_standard.linalg.sparse.marginal_variance(None, fQ = fQ_matrix, A = A_matrix) )
            # A_matrix = scipy.sparse.diags(normaliser, format = sparse_format) * A_matrix

        elif Q_matrix is not None:

            raise NotImplementedError
            # Taken from original prototype code but not yet implemented here due to reliance on cholmod:
            # normaliser = 1.0 / np.sqrt( advanced_standard.linalg.sparse.marginal_variance(Q_matrix, A = A_matrix) )
            # A_matrix = scipy.sparse.diags(normaliser, format = sparse_format) * A_matrix              
            
        return A_matrix
    

    def sigma_design_matrix_stationary(self):
        """Get sigma matrix for stationary process."""

        number_of_latent_variables = self.n_latent_variables()        
        sigma_design_matrix = np.zeros((number_of_latent_variables, 2))
        sigma_design_matrix[:,0] = 1.0
        return sigma_design_matrix

    def rho_design_matrix_stationary(self):
        """Get rho matrix for stationary process."""
        
        number_of_latent_variables = self.n_latent_variables()        
        rho_design_matrix = np.zeros((number_of_latent_variables, 2))
        rho_design_matrix[:,1] = 1.0
        return rho_design_matrix

    def build_Q_stationary(self, log_sigma, log_rho, alpha, H, sparse_format):

        Q_parameters = np.array([[ log_sigma, log_rho ]]).T           
        return self.build_Q(Q_parameters, self.sigma_design_matrix_stationary(), self.rho_design_matrix_stationary(), alpha, H, sparse_format)

    def build_dQdp_stationary(self, log_sigma, log_rho, alpha, H, parameter_index, sparse_format):

        Q_parameters = np.array([[ log_sigma, log_rho ]]).T           
        return self.build_dQdp(Q_parameters, self.sigma_design_matrix_stationary(), self.rho_design_matrix_stationary(), alpha, H, parameter_index, sparse_format)

    
    def build_Q(self, Q_parameters, sigma_design_matrix, rho_design_matrix, alpha, H, sparse_format):
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

        number_of_latent_variables = self.n_latent_variables()

        # Stationary constant sigma as default
        if Q_parameters is None:
            Q_parameters = np.array([[0.0, 0.0]]).T
        if sigma_design_matrix is None:
            sigma_design_matrix = np.zeros((number_of_latent_variables, 2))
            sigma_design_matrix[:,0] = 1.0
        if rho_design_matrix is None:
            rho_design_matrix = np.zeros((number_of_latent_variables, 2))
            rho_design_matrix[:,1] = 1.0
           
        
        # convert sigma, rho parameterisation into kappa, tau
        d = self.lattice.n_dims   # domain dimension
        nu = alpha - d / 2.0
        
        # Convert provided design matrices into appropriate design matrices for the internal (tau, kappa) representation.        
        b0_sigma = 0.0 # Equivalent to a unity base scale factor on sigma
        b0_rho   = 0.0 # Equivalent to a unity base scale factor on rho
        
        b0_kappa = np.log( 8.0 * nu ) / 2.0 - b0_rho
        b0_tau   = 0.5 * np.log( gamma(nu) / (gamma(alpha)* (4*np.pi)**(0.5*d) ) ) - b0_sigma - nu * b0_kappa
        
        kappa_design_matrix = -rho_design_matrix
        tau_design_matrix   = -sigma_design_matrix - nu * kappa_design_matrix
        
        # compute non-stationary kappa and tau vectors
        kappa = np.exp( b0_kappa + np.dot(kappa_design_matrix, Q_parameters) ).ravel()
        tau   = np.exp( b0_tau + np.dot(tau_design_matrix, Q_parameters) ).ravel()
        
        #print gamma(nu) / (gamma(alpha)* (4*np.pi)**(0.5*d)*kappa**(2*nu)*tau**2 )
        
        kappa = scipy.sparse.diags(kappa, format = sparse_format)
        tau = scipy.sparse.diags(tau, format = sparse_format)
        
        Q = self.build_Q_native(kappa = kappa, alpha= alpha, tau = tau, H = H, sparse_format = sparse_format)

        #import advanced_standard.linalg.sparse
        #print advanced_standard.linalg.sparse.marginal_variance(Q)
        #print advanced_standard.linalg.sparse.solve_A(Q, scipy.sparse.eye(number_of_latent_variables)).diagonal()
        
        return Q
    
    def build_dQdp(self, Q_parameters, sigma_design_matrix, rho_design_matrix, alpha=2, H = None, parameter_index = None, sparse_format = 'csc'):
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

        number_of_latent_variables = self.n_latent_variables()

        # convert sigma, rho parameterisation into kappa, tau
        d = self.lattice.n_dims   # domain dimension for manifold
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
        kappa = np.exp( b0_kappa + np.dot(kappa_design_matrix, Q_parameters) ).ravel()
        tau   = np.exp( b0_tau + np.dot(tau_design_matrix, Q_parameters) ).ravel()
    
        kappa = scipy.sparse.diags(kappa, format = sparse_format)
        tau = scipy.sparse.diags(tau, format = sparse_format)
        
        dkappa = scipy.sparse.diags(kappa_design_matrix[:,parameter_index] * kappa)
        dtau = scipy.sparse.diags(tau_design_matrix[:,parameter_index] * tau)
        
        dQdp = self.build_dQdp_native(alpha= alpha, kappa=kappa, tau = tau, dkappa=dkappa, dtau=dtau, H = H, sparse_format = sparse_format)        

        return dQdp
    
    def build_Q_native(self, kappa=0.5, tau = 1.0, alpha = 2, H = None, sparse_format = 'csc'):
        """Construct the precision matrix for the LatticeSPDE.
        
        Args:
        
        * kappa:
            Either a single value or a sparse matrix with values of kappa at the
            lattice node locations on the diagonal. The kappa parameter controls
            both the smoothness and the variance of the process.
            
        * tau:
            Either a single value or a sparse matrix with values of tau at the
            lattice node locations on the diagonal. The tau parameter controls
            the variance of the process. Note that tau and alpha will also influence
            the variance in this native parameterisation of the SPDE.
            
        * alpha (int):
            A single integer value. Controls the smoothness of the SPDE.
            The smoothness of the random process increase with increasing alpha.
            The following values of alpha are permitted:
                alpha = 1, an AR1 process
                alpha = 2, a convolution of two AR1 processes
        
        * H:
            A range scaling factor allowing that has the effect of adjusting the
            decorreation rate of the process. By default H takes a value of one.
            
            If H is provided as an array of length equal to the number of dimensions
            then H effectively scales the distances between nodes of the triangulation 
            allowing for limited anisotropy.
        
        Kwargs:
        
        * sparse_format (str):
            String description of sparse matrix format for output design
            matrix, e.g. 'csc', 'csr', 'dok', 'coo'.
            
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
        
        """
        
        if H is None:
            H = 1.0#self.lattice.node_spacing**2
        
        dim_scales = H / self.lattice.node_spacing**2
        
        Q_neighbour = scipy.sparse.dok_matrix( ( self.n_latent_variables(), self.n_latent_variables() ) )
        Q_diagonal = kappa * scipy.sparse.eye(self.n_latent_variables()) * kappa + scipy.sparse.eye(self.n_latent_variables()) * (np.sum(2.0 * dim_scales))
        
        for d in range(self.lattice.n_dims):
            
            # generate slices corresponding to indexing in each dimension
            slices_i = [slice(0, dim_length, None) for dim_length in self.lattice.number_of_nodes_per_dimension]
            slices_j = [slice(0, dim_length, None) for dim_length in self.lattice.number_of_nodes_per_dimension]
            
            # overwrite slice for dimension d so that i indices are offset in dimensions d by -1 in comparison to j indices
            if self.lattice.wrap_dimensions is not None and self.lattice.wrap_dimensions[d] is True:
                # slice for sequential indexing with additional wrapping so the last index maps to the first
                slices_i[d] = slice(-1, self.lattice.number_of_nodes_per_dimension[d]-1)
                slices_j[d] = slice(0, self.lattice.number_of_nodes_per_dimension[d]) 
            else:
                # only link sequential indexes
                slices_i[d] = slice(0, self.lattice.number_of_nodes_per_dimension[d]-1)
                slices_j[d] = slice(1, self.lattice.number_of_nodes_per_dimension[d]) 
            
            # construct mesh grid containing all indices to be used identifying neighours and unwrap the indices
            inds_i = np.mgrid[slices_i].T.reshape(-1,self.lattice.n_dims)
            inds_j = np.mgrid[slices_j].T.reshape(-1,self.lattice.n_dims)
            
            # normalise any negative indices to positive range
            inds_i = inds_i % self.lattice.number_of_nodes_per_dimension
            inds_j = inds_j % self.lattice.number_of_nodes_per_dimension
            
            # covert the n-dimensional indices into one dimensional indices that correspond to the indexing of Q
            inds_i = np.ravel_multi_index(inds_i.T, self.lattice.number_of_nodes_per_dimension)
            inds_j = np.ravel_multi_index(inds_j.T, self.lattice.number_of_nodes_per_dimension)
            
            # use the indexing to fill in linked nodes in Q
            Q_neighbour[inds_i, inds_j] = - dim_scales[d]
            Q_neighbour[inds_j, inds_i] = - dim_scales[d]
        
        Q = Q_diagonal + Q_neighbour
        
        if alpha == 2:
            Q = Q.T * Q
        elif alpha > 2:
            raise NotImplementedError('Values of alpha > 2 are not yet supported.')
        
        Q = tau * Q * tau

        s = np.prod(self.lattice.node_spacing)
        Q = s * Q

        return Q
    
    def build_dQdp_native(self, kappa= 1.0, tau = 1.0, alpha = 2, dkappa=None, dtau=None, level = None, H = None, sparse_format = 'csc'):
        """
        Build the precision matrix for the latent process vector following method of Lindgren and Rue.
        
        This routine uses the native kappa, tau parameterisation.
        
        """
        
        if H is None:
            H = H = 1.0#self.lattice.node_spacing**2
        
        dim_scales = H / self.lattice.node_spacing**2
        
        Q_neighbour = scipy.sparse.dok_matrix( ( self.n_latent_variables(), self.n_latent_variables() ) )
        Q_diagonal = kappa * scipy.sparse.eye(self.n_latent_variables()) * kappa + scipy.sparse.eye(self.n_latent_variables()) * (np.sum(2.0 * dim_scales))
        
        for d in range(self.lattice.n_dims):
            
            slices_i = [slice(0, dim_length, None) for dim_length in self.lattice.number_of_nodes_per_dimension]
            slices_i[d] = slice(0, self.lattice.number_of_nodes_per_dimension[d]-1)
            
            slices_j = [slice(0, dim_length, None) for dim_length in self.lattice.number_of_nodes_per_dimension]
            slices_j[d] = slice(1, self.lattice.number_of_nodes_per_dimension[d]) 
            
            inds_i = np.mgrid[slices_i].T.reshape(-1,self.lattice.n_dims)
            inds_j = np.mgrid[slices_j].T.reshape(-1,self.lattice.n_dims)

            inds_i = np.ravel_multi_index(inds_i.T, self.lattice.number_of_nodes_per_dimension)
            inds_j = np.ravel_multi_index(inds_j.T, self.lattice.number_of_nodes_per_dimension)
            
            Q_neighbour[inds_i, inds_j] = - dim_scales[d]
            Q_neighbour[inds_j, inds_i] = - dim_scales[d]
        
        #Q = Q_diagonal + Q_neighbour
        
        K = Q_diagonal + Q_neighbour        
        dKdp = dkappa * scipy.sparse.eye(self.n_latent_variables()) * kappa + kappa * scipy.sparse.eye(self.n_latent_variables()) * dkappa
        
        
        if alpha == 1:
            dQdp = dtau * K * tau + \
                    tau * dKdp * tau + \
                    tau * K * dtau  
        
        if alpha == 2:
            dQdp = dtau * K * K * tau + \
                    tau * dKdp * K * tau + \
                    tau * K * dKdp * tau + \
                    tau * K * K *dtau 
        elif alpha > 2:
            raise NotImplementedError('Values of alpha > 2 are not yet supported.')
        
        
        s = np.prod(self.lattice.node_spacing)
        dQdp = s * dQdp
        
        return dQdp

def demo():
    print "Running demo"
    my_lattice = LatticeMesh()
    
    locations = np.array([[0.5, 100.0],
                          [50.5, 500.0],
                          [0.0, 20.0],
                          [10.0, 30.0]])
    
    locations = np.array(my_lattice.axis_coordinates)
    
    N = 10000
    locations = np.random.uniform(0.0, 1.0, (N,2)) * my_lattice.domain_bounds[1,:]+my_lattice.domain_bounds[0,:]
    
    
    my_lattice.locate_nearby_basis_indices( locations )
    #print my_lattice.build_A(locations)

import unittest


class LatticeTest(unittest.TestCase):
    
    """
    
    Tests for evaluation of basis functions/design matrices on a LatticeMesh.
    
    """
    
    def setUp(self):
        """Generate lattices for use in tests."""
        
        dimension_specification = ((0., 1000., 2001),)
        
        self.lattice_1d = LatticeMesh( dimension_specification = dimension_specification,
                                    basis_function = WendlandC4Basis(basis_span = (2.5*1000/2001)) )
        
        dimension_specification = ((0., 10., 11),
                                   (0., 100., 11))
        
        self.lattice_2d = LatticeMesh( dimension_specification = dimension_specification,
                                    basis_function = WendlandC4Basis(basis_span = (2.5*1, 2.5*10)) )
        
        dimension_specification = ((0., 10., 4),
                                   (0., 100., 4),
                                   (0., 1000., 5))
        
        self.lattice_3d = LatticeMesh( dimension_specification = dimension_specification,
                                    basis_function = WendlandC4Basis(basis_span = (2.5*10/4, 2.5*100/4, 2.5*1000/5)) )
        
        
        self.lattices = [self.lattice_1d, self.lattice_2d, self.lattice_3d]
        
        
        dimension_specification = ((0., 100., 101),)
        
        self.lattice_1d_linear = LatticeMesh( dimension_specification = dimension_specification,
                                              basis_function = PiecewiseLinearBasis(basis_span = (1.0,)) )
        
        
        dimension_specification = ((0., 10., 11),
                                   (0., 100., 11))
        
        self.lattice_2d_linear = LatticeMesh( dimension_specification = dimension_specification,
                                    basis_function = WendlandC4Basis(basis_span = (1.0, 10.0)) )
        
        self.linear_lattices = [self.lattice_1d_linear, self.lattice_2d_linear]
        
        
        self.mini_1d = LatticeMesh( dimension_specification = ((0., 7., 8),),
                                    basis_function = PiecewiseLinearBasis(basis_span = (1.0, )) )
        
    def tearDown(self):
        
        self.lattices   = None
        self.lattice_1d = None
        self.lattice_2d = None
        self.lattice_3d = None
        
        self.lattice_1d_linear = None
        self.lattice_2d_linear = None
        self.linear_lattices = None
        
    def atest_A_matrix_at_nodes_1(self):
        """Assert that u == A_matrix(locations) * u when locations are set as the mesh node locations and the basis_span equals the node_spacing"""

        for lattice in self.linear_lattices:
        
            my_locations = lattice.nodes
            
            lattice.basis_function.basis_span = lattice.node_spacing

            for i in range(lattice.n_nodes):

                u = np.zeros((lattice.n_nodes,1))
                u[i,0] = 1.0

                #A = lattice.build_A(my_locations)
#                import matplotlib.pyplot as plt
#                plt.imshow(np.asarray(A.todense()))
#                plt.show()

                
                np.testing.assert_allclose(u, lattice.build_A(my_locations) * u)
            
            #import matplotlib.pyplot as plt
            #plt.imshow(np.asarray(A.todense()))
            #plt.show()
     
    def atest_A_matrix_at_nodes_2(self):
        """Assert that u != A_matrix(locations) * u when locations are set as the mesh node locations and the basis_span is greater than the node_spacing"""
        
        for lattice in self.lattices:
                    
            my_locations = lattice.nodes
            lattice.basis_function.basis_span = lattice.node_spacing * 1.1
            for i in range(lattice.n_nodes):
                u = np.zeros((lattice.n_nodes,1))
                u[i,0] = 1.0
                np.testing.assert_raises(AssertionError, np.testing.assert_allclose, u, lattice.build_A(my_locations) * u)
    
    def atest_A_matrix_basis_functions(self):
        """Assert that each row of A_matrix(locations) are equal to the basis function values when computed using the slow method."""

        for lattice in self.lattices:
            
            my_locations = lattice.nodes
            lattice.basis_function.basis_span = lattice.node_spacing * 25.0
    
            my_A_matrix = lattice.build_A(my_locations)
            
            for j in range(my_A_matrix.shape[0]):
                
                node_distances = my_locations[j,:] - lattice.nodes
    
                basis_values = lattice.basis_function.evaluate(node_distances)
    
                np.testing.assert_allclose(np.asarray(my_A_matrix[j,:].todense()).ravel(), basis_values)
    
    def test_Q_matrix(self):
        
        
        my_spde = LatticeSPDE( lattice = self.lattice_1d )
        #my_Q = my_spde.build_Q(kappa = 0.5, alpha = 1, H = 1.0)
        
        sigma_design_matrix = np.zeros( (my_spde.n_latent_variables(), 2) )
        sigma_design_matrix[:,0] = -np.linspace(0.0, 2.0, my_spde.n_latent_variables() )
        
        rho_design_matrix = np.zeros( (my_spde.n_latent_variables(), 2) )
        rho_design_matrix[:,1] = np.linspace(0.0, 4.0, my_spde.n_latent_variables() )
        
        my_Q = my_spde.build_Q( Q_parameters = np.log(np.array([[np.exp(1.0), np.exp(1.0)]]).T), sigma_design_matrix = sigma_design_matrix, rho_design_matrix = rho_design_matrix, alpha = 2, H = 1.0/my_spde.lattice.node_spacing )
        #my_Q = my_spde.build_Q(alpha = 2)
        #my_Q = my_spde.build_Q_nonstationary(alpha = 1, H = 1 )
        
        #my_Q = my_spde.build_Q()
        #print my_Q
        #print my_spde.lattice.nodes
        import matplotlib.pyplot as plt
        
        #plt.imshow(np.asarray(my_Q.todense()))
        #plt.imshow(np.asarray(scipy.sparse.linalg.spsolve(my_Q, scipy.sparse.eye(my_spde.n_latent_variables())).todense()))
        my_locations = my_spde.lattice.nodes
        my_locations = np.atleast_2d( np.linspace(0.0, 1000.0, 1000.0) ).T
        #print my_locations
        #import advanced_standard.linalg.sparse as sparse
        #u = sparse.sample_mvn(my_Q, 1)
        #process_values = my_spde.build_A(my_locations, Q_matrix = my_Q) * u
        #process_values = my_spde.build_A(my_locations) * u
        #plt.figure()
        #plt.plot(process_values.ravel())
        #plt.show()
    
    def test_Q_mini(self):
        
        my_spde = LatticeSPDE( lattice = self.mini_1d )
        
        print my_spde.build_Q(alpha = 1).todense()
        print my_spde.build_Q(alpha = 2).todense()
    
    def test_dQdtheta(self):
        """Assert that :class:`LatticeSPDE.build_dQdtheta` outputs a matrix that is close to a numerical approximation for small changes in Q parameters."""
        
        my_spde = LatticeSPDE( lattice = self.mini_1d )
        
        test_interval = 1e-8    # small change to apply to parameters for numerical approximation of gradient.
        
        alpha = 2
        
        # zeroth parameter
        d = my_spde.build_dQdtheta(Q_parameters = np.log(np.array([1.0, 1.0])), parameter_index=0, alpha = alpha, sparse_format='csc')
        a =  my_spde.build_Q(Q_parameters = np.log(np.array([1.0, 1.0])), alpha = alpha, sparse_format='csc')
        b =  my_spde.build_Q(Q_parameters = np.log(np.array([1.0, 1.0])) + np.array([test_interval, 0.0]), alpha = alpha, sparse_format='csc')
        
        np.testing.assert_allclose( ((b-a)/test_interval).todense(), d.todense(), atol = 5e-15 )
        
        # first parameter
        d = my_spde.build_dQdtheta(Q_parameters = np.log(np.array([1.0, 1.0])), parameter_index=1, alpha = alpha, sparse_format='csc')
        a =  my_spde.build_Q(Q_parameters = np.log(np.array([1.0, 1.0])), alpha = alpha, sparse_format='csc')
        b =  my_spde.build_Q(Q_parameters = np.log(np.array([1.0, 1.0])) + np.array([0.0, test_interval]), alpha = alpha, sparse_format='csc')
        
        np.testing.assert_allclose( ((b-a)/test_interval).todense(), d.todense(), atol = 5e-15 )
        
    
if __name__ == '__main__':
    """Run test suites"""

    import numpy
    import matplotlib.pyplot as plt
    
    dimension_specification = ((0., 10., 10),
                               (0., 10., 10),)
        
    lattice_2d = LatticeMesh( dimension_specification = dimension_specification,
                                basis_function = WendlandC4Basis(basis_span = 2.5),
                                wrap_dimensions = [True,True])
    
    my_spde = LatticeSPDE( lattice = lattice_2d )
    
    Q = my_spde.build_Q_stationary( np.log(1.0), np.log(3.0), 2.0, 1.0, sparse_format='csc')
    plt.imshow(numpy.asarray(Q.todense()))
    ###
    
    dimension_specification = ((0., 10., 10),)
        
    lattice_1d = LatticeMesh( dimension_specification = dimension_specification,
                                basis_function = WendlandC4Basis(basis_span = 2.5),
                                wrap_dimensions = [True,])

    obs_locations = numpy.array([[-0.5, 0.0, 0.5, 8.5, 9.0, 9.5]]).T

    print lattice_1d.axis_coordinates
    print obs_locations
    obs_map, node_map = lattice_1d.locate_nearby_basis_indices(obs_locations)
    print numpy.array( [obs_map, node_map] ).T
    print numpy.array( [obs_locations[obs_map], lattice_1d.axis_coordinates[node_map]] ).T
    #unittest.main()

    my_spde = LatticeSPDE( lattice = lattice_1d )

    A = my_spde.build_A(obs_locations)
    print A
    Q = my_spde.build_Q_stationary( np.log(1.0), np.log(3.0), 2.0, 1.0, sparse_format='csc')
    
    
    f, axes = plt.subplots(2)
    axes[0].imshow(numpy.asarray(A.todense()))
    axes[1].imshow(numpy.asarray(Q.todense()))
    
    
    
    
    L = scipy.linalg.cholesky(Q.todense(), lower=True)
    u = scipy.linalg.solve(L, numpy.random.normal(0.0, 1.0, (Q.shape[0], 5)))
    
    obs_locations = numpy.linspace(0.0, 10.0, 100)
    values = my_spde.build_A(obs_locations.reshape(-1,1)).dot(u)
    print values
    
    plt.figure()
    plt.plot(obs_locations, values)
    
    plt.show()