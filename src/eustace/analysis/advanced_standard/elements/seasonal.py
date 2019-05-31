"""
Seasonal model element.
"""

import scipy
import numpy
from datetime import datetime
from element import ElementLinear, ElementLinearDesign, ElementPrior
from element import Hyperparameters
from element import SPARSEFORMAT
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE

class SPDEScalarParameters(object):
    """The scalar parameters used to build SPDE."""

    def __init__(self, log_sigma, log_rho):
        """Build with specified parameters."""
        
        self.log_sigma = log_sigma
        self.log_rho = log_rho

    def get_dictionary(self):
        """Get as dictionary of named values."""
        
        return self.__dict__

class SeasonalHyperparameters(Hyperparameters):
    """Hyperparameters for seasonal element spatial harmonics.
       This has 2 parameters for each harmonic."""
    
    def __init__(self, n_spatial_components, common_log_sigma, common_log_rho):
        """Initialise with specified parameters the same for all harmonics."""

        self.harmonics = numpy.tile(numpy.array([ common_log_sigma, common_log_rho ], dtype=numpy.float64), n_spatial_components)
               
    def harmonic_parameters(self, harmonic_index):
        """The SPDE parameters associated with given harmonic index."""
        
        return SPDEScalarParameters(
            log_sigma=self.harmonics[2*harmonic_index],
            log_rho=self.harmonics[2*harmonic_index+1])

    def get_array(self):
        """Retrieve as 1D NumPy array."""

        return self.harmonics

    def set_array(self, values):
        """Set from NumPy array."""

        self.harmonics = values

class SeasonalElement(ElementLinear):
    """Seasonal model element."""

    def __init__(self, n_triangulation_divisions, n_harmonics, include_local_mean):
        """Initialise for given number of harmonics."""

        super(SeasonalElement, self).__init__()
        self.spde = SphereMeshSPDE(level=n_triangulation_divisions, sparse_format=SPARSEFORMAT)
        self.n_harmonics = n_harmonics
        self.include_local_mean = include_local_mean
        self.alpha = 2

    def element_design(self, observationstructure):
        """Build the SeasonalElementDesign instance associated with this element's parameters
           and the given observation structure."""

        return SeasonalElementDesign(observationstructure, self.spde, self.n_harmonics, self.include_local_mean)

    def element_prior(self, hyperparameters):

        return SeasonalElementPrior(hyperparameters, self.spde, self.n_harmonics, self.include_local_mean, self.alpha)


class SeasonalElementDesign(ElementLinearDesign):


    def __init__(self, observationstructure, spde, n_harmonics, include_local_mean):

        self.observationstructure = observationstructure
        self.spde = spde
        self.n_harmonics = n_harmonics
        self.include_local_mean = include_local_mean

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        time_count = 2*self.n_harmonics + int(self.include_local_mean)
        space_count = self.spde.n_latent_variables()
        return time_count * space_count

    def build_A_time(self):
        """
        Build a temporal component for a climatology design matrix.
        
        Builds array of a set of harmonic functions evaluated at the given datetime. The output array is structure such that
        output A_time_out[i,0] corresponds to harmonic function i evaluated at t_datetime.
        
        If self.include_mean is set to True then a constant term is included.
        The number of pairs of harmonic functions to be evaluated are given by self.n_harmonics. 
        
        Args:
        
        * t_datetime (datetime):
            The datetime instance for which to build the design matrix.
        
        Returns:
            2D numpy.ndarray of temporal basis function values.
        
        """
        
        decimal_year = datetime_to_decimal_year(self.observationstructure.time_datetime())
        print "time decimal and full:", decimal_year, self.observationstructure.time_datetime()
        A_list = [ ]
        
        if self.include_local_mean:
            A_list.append(1.0)
        
        for harmonic_order in range(1, self.n_harmonics+1):
            A_list.append(sine_seasonal_harmonic(decimal_year, harmonic_order))
            A_list.append(cosine_seasonal_harmonic(decimal_year, harmonic_order))
        print numpy.array([ A_list ]).shape, numpy.array([ A_list ])
        return numpy.array([ A_list ])

        
    def build_A_space(self):
        '''
        Build a spatial component for a climatology design matrix.
        
        Builds a sparse matrix in the basis functions that are centered on the nodes of self.triangulation are evaluated at locations
        latitudes and longitudes (specified in radians). A_space is structured such that A_space[i,j] takes the value of basis function
        i, at centered at triangulation node i, at the location latitudes[j], longitudes[j] projected onto the piecewise linear surface
        mesh of the triangulated sphere self.triangulation. 
        
        '''
        
        A_space = self.spde.build_A(self.observationstructure.location_polar_coordinates())

        return A_space

    def design_matrix(self):
        """Compute design matrix that maps model parameters to observation space."""

        A_time = self.build_A_time()
        A_space = self.build_A_space()
        return scipy.sparse.kron(A_time, A_space, format=SPARSEFORMAT)

class SeasonalElementPrior(ElementPrior):

    def __init__(self, hyperparameters, spde, n_harmonics, include_local_mean, alpha):

        self.hyperparameters = hyperparameters
        self.spde = spde
        self.n_harmonics = n_harmonics
        self.include_local_mean = include_local_mean
        self.alpha = alpha


    def prior_number_of_state_parameters(self):
        """Return number of state parameters represented by this prior."""

        return self.n_spatial_components() * self.spde.n_latent_variables()

    def n_spatial_components(self):
        """Number of spatial components (each having distinct prior in hyperparameters)."""
        
        return 2*int(self.n_harmonics) + int(self.include_local_mean)


    def harmonic_index_for_spatial_component(self, spatial_component_index):
        """The index of the harmonic corresponding to specified spatial component block."""
        
        return ( int(spatial_component_index) + int(self.include_local_mean) ) / 2


    def prior_precision(self):
        """Build the prior precision matrix for the seasonal model.
        
        This precision matrix is built as a block diagonal matrix of 
        self.n_components() blocks.  Settings are shared for each harmonic order
        such that the first block corresponds to a zeroth ordered harmonic, if one exists
        in the model, folllowed by pairs of blocks that share parameters for sine-cosine
        pairs.
        
        
        Args:
        
            Q_parameters (list):
                * A concatenated list of hyperparameter values for each harmonic order.
        
            Q_settings (list):
                * A list of pairs of self.n_harmonics() dictionaries of kwargs to pass 
                for each block.
        
        """
        
        Q_list = []
        
        for spatial_component_index in range(self.n_spatial_components()):
            
            harmonic_index = self.harmonic_index_for_spatial_component(spatial_component_index)
            spde_parameters = self.hyperparameters.harmonic_parameters(harmonic_index)
            block_Q = self.spde.build_Q_stationary(alpha=self.alpha, **spde_parameters.get_dictionary())
            Q_list.append(block_Q)
          
        return scipy.sparse.block_diag(Q_list, format=SPARSEFORMAT)

    def prior_precision_derivative(self, parameter_index):
        """Compute prior covariance derivative with respect to hyperparameters."""
        
        dQdp_list = []
        
        for spatial_component_index in range(self.n_spatial_components()):
            
            harmonic_index = self.harmonic_index_for_spatial_component(spatial_component_index)
            
            if (parameter_index / 2) == harmonic_index:
                                
                spde_parameters = self.hyperparameters.harmonic_parameters(harmonic_index)                
                component_parameter_index = parameter_index % 2
                dQdp_block = self.spde.build_dQdp_stationary(alpha=self.alpha, parameter_index=component_parameter_index, **spde_parameters.get_dictionary())
        
            else:
                subQdim = self.spde.n_latent_variables()
                dQdp_block = scipy.sparse.coo_matrix( (subQdim, subQdim) )  # create empty matrix for this sub-block
          
            dQdp_list.append(dQdp_block)
    
        dQdp_matrix = scipy.sparse.block_diag(dQdp_list, format=SPARSEFORMAT)
        
        return dQdp_matrix

def cosine_seasonal_harmonic(decimal_year, order):
    
    return numpy.cos( order * 2.0 * numpy.pi * decimal_year )

def sine_seasonal_harmonic(decimal_year, order):
    
    return numpy.sin( order * 2.0 * numpy.pi * decimal_year )
    
def datetime_to_decimal_year(t_datetime):
    """Convert a datetime to decimal year."""

    # The datetime object at start of year
    year_start = datetime(t_datetime.year, 1, 1)

    # The datetime object at start of next year
    year_end = datetime(t_datetime.year + 1, 1, 1)

    # Seconds since start of year
    t_seconds = numpy.float64( (t_datetime - year_start).total_seconds() )

    # Total seconds in this year
    total_seconds = numpy.float64( (year_end - year_start).total_seconds() )

    # Compute decimal
    return numpy.float64(t_datetime.year) + (t_seconds / total_seconds)
