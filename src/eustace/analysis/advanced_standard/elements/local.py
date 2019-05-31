"""Local model."""

import numpy as np
import scipy
from element import SPARSEFORMAT
from element import ElementLinear, ElementPrior, ElementLinearDesign
from element import Hyperparameters
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE

class LocalHyperparameters(Hyperparameters):
    """Hyperparameters for local element."""
    
    def __init__(self, log_sigma, log_rho):
        """Initialise with specified parameters."""

        self.log_sigma = log_sigma
        self.log_rho = log_rho
       
    def get_array(self):
        """Retrieve as 1D NumPy array."""
        
        valueslist = [ self.log_sigma, self.log_rho ]
        return np.array(valueslist)
        
    def set_array(self, values):
        """Set from NumPy array."""

        self.log_sigma = values[0]
        self.log_rho = values[1]
                
    def get_dictionary(self):
        """Get as dictionary of named values."""
        
        return self.__dict__
        
class LocalElement(ElementLinear):    
    """
    Class for evaluating design and precision matrices for SPDE models on a triangulated sphere.    
    """
    
    def __init__(self, n_triangulation_divisions):
        """
        Setup triangulation of sphere for :class:`LocalModel`.

        Args:
        
        * n_triangulation_divisions (int):
            Number of subdivisions of icosohedron for triangulation of sphere.
        
        """
                    
        super(LocalElement, self).__init__()
        self.spde = SphereMeshSPDE(level=n_triangulation_divisions, sparse_format=SPARSEFORMAT)

    def element_design(self, observationstructure):

        return LocalDesign(observationstructure, self.spde)

    def element_prior(self, hyperparameters):
        """
        Return LocalPrior instance for local model, using given hyperparameters.
        """
    
        return LocalPrior(hyperparameters, self.spde)


class LocalDesign(ElementLinearDesign):

    def __init__(self, observationstructure, spde):

        self.observationstructure = observationstructure
        self.spde = spde

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        return self.spde.n_latent_variables()
    
    def design_matrix(self):
        """
        Build the spatial design matrix for the local model.
        
        Builds a sparse matrix in the basis functions that are centered on the nodes of self.triangulation are evaluated at locations
        latitudes and longitudes (specified in radians). A_space is structured such that A_space[i,j] takes the value of basis function
        i, at centered at triangulation node i, at the location latitudes[j], longitudes[j] projected onto the piecewise linear surface
        mesh of the triangulated sphere self.triangulation.
        
            
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
            
        """        
    
        return self.spde.build_A(self.observationstructure.location_polar_coordinates())

class LocalPrior(ElementPrior):

    def __init__(self, hyperparameters, spde):

        self.hyperparameters = hyperparameters
        self.spde = spde
        self.alpha = 2
    
    def prior_number_of_state_parameters(self):
        """Return number of state parameters represented by this prior."""

        return self.spde.n_latent_variables()

    def prior_precision(self):
        """
        Builds a precision matrix for the local model.
        
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
               
        """
        
        return self.spde.build_Q_stationary(alpha=self.alpha, **self.hyperparameters.get_dictionary())
    
    def prior_precision_derivative(self, parameter_index):
        """
        Builds the derivative of the precision matrix for the local model with respect to the precision matrix hyperparmaters.
        
        Kwargs:
        
        * parameter_index:
            An index to the precision matrix parameter that the derivative is to be computed with 
            respect to.
                    
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
                
        """
        
        return self.spde.build_dQdp_stationary(alpha=self.alpha, parameter_index=parameter_index, **self.hyperparameters.get_dictionary())
