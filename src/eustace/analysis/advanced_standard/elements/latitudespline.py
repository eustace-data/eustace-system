"""

Latitude spline SPDE element.

"""

import scipy
import numpy
from element import ElementLinear, ElementLinearDesign, ElementPrior
from element import SPARSEFORMAT
from eustace.analysis.advanced_standard.stats.spde.lattice import LatticeSPDE, WendlandC4Basis

HALF_SPHERE_DEGREES = 180.0
PADDING_NODES = 180

class LatitudeSplineElement(ElementLinear):        
    
    def __init__(self, alpha, n_nodes, overlap_factor, H):
        """Initialise SPDE definitions.
        
        1D SPDE model in latitude with n_nodes spanning that range [-90, 90] degrees.
        
        """
        
        super(LatitudeSplineElement, self).__init__()
        
        node_separation = HALF_SPHERE_DEGREES / (n_nodes-1)
        n_nodes_padded = n_nodes + PADDING_NODES * 2
        padding_distance = PADDING_NODES * node_separation
        
        self.alpha = alpha
        self.spde = LatticeSPDE.construct(
            dimension_specification = [(-90.0-padding_distance, 90.0+padding_distance, n_nodes_padded)],
            basis_function = WendlandC4Basis(),
            overlap_factor = overlap_factor)
        self.H = H
    
    def element_design(self, observationstructure):
        """Build the SeasonalElementDesign instance associated with this element's parameters
           and the given observation structure."""

        return LatitudeSplineElementDesign(observationstructure, self.spde, self.alpha, self.H)

    def element_prior(self, hyperparameters):
        """Build prior instance associated with given hyperparameters."""

        return LatitudeSplineElementPrior(hyperparameters, self.spde, self.alpha, self.H)


class LatitudeSplineElementDesign(ElementLinearDesign):
    """Base class for SPDE design. Inherit from this together with ElementLinearDesign or ElementNonlinearDesign as required."""

    def __init__(self, observationstructure, spde, alpha, H):
        """Initialise space and time SPDE designs."""

        super(LatitudeSplineElementDesign, self).__init__()

        self.spde = spde

        self.alpha = alpha
        self.H = H
        
        self.latitudes = observationstructure.location_polar_coordinates()[:,0]

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        return self.spde.n_latent_variables()

    def design_matrix(self):
        """
        Build the design matrix for the latitude model.
            
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
            
        """        
    
        return self.spde.build_A(self.latitudes)

    
class LatitudeSplineElementPrior(ElementPrior):

    def __init__(self, hyperparameters, spde, alpha, H):
        
        super(LatitudeSplineElementPrior, self).__init__()
        
        self.hyperparameters = hyperparameters
        self.spde = spde
        self.alpha = alpha
        self.H = H

    def prior_number_of_state_parameters(self):

        return self.spde.n_latent_variables()
    
    def prior_precision(self):
        """
        Builds a precision matrix for the SPDE model.
        
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
               
        """
        
        return self.spde.build_Q_stationary(log_sigma = self.hyperparameters.log_sigma, log_rho = self.hyperparameters.log_rho, alpha=self.alpha, H = self.H, sparse_format = SPARSEFORMAT)

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
        
        return self.spde.build_dQdp_stationary(log_sigma = self.hyperparameters.log_sigma, log_rho = self.hyperparameters.log_rho, alpha=self.alpha, H = self.H, parameter_index=parameter_index, sparse_format = SPARSEFORMAT)