"""Elements to model observational bias."""

import numpy as np
import scipy
from element import ElementLinear
from element import ElementLinearDesign
from element import SPARSEFORMAT
from local import LocalPrior
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE

class SpatialBiasElement(ElementLinear):    
    """
    Class for evaluating design and precision matrices for SPDE models on a triangulated sphere.    
    """

    def __init__(self, groupname, n_triangulation_divisions):
        """
        Setup triangulation of sphere for LocalBias model.

        Args:
        
        * groupname:
            The name of the bias group that this bias corresponds to. Not currently used.
        
        * n_triangulation_divisions (int):
            Number of subdivisions of icosohedron for triangulation of sphere.
        
        """
                    
        super(SpatialBiasElement, self).__init__()
        
        self.groupname = groupname
        
        self.spde = SphereMeshSPDE(level=n_triangulation_divisions, sparse_format=SPARSEFORMAT)
        self.number_of_biases = self.spde.n_latent_variables()

    def element_design(self, observationstructure):

        return SpatialBiasDesign(observationstructure, self.spde, self.groupname)

    def element_prior(self, hyperparameters):
        """
        Return LocalPrior instance for local model, using given hyperparameters.
        """
    
        return LocalPrior(hyperparameters, self.spde)

class SpatialBiasDesign(ElementLinearDesign):

    def __init__(self, observationstructure, spde, groupname):
        
        self.observationstructure = observationstructure
        self.spde = spde
        self.groupname = groupname
        
        # Count total observations (including those not affected by this bias)
        self.number_of_observations = observationstructure.number_of_observations()

        # Look for this effect on observation structure.
        # The bias will map onto the observations if the effect is not None.
        self.effect = self.observationstructure.covariate_effect(self.groupname)

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        return self.spde.n_latent_variables()
    
    def design_matrix(self):
        """
        Build the spatial design matrix for the local bias model.
            
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format in SPARSEFORMAT.
            
        """        
    
        if self.effect is not None:
            # We do not actually use the effect information as all observations with a non-None effect are expected
            # to be affected based on their spatial location. In future the effect information it could be used to 
            # indicate the strength of an interaction, e.g. for locally correlated uncertainties of varying magnitude.
            
            A = self.spde.build_A(self.observationstructure.location_polar_coordinates())

        else:
            A = scipy.sparse.dok_matrix((self.number_of_observations, self.design_number_of_state_parameters()))

        return A.asformat(SPARSEFORMAT)

