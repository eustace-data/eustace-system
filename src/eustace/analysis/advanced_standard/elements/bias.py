"""Elements to model observational bias."""

import numpy as np
import scipy
from element import ElementLinear
from element import ElementLinearDesign
from element import SPARSEFORMAT
from covariate import CovariatePrior

class BiasElement(ElementLinear):
    """
    
    Represent covariate-like observational biases as a linear model.
    
    """
    
    def __init__(self, groupname, number_of_biases):
        """
        
        A bias group name names a group of covariates that share prior hyperparameters,
        e.g. relating to a given sensor type.

        The indices of any observations affected can be retrieved using the
        covariate_effect method of the observation structure.

        """
        
        super(BiasElement, self).__init__()
        self.groupname = groupname
        self.number_of_biases = number_of_biases

    def element_design(self, observationstructure):
        """Build the BiasDesign instance associated with this element's parameters
           and the given observation structure."""

        return BiasDesign(observationstructure, self.groupname, self.number_of_biases)

    def element_prior(self, hyperparameters):
        """
        Build CovariatePrior instance which produces diagnonal precision matrices with a
        constant value for the biases affected by this design.
        """
        
        return CovariatePrior(hyperparameters, number_of_state_parameters=self.number_of_biases)

    
class BiasDesign(ElementLinearDesign):
    """

    Build design matrices for given covariate group and observation structure.

    """

    def __init__(self, observationstructure, groupname, number_of_biases):

        self.groupname = groupname
        self.number_of_biases = number_of_biases

        # Count total observations (including those not affected by this bias)
        self.number_of_observations = observationstructure.number_of_observations()

        # Look for this effect on observation structure        
        self.effect = observationstructure.covariate_effect(self.groupname)
        print self.effect


    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        return self.number_of_biases
                        
    def design_matrix(self):
        """Construct a design matrix for the covariate model.
        
        Build design matrix for a set of n_observations for the specified
        covariate group name.
            
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
        
        """

        A = scipy.sparse.dok_matrix((self.number_of_observations, self.number_of_biases))

        if self.effect is not None:

            # Set design elements to 1.0
            A[self.effect[:,0], self.effect[:,1]] = 1.0

        return A.asformat(SPARSEFORMAT)
