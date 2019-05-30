"""Elements to model insitu-land observational bias."""

import numpy as np
import scipy

from bias import BiasElement
from bias import BiasDesign
from element import ElementLinearDesign
from element import SPARSEFORMAT
from eustace.preprocess.fileio.insitu_land_breakpoints import ObservationBreakPointSourceInsituLand

class InsituLandBiasElement(BiasElement):
    """
    
    Represent covariate-like insitu-land observational biases as a linear model.

    """

    GROUPNAME = 'insitu_land'
    OBSERVABLE = 'merged_break'

    def __init__(self, filename, apply_policy = False, cut_value = 0):
        """
        
        A bias group name names a group of covariates that share prior hyperparameters,
        e.g. relating to a given sensor type.

        The indices of any observations affected can be retrieved using the
        covariate_effect method of the observation structure.

	For insitu-land sources the number of covariates is determined from the number of stations breakpoints, 
        """
        
	breakpoints_reader = ObservationBreakPointSourceInsituLand(filename)
	self.observed_breakpoints = breakpoints_reader.observations(InsituLandBiasElement.OBSERVABLE)
	breakpoints_reader.dataset.close()
	if apply_policy:
	  self.observed_breakpoints.apply_policy('HARD_CUTOFF', threshold = cut_value)
	 
	print("Number of BIAS latent variables = ", self.observed_breakpoints.break_time.size)
        super(InsituLandBiasElement, self).__init__(InsituLandBiasElement.GROUPNAME, self.observed_breakpoints.break_time.size)
        

    def element_design(self, observationstructure):	
        """Build the BiasDesign instance associated with this element's parameters
           and the given observation structure."""

        return InsituLandBiasDesign(observationstructure, self.groupname, self.number_of_biases, self.observed_breakpoints)

class InsituLandBiasDesign(ElementLinearDesign):
    """

    Build design matrices for given covariate group and observation structure.

    """

    def __init__(self, observationstructure, groupname, number_of_biases, observed_breakpoints):

        self.groupname = groupname
        self.number_of_biases = number_of_biases
        self.observed_breakpoints = observed_breakpoints

        # Count total observations (including those not affected by this bias)
        self.number_of_observations = observationstructure.number_of_observations()

        # Look for this effect on observation structure        
        self.effect = observationstructure.covariate_effect(self.groupname, breakpoints = observed_breakpoints)


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