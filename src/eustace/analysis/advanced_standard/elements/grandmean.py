"""Include a state parameter which is an offset to all parameters."""

from element import ElementLinear
from element import ElementLinearDesign
from element import SPARSEFORMAT
from covariate import CovariatePrior
import numpy
import scipy

class GrandMeanElement(ElementLinear):

    def __init__(self):
        """Constuct."""

        super(GrandMeanElement, self).__init__()

    def element_design(self, observationstructure):

        number_of_observations = observationstructure.location_polar_coordinates().shape[0]
        return GrandMeanDesign(number_of_observations)

    def element_prior(self, hyperparameters):
        """Return a CovariatePrior object with one state parameter corresponding to the grand mean."""

        return CovariatePrior(hyperparameters, number_of_state_parameters=1)

class GrandMeanDesign(ElementLinearDesign):
    """Compute design matrix row for grand mean of observations."""

    def __init__(self, number_of_observations):
        """Initialise for specified number of observations."""

        self.number_of_observations = number_of_observations

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector (always just one for grand mean)."""

        return 1

    def design_matrix(self):
        """Get a row of ones for all observations."""

        return scipy.sparse.csc_matrix(numpy.ones((self.number_of_observations,1))).asformat(SPARSEFORMAT)

