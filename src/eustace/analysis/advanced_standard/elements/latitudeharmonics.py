"""
Latitude harmonics element.
"""

import scipy
import numpy
from element import ElementLinear, ElementLinearDesign, ElementPrior
from element import Hyperparameters
from element import SPARSEFORMAT
from covariate import CovariatePrior
from combination import CombinationPrior

class LatitudeFunction(object):
    """Wrapper for any function over latitude."""

    def __init__(self, numpymethod, order):
        """Construct to run specified NumPy method after multiplying latitude by order."""

        self.numpymethod = numpymethod
        self.order = order 

    def compute(self, latitudes):
        """Compute result"""

        resultvector = self.numpymethod(numpy.radians(latitudes) * self.order)
        return resultvector.reshape((latitudes.shape[0], 1))

class LatitudeHarmonicsElement(ElementLinear):
    """Latitude harmonics model element."""

    """Fixed list of all harmonic functions to use."""
    HARMONICS = [ 
        LatitudeFunction(numpy.cos, 2.0),
        LatitudeFunction(numpy.sin, 2.0),
        LatitudeFunction(numpy.cos, 4.0),
        LatitudeFunction(numpy.sin, 4.0),
    ]

    def __init__(self):
        """Initialise for given number of harmonics (including zeroth)"""

        super(LatitudeHarmonicsElement, self).__init__()
        
    def element_design(self, observationstructure):
        """Build the SeasonalElementDesign instance associated with this element's parameters
           and the given observation structure."""

        return LatitudeHarmonicsElementDesign(observationstructure)

    def element_prior(self, hyperparameters):
        """Build prior instance associated with given hyperparameters."""

        return LatitudeHarmonicsPrior(hyperparameters)

class LatitudeHarmonicsElementDesign(ElementLinearDesign):

    def __init__(self, observationstructure):
        """Construct using latitudes from given observation structure"""

        self.latitudes = observationstructure.location_polar_coordinates()[:,0]

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        return len(LatitudeHarmonicsElement.HARMONICS)

    def design_matrix(self):
        """Compute design matrix that maps model parameters to observation space."""

        matrixlist = [ scipy.sparse.coo_matrix(harmonic.compute(self.latitudes)) for harmonic in LatitudeHarmonicsElement.HARMONICS ]
        return scipy.sparse.hstack(matrixlist, format=SPARSEFORMAT)

class LatitudeHarmonicsPrior(CombinationPrior):
    """Subclass of CombinationPrior to give independent prior for each harmonic."""

    def __init__(self, hyperparameters):

        # Compute expected number
        expected_dimension = len(LatitudeHarmonicsElement.HARMONICS)

        # Get parameters themselves (assumes CombinationHyperparameters)
        hlist = hyperparameters.elementparameters

        # Check number of hyperparameters is as expected
        if len(hlist) != expected_dimension:
            raise ValueError('Expected {0} hyperparameters but got {1}'.format(
                expected_dimension,
                len(hlist)))

        # Construct combination
        super(LatitudeHarmonicsPrior, self).__init__(
            hyperparameters, 
            [ CovariatePrior(h, 1) for h in hlist ])
