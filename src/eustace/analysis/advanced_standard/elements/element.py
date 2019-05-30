"""
Define an analysis observation model.
"""

import numpy

SPARSEFORMAT = 'csr'
"""Format used for sparse matrices."""

class Element(object):
    """
    An analysis element defines constants required for computation of linear model design matrix,
    and coresponding prior precision matrix. It can build an ElementDesign instance which
    can calculate matrices used to solve the system, and can produce instances needed to
    build prior precision matrices from hyperparameters.
    Abstract base class.
    """

    def __init__(self):
        """Initialise."""

        pass

    def isnonlinear(self):
        """Determine whether the element design will represent an exact linear system
           or linear approximations to a nonlinear system."""

        raise NotImplementedError

    def element_design(self, observationstructure):
        """Build the ElementDesign instance associated with this element's parameters
           and the given observation structure."""

        raise NotImplementedError

    def element_prior(self, hyperparameters):
        """
        Build ElementPrior instance which can produce prior precision 
        matrices corresponding to this element for the given set of
        hyperparameters.
        """

        raise NotImplementedError

class ElementLinear(object):
    """Abstract base class for linear elements."""

    def isnonlinear(self):
        """Always returns False for linear systems."""

        return False

class ElementNonlinear(object):
    """Abstract base class for non-linear elements."""

    def isnonlinear(self):
        """Always returns True for non-linear systems."""

        return True

class ElementDesign(object):
    """
    This brings together element parameters with a given observation structure
    to compute design matrices.  It can produce an ElementPrior instance
    which can calculate prior precision matrices for a given set of hyperparameters.
    """

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        raise NotImplementedError
        

    def design_function(self, currentstate):
        """
        Compute function f that maps model parameters to observation space,
        at given model state.
        """

        raise NotImplementedError

        
    def design_jacobian(self, currentstate):
        """
        Compute jacobian matrix J = df / dx of the function f that maps model parameters to 
        observation space, evaluated at given model state.
        y = f(x) + error terms
        """
    
        raise NotImplementedError


    def isnonlinear(self):
        """Determine whether this is a linear or nonlinear design element.
           Linear designs should implement the design_matrix method.
           Nonlinear designs should implement the design_jacobian method,
           design_function, and design_initial state methods."""

        raise NotImplementedError


    def design_initial_state(self):
        """Initial state of nonlinear optimisation, or None for linear systems."""

        raise NotImplementedError

       
class ElementLinearDesign(ElementDesign):
    """
    Child of ElementDesign to implement a linear system.
    """

    def isnonlinear(self):
        """Returns False in this class, as this represents a linear system."""

        return False

    def design_initial_state(self):
        """Returns None as there is no initial state for linear systems (they are solved directly)."""

        return None

    def design_matrix(self):
        """
        Compute design matrix A that maps model parameters to observation space
        y = A.x + error terms
        """
    
        raise NotImplementedError

    def design_jacobian(self, currentstate):
        """For linear systems this is equivalent to design_matrix method and the currentstate parameter is ignored."""

        return self.design_matrix()

    def design_function(self, currentstate, rectified=False):
        """Compute y = A.x."""

        if rectified:
            return numpy.abs(self.design_matrix()).dot(currentstate)
        else:
            return self.design_matrix().dot(currentstate)


class ElementNonlinearDesign(ElementDesign):
    """
    Child of ElementDesign to implement a nonlinear system.
    """

    def isnonlinear(self):
        """Returns True in this class, as this represents a nonlinear system."""

        return True

class ElementPrior(object):

    def prior_number_of_state_parameters(self):
        """Return number of state parameters represented by this prior."""

        raise NotImplementedError

    def prior_precision(self):
        """Compute prior precision for given hyperparameters.
           This might also depend on the observation structure,
           for example if one bias parameter is used per observation."""

        raise NotImplementedError

    def prior_precision_derivative(self, parameter_index):
        """Compute prior precision partial derivatives at given hyperparameters."""

        raise NotImplementedError
               
class ObservationStructure(object):
    """
    An observation structure describes observations at a specific time.
    Abstract base class.
    """
    
    def time_index(self):
        """
        A discretised time index at appropriate resolution for input.
        In EUSTACE fullstace system these correspond to day numbers since
        01/01/1850.
        """
        
        raise NotImplementedError

    def time_datetime(self):
        """The absolute datetime of these observations
        (needed for seasonal model and factor analysis which 
        must compute a fractional year value).
        """

        raise NotImplementedError

    def number_of_observations(self):
        """Return total observations at this time index."""

        raise NotImplementedError
        
    def location_polar_coordinates(self):
        """Array of polar coordinates of observation location on unit sphere."""

        raise NotImplementedError        
                
    def covariate_effect(self, groupname, **kwargs):
        """Retrieve details of bias that affects these observations.
           NumPy matrix of integers with 2 columns.  
           First column gives indices of observations affected.
           Second column is index of corresponding bias parameter affecting it.
           None if no effect defined."""
        
        return None

    def observation_vector(self):
        """Array of observations."""

        raise NotImplementedError

    def observation_precision(self):
        """Observation precision matrix (sparse)."""

        raise NotImplementedError
       
class Hyperparameters(object):
    """
    Hyperparameters are parameters used to compute prior precision matrices,
    which may be optimised to best explain the data.
    Abstract base class.
    """
        
    def get_array(self):
        """Retrieve as 1D NumPy array."""
        
        raise NotImplementedError
        
        
    def set_array(self, values):
        """Set from NumPy array."""
        
        raise NotImplementedError
