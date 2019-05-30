"""Combination of elements."""

from element import Element
from element import ElementDesign
from element import ElementPrior
from element import Hyperparameters
from element import SPARSEFORMAT
import numpy
import scipy

class CombinationHyperparameters(Hyperparameters):
    """Hyperparameter object built using hyperparameters from
       multiple elements."""
    
    def __init__(self, elementparameters):
        
        super(CombinationHyperparameters, self).__init__()
        self.elementparameters = elementparameters
            
    def get_array(self):
        """Retrieve as 1D NumPy array."""
        
        vectors = [ parameters.get_array() for parameters in self.elementparameters ]
        return numpy.hstack(vectors)
        
    def set_array(self, values):
        """Set from NumPy array."""
        
        # List with one range for each element        
        elementranges = self.get_element_ranges()
        
        # Total count is one greater than last index of last element
        totalcount = elementranges[-1][-1] + 1

        # Check total count matches number of values being set
        if values.shape != (totalcount,):
            raise ValueError('wrong vector shape: expected {0} got {1}', (totalcount,), values.shape)

        # Set values        
        for index, parameters in enumerate(self.elementparameters):            
            parameters.set_array(values[elementranges[index]])

    def get_element_ranges(self):
        """Retrieve list containing ranges for each element within vector of all parameters."""
        
        ranges = [ ]
        startindex = 0
        for parameters in self.elementparameters:            
            count = parameters.get_array().shape[0]
            nextstartindex = startindex + count
            ranges.append(range(startindex, nextstartindex))
            startindex = nextstartindex
        return ranges
       

class CombinationElement(Element):
    """A combination element comprises several individual analysis elements."""
           
    def __init__(self, combination):
        """Initialise form list of component elements."""
        
        super(CombinationElement, self).__init__()
        self.combination = combination

    def isnonlinear(self):
        """Will return True if any individual element represents a nonlinear system, False otherwise."""

        return any([ element.isnonlinear() for element in self.combination ])

    def element_design(self, observationstructure):

        designlist = [ element.element_design(observationstructure) for element in self.combination ]
        return CombinationDesign(designlist)

    def element_prior(self, hyperparameters):
        """
        Return a combination prior.
        """

        priorlist = [ element.element_prior(hyperparameters.elementparameters[index]) for index, element in enumerate(self.combination) ]
        return CombinationPrior(hyperparameters, priorlist)


class CombinationDesign(ElementDesign):
    """Build design matrices for a combination of elements."""

        
    def __init__(self, designlist):

        self.designlist = designlist


    def isnonlinear(self):
        """Returns True if any element of the combination is nonlinear, false otherwise."""

        return any([ design.isnonlinear() for design in self.designlist ])


    def design_number_of_state_parameters(self):
        """Number of state parameters.  For combination element this is the sum of state parameters from all sub-elements."""

        return sum([ element_design.design_number_of_state_parameters() for element_design in self.designlist ])


    def design_initial_state(self):
        """Initial state. Uses initial substate information from nonlinear subdesigns where applicable."""

        initialstate = None
        
        if self.isnonlinear():

            # Initialise to zero
            initialstate = numpy.zeros((self.design_number_of_state_parameters(),))

            # Populate substates from element initial state values
            state_startindex = 0
            for subdesign in self.designlist:
                state_count = subdesign.design_number_of_state_parameters()
                state_nextstartindex = state_startindex + state_count
                if subdesign.isnonlinear():
                    initialstate[state_startindex:state_nextstartindex] = subdesign.design_initial_state()
                state_startindex = state_nextstartindex
        
        return initialstate


    def element_states(self, currentstate):
        """List of vectors containing the subsets of state vector that correspond to each element."""

        if currentstate is None:

            # If no current state specified this implies it's a fully linear system 
            element_state_list = [ None for element in self.designlist ]

        else:

            # We have a state specified - split it into element states
            element_state_list = [ ]
            state_startindex = 0
            for element in self.designlist:

                state_count = element.design_number_of_state_parameters()
                state_nextstartindex = state_startindex + state_count
                element_state = currentstate[state_startindex:state_nextstartindex]
                element_state_list.append(element_state)
                state_startindex = state_nextstartindex

        return element_state_list

    def design_function(self, currentstate, rectified=False):
        """
        Compute function f that maps model parameters to observation space,
        at given model state.
        """

        states = self.element_states(currentstate)
        vectorlist = [ element.design_function(states[elementindex], rectified) for elementindex, element in enumerate(self.designlist)  ]
        result = numpy.zeros(vectorlist[0].shape)
        for vector in vectorlist:
            result += vector
        return result


    def design_jacobian(self, currentstate):
        """
        Compute jacobian matrix J = df / dx of the function f that maps model parameters to 
        observation space, evaluated at given model state.
        y = f(x) + error terms
        """

        states = self.element_states(currentstate)
        matrixlist = [ element.design_jacobian(states[elementindex]) for elementindex, element in enumerate(self.designlist)  ]
        return scipy.sparse.hstack(matrixlist, format=SPARSEFORMAT) 

    def design_matrix(self):
        """Return the joint design matrix for the combination of elements"""
        
        #return scipy.sparse.hstack( [design.design_matrix() for design in self.designlist], format=SPARSEFORMAT)
        design_matrix_list = []
        for design in self.designlist:
            print design
            design_matrix_list.append( design.design_matrix() )
            
        return scipy.sparse.hstack(design_matrix_list , format=SPARSEFORMAT)
        

class CombinationPrior(ElementPrior):
    """Build prior precision matrices for a combination of elements."""

    def __init__(self, hyperparameters, priorlist):
        """Initialise for given design size and hyperparameters."""

        self.hyperparameters = hyperparameters
        self.priorlist = priorlist

    def prior_number_of_state_parameters(self):
        """Return number of state parameters represented by this prior."""

        return sum([ prior.prior_number_of_state_parameters() for prior in self.priorlist ])

    def prior_precision(self):
        """Compute prior precision for given hyperparameters."""

        # List of priors
        precision_list = [ elementprior.prior_precision() for elementprior in self.priorlist ]

        # Form block diagonal result
        return scipy.sparse.block_diag(precision_list, format=SPARSEFORMAT)

    def prior_precision_derivative(self, parameter_index):
        """Compute prior precision partial derivatives at given hyperparameters."""

        # Hyperparameter indices corresponding to each element
        elementranges = self.hyperparameters.get_element_ranges()
        
        # List to buid
        derivative_list = [ ]        
        
        # Zeros for everything except the specified range
        for elementindex, elementrange in enumerate(elementranges):
        
            # The requested parameter will only apply to one element
            if parameter_index in elementrange:

                # Hyperparameters for this element
                elementparameters = self.hyperparameters.elementparameters[elementindex]

                # The index relative to this element
                relativeindex = parameter_index - elementrange[0]
                
                # Get derivative
                derivative = self.priorlist[elementindex].prior_precision_derivative(relativeindex)

            else:
                
                # Zeros for elements not affected by the parameter
                dimension = self.priorlist[elementindex].prior_number_of_state_parameters()
                derivative = scipy.sparse.csc_matrix((dimension,dimension))
                
            derivative_list.append(derivative)


        # Block diagonal result
        return scipy.sparse.block_diag(derivative_list, format=SPARSEFORMAT)


    def element_states(self, currentstate):
        """List of vectors containing the subsets of state vector that correspond to each element."""

        if currentstate is None:

            # If no current state specified this implies it's a fully linear system 
            element_state_list = [ None for prior in self.priorlist ]

        else:

            # We have a state specified - split it into element states
            element_state_list = [ ]
            state_startindex = 0
            for prior in self.priorlist:

                state_count = prior.prior_number_of_state_parameters()
                state_nextstartindex = state_startindex + state_count
                element_state = currentstate[state_startindex:state_nextstartindex]
                element_state_list.append(element_state)
                state_startindex = state_nextstartindex

        return element_state_list