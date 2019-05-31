"""Define calculation of priors for covariates."""

import importlib
import json
import numpy as np
import scipy

from combination import CombinationHyperparameters
from element import ElementPrior
from element import Hyperparameters
from element import SPARSEFORMAT


class CovariateHyperparameters(Hyperparameters):
    """By definition, parameters affected by a covariate element share a single hyperparameter."""
        
    def __init__(self, value):
        """Initialise with given single value."""
        
        self.value = value
                
    def get_array(self):
        """Retrieve as 1D NumPy array."""
        
        return np.array([self.value])
                
    def set_array(self, values):
        """Set from NumPy array."""
        
        self.value = values[0]
    
    
class CovariatePrior(ElementPrior):
    
    def __init__(self, hyperparameters, number_of_state_parameters):

        self.hyperparameters = hyperparameters
        self.number_of_state_parameters = number_of_state_parameters

    def prior_number_of_state_parameters(self):
        """Return number of state parameters represented by this prior."""

        return self.number_of_state_parameters

    def prior_precision(self):
        """Return prior precision matrix.
        
        Allows computation of precision matrices in which specified blocks of the
        precision matrix are scaled in amplitude by 1 / exp(Q_parameters)**2. With this
        parameterisation, exp(Q_parameters)**2 effectively controls the prior variance
        for specified blocks of the precision matrix.
                
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.                
        """
  
        prior = np.exp(-2.0*self.hyperparameters.value)
        return scipy.sparse.diags(np.repeat(prior, self.number_of_state_parameters), format=SPARSEFORMAT)
    
    def prior_precision_derivative(self, parameter_index):
        """Compute derivative of Q matrix by Q_parameters[parameter_index]
                                
            * parameter_index:

                Expected to be zero as there is only one hyperparameter.
        
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.                
        """
        
        if parameter_index == 0:

            derivative = -2.0 * np.exp(-2.0*self.hyperparameters.value)
            return scipy.sparse.diags(np.repeat(derivative, self.number_of_state_parameters), format=SPARSEFORMAT)

        else:

            raise ValueError('parameter index out of range')

class LoadCovariateElement(object):
      """Helper class for building a covariate element according to user instructions. It is useful for adding fields to a given model component"""

      AVAILABLE_COVARIATES = ['latitude_harmonics', 'altitude', 'coastal_influence']

      def __init__(self, json_descriptor):
	  with open(json_descriptor) as data_file:    
	      self.data = json.load(data_file)
	      data_file.close()

      def check_keys(self):
	  """Checks that the user given covariates match the list of available ones"""
	  for key in self.data.keys():
	      if key not in LoadCovariateElement.AVAILABLE_COVARIATES:
		  error_message = 'Covariate {0} currently not implemented in the system, please remove it from the json descriptor.'.format(key)
		  raise ValueError(error_message)
      
      def load_covariates_and_hyperparameters(self):
	  """Load list covariates elements and hyperparameters according to user given instrunctions"""
          
	  list_of_covariate_elements = []
	  list_of_covariate_hyperparameters = []
          for key in self.data.keys():
	      list_of_covariate_elements.append(self.load_covariate(key))
	      list_of_covariate_hyperparameters.append(self.load_hyperparameters(key))

	  return list_of_covariate_elements, list_of_covariate_hyperparameters

      def load_covariate(self, covariate):
          """Load a single covariate object"""

	  [ modulename, classname ] = self.data[covariate]['element']['python_class'].rsplit('.', 1)
	  run_module = importlib.import_module(modulename)
	  element_class = getattr(run_module, classname)
	  if 'parameters' in self.data[covariate]['element'].keys():
	    instance = element_class(**self.data[covariate]['element']['parameters'])
	    if(hasattr(instance,'load')):
		instance.load()
	    return instance
	  else:
	    return element_class()

      def load_hyperparameters(self, covariate):
	  """Load a collection of hyperparameters objects"""

	  list_of_hyperparameters_objects = []
	  for value in self.data[covariate]['hyperparameters'].values():
	      list_of_hyperparameters_objects.append(CovariateHyperparameters(np.log(value)))

	  if len(list_of_hyperparameters_objects) > 1:
	    return CombinationHyperparameters(list_of_hyperparameters_objects)
	  else:
	    return list_of_hyperparameters_objects[0]

	      