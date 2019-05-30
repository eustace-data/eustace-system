"""Base class for analysis components."""

from eustace.analysis.advanced_standard.linalg.linearsystem import MeasurementUpdate

class Component(object):
    """Associate the construction of a component solution (that can do the maths)
       with a suitable storage class for system parameters and hyperparameters."""

    def __init__(self, storage, solutionstorage):

        # This is used to load element definition and hyperparameters
        self.storage = storage

        # This is used to maintain object state during solution phase
        self.solutionstorage = solutionstorage

    def component_solution(self):
        """Construct a ComponentSolution instance that can be updated with observations,
           and provide it with element and hyperparameters."""

        raise NotImplementedError

class ComponentStorage(object):
    """Ability to store element and hyperparameter definition for a component.
       Sublassed by component implementations to maintain state between timesteps
       and iterations."""

    def element_read(self):
        """Retrieve Element instance."""
        
        raise NotImplementedError

    def hyperparameters_read(self):
        """Retrieve current hyperparameters."""

        raise NotImplementedError

class ComponentSolutionStorage(object):
    """Ability to store/load info needed during solution process.
       Will be different for each solution."""

    pass

class ComponentSolution(object):
    """A component solution is the linear system (or array of linear systems) for a component."""

    def __init__(self, element, hyperparameters, solutionstorage):

        self.element = element
        self.hyperparameters = hyperparameters
        self.solutionstorage = solutionstorage

    def process_observations(self, observations, subtract_offset):

        raise NotImplementedError

    def solution_observation_expected_value(self, observationstructure):

        raise NotImplementedError        

    def update_time_step(self):

        raise NotImplementedError

    def update(self):

        raise NotImplementedError

class HyperparameterEvaluation(object):

    def __init__(self, system_mean, prior, initial_hyperparameters, observation_precision):

        self.system_mean = system_mean
        self.observation_precision = observation_precision
        self.prior = prior
        self.hyperparameters = initial_hyperparameters

    def negative_log_liklihood(self):
        """Negative log liklihood at current state (used as cost function to optimise hyperparameters)."""

        pass

    def negative_log_liklihood_derivative(self):
        """Vector of derivatives with respect to each hyperparameter."""

        pass
