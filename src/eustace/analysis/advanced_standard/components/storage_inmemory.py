
from component import ComponentStorage
from spatial import SpatialComponentSolutionStorage
from spacetime import SpaceTimeComponentSolutionStorage
import numpy

class ComponentStorage_InMemory(ComponentStorage):
    """An in-memory implementation of element and hyperparameters."""

    def __init__(self, element, hyperparameters):

        self.element = element
        self.hyperparameters = hyperparameters

    def element_read(self):

        return self.element

    def hyperparameters_read(self):

        return self.hyperparameters

class SpaceTimeComponentSolutionStorage_InMemory(SpaceTimeComponentSolutionStorage):

    def __init__(self):

        super(SpaceTimeComponentSolutionStorage_InMemory, self).__init__()
        self.state = None
        self.state_marginal_std = None
        self.state_sample = None
        self.state_prior_sample = None
        self.measurement = { }

    def state_write(self, state):

        self.state = state

    def state_marginal_std_write(self, state_marginal_std):

        self.state_marginal_std = state_marginal_std

    def state_sample_write(self, sample):

        self.state_sample = sample

    def state_prior_sample_write(self, sample):

        self.state_prior_sample = sample

    def state_read(self):

        return self.state

    def state_marginal_std_read(self):

        return self.state_marginal_std

    def state_sample_read(self):

        return self.state_sample
        
    def state_prior_sample_read(self):

        return self.state_prior_sample

    def measurement_write(self, measurement, time_index):

        self.measurement[time_index] = measurement

    def measurement_read(self, time_index):

        return self.measurement[time_index]

    def timeindices_read(self):

        return self.measurement.keys()

class SpatialComponentSolutionStorage_InMemory(SpatialComponentSolutionStorage):

    def __init__(self):

        super(SpatialComponentSolutionStorage_InMemory, self).__init__()
        self.state_at_time = { }
        self.state_marginal_std_at_time = { }
        self.state_sample_at_time = {}
        self.state_prior_sample_at_time = {}
        self.measurement_at_time = { }
        
    def partial_state_write(self, state, time_index):

        self.state_at_time[time_index] = state

    def partial_state_read(self, time_index):

        try:

            return self.state_at_time[time_index]

        except KeyError:

            return None

    def partial_state_marginal_std_write(self, state_marginal_std, time_index):

        self.state_marginal_std_at_time[time_index] = state_marginal_std

    def partial_state_marginal_std_read(self, time_index):

        try:

            return self.state_marginal_std_at_time[time_index]

        except KeyError:

            return None

    def partial_state_sample_write(self, sample, time_index):

        self.state_sample_at_time[time_index] = sample

    def partial_state_sample_read(self, time_index):

        try:

            return self.state_sample_at_time[time_index]

        except KeyError:

            return None
            
    def partial_state_prior_sample_write(self, sample, time_index):

        self.state_prior_sample_at_time[time_index] = sample

    def partial_state_prior_sample_read(self, time_index):

        try:

            return self.state_prior_sample_at_time[time_index]

        except KeyError:

            return None

    def measurement_write(self, measurement, time_index):

        self.measurement_at_time[time_index] = measurement

    def measurement_read(self, time_index):

        return self.measurement_at_time[time_index]

    def timeindices_read(self):

        return self.measurement_at_time.keys()