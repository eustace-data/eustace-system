"""A space-time component has state parameters which couple space and time.
   This is used for the large-scale and climatology analysis components."""

import numpy

from eustace.analysis.advanced_standard.components.component import Component
from eustace.analysis.advanced_standard.components.component import ComponentStorage
from eustace.analysis.advanced_standard.components.component import ComponentSolution
from eustace.analysis.advanced_standard.components.component import ComponentSolutionStorage
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystem
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_Cholesky
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_CG
from eustace.analysis.advanced_standard.linalg.linearsystem import MeasurementUpdate

import os
import psutil
import sys

class SpaceTimeComponent(Component):

    def __init__(self, storage, solutionstorage, printstats=False, compute_uncertainties = False, method='EXACT', compute_sample=False, sample_size=1, compute_prior_sample=False):

        super(SpaceTimeComponent, self).__init__(storage, solutionstorage)
        self.printstats = printstats
        self.compute_uncertainties = compute_uncertainties
        self.method=method
        self.compute_sample=compute_sample
        self.sample_size=sample_size
        self.compute_prior_sample = compute_prior_sample
        
    def component_solution(self):

        return SpaceTimeComponentSolution(self.storage.element_read(), self.storage.hyperparameters_read(), 
                                          self.solutionstorage, printstats=self.printstats, 
                                          compute_uncertainties=self.compute_uncertainties, method=self.method,
                                          compute_sample=self.compute_sample, sample_size=self.sample_size,
                                          compute_prior_sample = self.compute_prior_sample)

class SpaceTimeComponentSolutionStorage(ComponentSolutionStorage):
    """Provide persistent state information as required."""

    def state_write(self, state):

        raise NotImplementedError

    def state_marginal_std_write(self, marginal_std):
      
        raise NotImplementedError

    def state_sample_write(self, sample):
      
        raise NotImplementedError

    def state_read(self):

        raise NotImplementedError

    def state_marginal_std_read(self):
      
        raise NotImplementedError

    def state_sample_read(self):
      
        raise NotImplementedError
        
    def state_prior_sample_write(self, sample):

        raise NotImplementedError

    def state_prior_sample_read(self):

        raise NotImplementedError

    def measurement_write(self, measurement, time_index):

        raise NotImplementedError

    def measurement_read(self, time_index):

        raise NotImplementedError

    def timeindices_read(self):

        raise NotImplementedError

class SpaceTimeComponentSolution(ComponentSolution):

    def __init__(self, element, hyperparameters, solutionstorage, printstats=False, compute_uncertainties=False , method='EXACT', compute_sample=False, sample_size=1, compute_prior_sample=False):

        super(SpaceTimeComponentSolution, self).__init__(element, hyperparameters, solutionstorage)
        self.measurement = None
        self.measurement_time_index = None
        self.printstats = printstats
        self.compute_uncertainties = compute_uncertainties
        self.method = method
        self.compute_sample = compute_sample
        self.sample_size = sample_size
        self.compute_prior_sample = compute_prior_sample
        
    def process_observations(self, observations, subtract_offset):

        # Design associated with this element
        design = self.element.element_design(observations)

        # Compute relative observations
        relative_observations = observations.observation_vector() - subtract_offset
        # Load current state estimate (or None if we haven't yet computed the first one)
        current_state_estimate = self.solutionstorage.state_read()

        # Nonlinear systems require an initial estimate to establish jacobian
        if self.element.isnonlinear():

            # Check whether this is first iteration (at which current estimate will be None)
            if current_state_estimate is None:
                
                # The first time this is run there will be no state estimate - use initial one
                current_state_estimate = design.design_initial_state()

                # And store it
                self.solutionstorage.state_write(current_state_estimate)

            # Subtract evaluation at linearisation point when working on nonlinear system
            relative_observations -= design.design_function(current_state_estimate)

        # Evaluate jacobian (using linearisation point if this is a nonlinear system)
        design_jacobian = design.design_jacobian(current_state_estimate)
        
        # Compute update
        update = MeasurementUpdate.create_from_system_observations(design_jacobian, observations.observation_precision(), relative_observations)

        # Check that the time index has been reset to None if processing observations from a next time index
        if self.measurement_time_index is None:            
            self.measurement_time_index = observations.time_index()
        else:
            if self.measurement_time_index != observations.time_index():
                raise ValueError('expect observations all from same time index until next call to update_time_step')

        # Initialise or accumulate the measurements
        if self.measurement is None:
            # Initialise if this is the first
            self.measurement = update
        else:
            # Append to others if a measurement update already exists
            self.measurement.append_measurement_update(update)

    def update_time_step(self):

        self.solutionstorage.measurement_write(self.measurement, self.measurement_time_index)
        self.measurement = None
        self.measurement_time_index = None

    def update(self):

        # Report to user
        print 'Pre-posterior increment: memory {0}MB\n'.format(
                memory_usage_megabytes())
        
        # Build system
        spacetimesystem = LinearSystem(self.element.element_prior(self.hyperparameters))

        if self.compute_prior_sample:
            # only compute prior sample, do not solve for posterior

            # Report to user
            print 'Pre-prior sample: memory {0}MB\n'.format(memory_usage_megabytes())
            
            prior_system_solution = LinearSystemSolution_Cholesky( LinearSystem(self.element.element_prior(self.hyperparameters)), printstats=self.printstats)
            
            # Report to user
            print 'Post-prior sample: memory {0}MB\n'.format(memory_usage_megabytes())
            
            variate = numpy.random.normal(0.0, 1.0, (prior_system_solution.number_of_state_parameters, self.sample_size)  )
            sample = prior_system_solution.white_noise_to_posterior_distribution( variate )
            self.solutionstorage.state_prior_sample_write(sample)
        else:
            # Append all data from cache to system and solve posterior
            for time_index in self.solutionstorage.timeindices_read():

                if self.printstats:
                    print '{0}: loading time_index {1}'.format(__name__, time_index)

                # Retrieve from cache
                update = self.solutionstorage.measurement_read(time_index)

                # Append to system
                spacetimesystem.append_measurement_update(update)
            
            # Report to user
            print 'Post-posterior increment: memory {0}MB\n'.format(memory_usage_megabytes())

            # Solve system
            
            # Report to user
            print 'Pre-solution: memory {0}MB\n'.format(memory_usage_megabytes())
            
            solution = LinearSystemSolution_Cholesky(spacetimesystem, printstats=self.printstats)
            
            # Report to user
            print 'Post-solution: memory {0}MB\n'.format(memory_usage_megabytes())
            
            # self.solution = LinearSystemSolution_CG(spacetimesystem)

            # Get expected value
            state_update = solution.maximum_a_posteriori()
                    
            if self.element.isnonlinear():

                # For nonlinear system this is an offset versus linearisation point
                new_state_estimate = state_update + self.solutionstorage.state_read()

            else:

                # Otherwise just use directly
                new_state_estimate = state_update

            # Store new estimate and clear the system
            self.solutionstorage.state_write(new_state_estimate)
            
            if self.compute_uncertainties:
                if self.method == 'APPROXIMATED':
                    marginal_std = solution.approximated_marginal_std()
                    self.solutionstorage.state_marginal_std_write(marginal_std)
                elif self.method == 'EXACT':
                    marginal_std = solution.marginal_std()
                    self.solutionstorage.state_marginal_std_write(marginal_std)

            if self.compute_sample:
                # Report to user
                print 'Pre-posterior sample: memory {0}MB\n'.format(memory_usage_megabytes())
                
                variate = numpy.random.normal(0.0, 1.0, (solution.number_of_state_parameters, self.sample_size)  )
                sample = solution.white_noise_to_posterior_distribution( variate )
                
                self.solutionstorage.state_sample_write(sample)
            
                # Report to user
                print 'Post-posterior sample: memory {0}MB\n'.format(memory_usage_megabytes())


    def merge_updates(self):
        """Combine measurement updates for multiple time indices into a single measurement update"""
        merged_updates = None
        for time_index in self.solutionstorage.timeindices_read():

            if self.printstats:
                print '{0}: loading time_index {1}'.format(__name__, time_index)

            # Retrieve from cache
            update = self.solutionstorage.measurement_read(time_index)

            if merged_updates is None:
                merged_updates = update
            else:
                merged_updates.append_measurement_update(update)

        return merged_updates

        
    def solution_observation_expected_value(self, observationstructure):

        # Design element for function evaluation
        design = self.element.element_design(observationstructure)

        # Load current state estimate (or None if we haven't yet computed the first one)
        current_state_estimate = self.solutionstorage.state_read()

        # Get state estimate (might be None before first iteration)
        if current_state_estimate is None:

            if self.element.isnonlinear():

                # Evaluate nonlinear function at initial state value
                return design.design_function(design.design_initial_state())

            else:

                # Initial state of linear system is zero
                return numpy.zeros((observationstructure.number_of_observations(),))

        else:

            # We have an estimate from previous evaluation - use this
            function_value = design.design_function(current_state_estimate)
            return function_value
            

    def solution_observation_expected_uncertainties(self, observationstructure):
        
        # Design element for function evaluation
        design = self.element.element_design(observationstructure)
        
        # Load current state estimate (or None if we haven't yet computed the first one)
        current_state_estimate = self.solutionstorage.state_marginal_std_read()
        
        # Get state estimate (might be None before first iteration)
        if current_state_estimate is None:
                
                # Initial state of linear system is zero
                return numpy.zeros((observationstructure.number_of_observations(),))
        else:
            
            # We have an estimate from previous evaluation - use this
            #function_value = design.design_function(current_state_estimate)
            function_value = design.design_function(current_state_estimate, rectified=True)
            return function_value

    def solution_observation_projected_sample(self, observationstructure):

        # Design element for function evaluation
        design = self.element.element_design(observationstructure)

        # Load current state estimate (or None if we haven't yet computed the first one)
        current_state_estimate = self.solutionstorage.state_sample_read()

        # Get state estimate (might be None before first iteration)
        if current_state_estimate is None:

            if self.element.isnonlinear():

                # Evaluate nonlinear function at initial state value
                return design.design_function(design.design_initial_state())

            else:

                # Initial state of linear system is zero
                return numpy.zeros((observationstructure.number_of_observations(),self.sample_size))

        else:

            # We have an estimate from previous evaluation - use this
            function_value = design.design_function(current_state_estimate)
            return function_value

def memory_usage_megabytes():
    """Utility to find current memory usage in megabytes."""

    return psutil.Process(os.getpid()).memory_info().rss / (1024*1024)
