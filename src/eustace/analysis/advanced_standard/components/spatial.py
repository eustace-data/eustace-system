"""A spatial component is one which has an independent state vector at every time increment
   This is used for the local analysis component."""

from eustace.analysis.advanced_standard.components.component import Component
from eustace.analysis.advanced_standard.components.component import ComponentStorage
from eustace.analysis.advanced_standard.components.component import ComponentSolution
from eustace.analysis.advanced_standard.components.component import ComponentSolutionStorage
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystem
from eustace.analysis.advanced_standard.linalg.linearsystem import MeasurementUpdate
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_Cholesky
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_CG

import numpy

class SpatialComponent(Component):

    def __init__(self, storage, solutionstorage, compute_uncertainties = False, method='EXACT', compute_sample=False, sample_size=1, compute_prior_sample=False):

        super(SpatialComponent, self).__init__(storage, solutionstorage)
        self.compute_uncertainties = compute_uncertainties
        self.method = method
        self.compute_sample=compute_sample
        self.sample_size=sample_size
        self.compute_prior_sample = compute_prior_sample
	
    def component_solution(self):

        return SpatialComponentSolution(self.storage.element_read(), self.storage.hyperparameters_read(), 
                                        self.solutionstorage, compute_uncertainties=self.compute_uncertainties, method=self.method,
                                        compute_sample=self.compute_sample, sample_size=self.sample_size, compute_prior_sample=self.compute_prior_sample)

class SpatialComponentSolutionStorage(ComponentSolutionStorage):
    """Extend storage class to maintain persistent partial solution information
       at each timestep."""

    def partial_state_write(self, state, time_index):

        raise NotImplementedError

    def partial_state_marginal_std_write(self, marginal_std, time_index):
      
        raise NotImplementedError

    def partial_state_sample_write(self, sample, time_index):
      
        raise NotImplementedError

    def partial_state_read(self, time_index):

        raise NotImplementedError

    def partial_state_marginal_std_read(self, time_index):

        raise NotImplementedError

    def partial_state_sample_read(self, time_index):

        raise NotImplementedError
        
    def partial_state_prior_sample_write(self, time_index):

        raise NotImplementedError

    def partial_state_prior_sample_read(self, time_index):

        raise NotImplementedError


class SpatialComponentSolution(ComponentSolution):

    def __init__(self, element, hyperparameters, solutionstorage, compute_uncertainties=False, method='EXACT', compute_sample=False, sample_size=1, compute_prior_sample=False):

        super(SpatialComponentSolution, self).__init__(element, hyperparameters, solutionstorage)

        # Pre-compute prior
        self.prior = self.element.element_prior(self.hyperparameters)

        # Cache for linear system (current processing time only)
        self.localsystem_cache = None
        self.localsystem_time_index = None
        self.compute_uncertainties = compute_uncertainties
        self.method = method
        self.compute_sample = compute_sample
        self.sample_size = sample_size
        self.compute_prior_sample = compute_prior_sample

    def process_observations(self, observations, subtract_offset):

        # Build the system
        design = self.element.element_design(observations)
        design_jacobian = design.design_jacobian(currentstate=None)

        # Compute relative observations
        relative_observations = observations.observation_vector() - subtract_offset

        # Build the measurement update
        update = MeasurementUpdate.create_from_system_observations(design_jacobian, observations.observation_precision(), relative_observations)

        # Build linear system for this time step if it doesn't exist
        if self.localsystem_time_index is None:

            # Need new system for this time step
            self.localsystem_time_index = observations.time_index()
            self.localsystem_cache = LinearSystem(self.prior)

        else:

            if self.localsystem_time_index != observations.time_index():
                print self.localsystem_time_index, observations.time_index()
                raise ValueError('expect observations all from same time index until next call to update_time_step')                

        # Append this measurement update
        self.localsystem_cache.append_measurement_update(update)

    def update_time_step(self):
        """Perform update after processing one time step of observations."""
        
        if self.localsystem_cache is None or self.localsystem_time_index is None:
            print "Linear system is not initialised. Component solution will not be processed"
        else:
            # Solve for this time step
            solution = LinearSystemSolution_Cholesky(self.localsystem_cache)
            # solution = LinearSystemSolution_CG(self.localsystem_cache)
            self.solutionstorage.partial_state_write(solution.maximum_a_posteriori(), self.localsystem_time_index)
            if self.compute_uncertainties:
                if self.method == 'APPROXIMATED':
                    self.solutionstorage.partial_state_marginal_std_write(solution.approximated_marginal_std(), self.localsystem_time_index)
                elif self.method == 'EXACT':
                    self.solutionstorage.partial_state_marginal_std_write(solution.marginal_std(), self.localsystem_time_index)
                    
            if self.compute_sample:
                variate = numpy.random.normal(0.0, 1.0, (solution.number_of_state_parameters, self.sample_size)  )
                sample = solution.white_noise_to_posterior_distribution( variate )
                self.solutionstorage.partial_state_sample_write(sample, self.localsystem_time_index)
                
                if self.compute_prior_sample:
                    prior_system_solution = LinearSystemSolution_Cholesky( LinearSystem(self.element.element_prior(self.hyperparameters)))
                    sample = numpy.squeeze( prior_system_solution.white_noise_to_posterior_distribution( variate ) )
                    self.solutionstorage.partial_state_prior_sample_write(sample, self.localsystem_time_index)
                    #self.solutionstorage.partial_state_precision_write(solution.posterior_precision, self.localsystem_time_index) # TODO remove this line to test size of saved precision matrices
            
        self.localsystem_cache = None
        self.localsystem_time_index = None

    def update(self):
        """Perform any update after processing all observations.
        In SpatialComponentSolution this performs no additional operations."""
        
        pass

    def solution_observation_expected_value(self, observations):

        # Time index to look up
        time_index = observations.time_index()

        # Attempt to retrieve solution
        partial_state = self.solutionstorage.partial_state_read(time_index)

        # Return zeros if can't find it
        if partial_state is None:
            return numpy.zeros((observations.number_of_observations(),))

        # Build new design for this observation structure
        design = self.element.element_design(observations)

        # Evaluate at the computed state
        return design.design_function(partial_state)
        
    def solution_observation_expected_uncertainties(self, observations):

        # Time index to look up
        time_index = observations.time_index()

        # Attempt to retrieve solution
        partial_state = self.solutionstorage.partial_state_marginal_std_read(time_index)
        
        # Return zeros if can't find it
        if partial_state is None:
            return numpy.zeros((observations.number_of_observations(),))

        # Build new design for this observation structure
        design = self.element.element_design(observations)

        # Evaluate at the computed state: correct with absolute value to compensate zero crossing from harmonic time components
        return design.design_function(partial_state, rectified=True)
        
    def solution_prior_std(self, number_of_samples=1000):
        """Extract marginal std from prior distribution"""
        
        linear_system = LinearSystem(self.prior)
        solution = LinearSystemSolution_Cholesky(linear_system, printstats=False)
        return solution.approximated_marginal_std(number_of_samples)

    def solution_observation_prior_uncertainties(self, observations):
        """Project prior marginal std onto observation grid, usign desgin matrix"""
        
        # Build new design for this observation structure
        design = self.element.element_design(observations)

        # Evaluate at the computed state
        return design.design_function(self.solution_prior_std())

    def solution_observation_projected_sample(self, observations):

        # Time index to look up
        time_index = observations.time_index()

        # Attempt to retrieve solution
        partial_state = self.solutionstorage.partial_state_sample_read(time_index)

        # Return zeros if can't find it
        if partial_state is None:
            return numpy.zeros((observations.number_of_observations(),self.sample_size))

        # Build new design for this observation structure
        design = self.element.element_design(observations)

        # Evaluate at the computed state
        return design.design_function(partial_state)
