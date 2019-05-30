"""A component for multiple replicates of a spatial process that are processed
in a delayed fashion. The component splits obervation updates into two steps:
first the measurements are processed for each replicate and then, once all are
processed, each replicate can the be solved. 

This differs from the regular SpaceTimeComponentSolution in having one linear 
system solve per time step, rather than having one solve for a single system 
that encompases all time steps. This approach separates observation 
preprocessing and the delayed system solves and allows the pre-processed 
observations to be reused multiple times in component hyperparameter 
optimisation.

"""

from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponentSolutionStorage
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponentSolution
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystem
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_Cholesky
from eustace.analysis.advanced_standard.linalg.linearsystem import LinearSystemSolution_CG
from eustace.analysis.advanced_standard.linalg.linearsystem import MeasurementUpdate

from eustace.analysis.advanced_standard.linalg.costfunction import MarginalGaussianReplicateCF_InMemory
from eustace.analysis.advanced_standard.linalg.optimize import Optimizer

import numpy 

class DelayedSpatialComponent(SpaceTimeComponent):
    
    def __init__(self, storage, solutionstorage, printstats=False, compute_uncertainties = False, method='EXACT', compute_sample=False, sample_size=1):

        super(DelayedSpatialComponent, self).__init__(storage, solutionstorage, printstats = printstats, compute_uncertainties = compute_uncertainties, method=method)
        
    def component_solution(self):
        
        return DelayedSpatialComponentSolution(self.storage.element_read(), self.storage.hyperparameters_read(), 
                                               self.solutionstorage, printstats=self.printstats, 
                                               compute_uncertainties=self.compute_uncertainties, method=self.method,
                                               compute_sample=self.compute_sample, sample_size=self.sample_size)
        
class DelayedSpatialComponentSolutionStorage(SpaceTimeComponentSolutionStorage):
    """Template super class for storage replicate measurement/state storage"""
    
    def cost_function_read(self):
        
        raise NotImplementedError
    
    def cost_function_write(self):
        
        raise NotImplementedError

class DelayedSpatialComponentSolution(SpaceTimeComponentSolution):
    
    def __init__(self, element, hyperparameters, solutionstorage, printstats=False, compute_uncertainties=False , method='EXACT', compute_sample=False, sample_size=1):

        super(DelayedSpatialComponentSolution, self).__init__(element, hyperparameters, solutionstorage)
        self.measurement = None
        self.measurement_time_index = None
        self.printstats = printstats
        self.compute_uncertainties = compute_uncertainties
        self.method = method
        self.compute_sample = compute_sample
        self.sample_size = sample_size
    
    def process_observations(self, observations, subtract_offset):
        
        # Get index for current time from observations
        if self.measurement_time_index is None:

            # Initialise if this is the first
            self.measurement_time_index = observations.time_index()
            
        else:

            if self.measurement_time_index != observations.time_index():
                raise ValueError('expect observations all from same time index until next call to update_time_step')

        
        # Design associated with this element
        design = self.element.element_design(observations)

        # Compute relative observations
        relative_observations = observations.observation_vector() - subtract_offset
        
        # Load current state estimate (or None if we haven't yet computed the first one)
        current_state_estimate = self.solutionstorage.partial_state_read( self.measurement_time_index )

        # Nonlinear systems require an initial estimate to establish jacobian
        if self.element.isnonlinear():

            # Check whether this is first iteration (at which current estimate will be None)
            if current_state_estimate is None:
                
                # The first time this is run there will be no state estimate - use initial one
                current_state_estimate = design.design_initial_state()

                # And store it
                self.solutionstorage.partial_state_write(current_state_estimate, self.measurement_time_index)

            # Subtract evaluation at linearisation point when working on nonlinear system
            relative_observations -= design.design_function(current_state_estimate)

        # Evaluate jacobian (using linearisation point if this is a nonlinear system)
        design_jacobian = design.design_jacobian(current_state_estimate)
                    
        # Compute update
        update = MeasurementUpdate.create_from_system_observations(design_jacobian, observations.observation_precision(), relative_observations)

        if self.measurement is None:
            self.measurement = update
        else:
            # Append to others
            self.measurement.append_measurement_update(update)

    def update(self):

        # Append all data from cache to system
        for time_index in self.solutionstorage.timeindices_read():

            # Build system
            spacetimesystem = LinearSystem(self.element.element_prior(self.hyperparameters))

            if self.printstats:
                print '{0}: loading time_index {1}'.format(__name__, time_index)

            # Retrieve from cache
            update = self.solutionstorage.measurement_read(time_index)

            # Append to system
            spacetimesystem.append_measurement_update(update)

            # Solve system
            solution = LinearSystemSolution_Cholesky(spacetimesystem, printstats=self.printstats)

            # Get expected value
            state_update = numpy.squeeze(solution.maximum_a_posteriori())
                    
            if self.element.isnonlinear():

                # For nonlinear system this is an offset versus linearisation point
                new_state_estimate = state_update + numpy.squeeze(self.solutionstorage.partial_state_read())

            else:

                # Otherwise just use directly
                new_state_estimate = state_update

            # Store new estimate and clear the system
            self.solutionstorage.partial_state_write(new_state_estimate, time_index)
            
            if self.compute_uncertainties:
                if self.method == 'APPROXIMATED':
                    marginal_std = numpy.squeeze(solution.approximated_marginal_std())
                    self.solutionstorage.partial_state_marginal_std_write(marginal_std, time_index)
                elif self.method == 'EXACT':
                    marginal_std = numpy.squeeze(solution.marginal_std())
                    self.solutionstorage.partial_state_marginal_std_write(marginal_std, time_index)
                    
            if self.compute_sample:
                sample = numpy.squeeze(solution.random_sample(self.sample_size))
                self.solutionstorage.state_sample_write(sample, time_index)
  
    def solution_observation_expected_value(self, observations):

        # Time index to look up
        time_index = observations.time_index()

        # Attempt to retrieve solution
        state = self.solutionstorage.partial_state_read(time_index)

        # Return zeros if can't find it
        if state is None:
            return numpy.zeros((observations.number_of_observations(),))

        # Build new design for this observation structure
        design = self.element.element_design(observations)

        # Evaluate at the computed state
        return design.design_function(state)
        
    def solution_observation_expected_uncertainties(self, observations):

        # Time index to look up
        time_index = observations.time_index()

        # Attempt to retrieve solution
        state = self.solutionstorage.partial_state_marginal_std_read(time_index)
        
        # Return zeros if can't find it
        if state is None:
            return numpy.zeros((observations.number_of_observations(),))

        # Build new design for this observation structure
        design = self.element.element_design(observations)

        # Evaluate at the computed state: correct with absolute value to compensate zero crossing from harmonic time components
        return design.design_function(state, rectified=True)

    def solution_observation_projected_sample(self, observations):

        # Time index to look up
        time_index = observations.time_index()

        # Attempt to retrieve solution
        state = self.solutionstorage.state_sample_read(time_index)
        
        # Return zeros if can't find it
        if state is None:
            return numpy.zeros((observations.number_of_observations(),self.sample_size))

        # Build new design for this observation structure
        design = self.element.element_design(observations)

        # Return the state samples projected onto observation locations
        return design.design_function(state)

    def optimize(self):
        # This version loads all measurement and state information into memory before passing to optimization algorithm
        
        # Accumulate the required information from storage for each time step
        information_increments = []
        precision_increments = []
        prior_means = []
        set_points = []
        
        for time_index in self.solutionstorage.timeindices_read():
            
            # Retrieve measurements from cache
            measurement_update_at_time = self.solutionstorage.measurement_read(time_index)
            information_increments.append( numpy.squeeze(measurement_update_at_time.information_vector) )
            precision_increments.append(measurement_update_at_time.information_precision)
            
            # Retrieve measurements from cache
            state_at_time = self.solutionstorage.partial_state_read(time_index)
            set_points.append(state_at_time)
            prior_means.append(numpy.zeros(state_at_time.shape))
        
        # Setup the MarginalGaussianCF cost function object and its arguments
        initial_hyperparameters = self.hyperparameters.get_array()
        
        # Setup wrappers for the prior calculation and its derivative to standardise the argument list between the two
        def Q_function(hyperparameter_values):
            self.hyperparameters.set_array(hyperparameter_values)
            return self.element.element_prior(self.hyperparameters).prior_precision()
        
        def dQ_function(hyperparameter_values, parameter_index):
            self.hyperparameters.set_array(hyperparameter_values)
            return self.element.element_prior(self.hyperparameters).prior_precision_derivative(parameter_index)
        
        # Setup the MarginalGaussianReplicateCF_InMemory cost function object and its arguments
        Qargs = ()
        print information_increments
        cost_function = MarginalGaussianReplicateCF_InMemory( Q_function, dQ_function)
        cost_function_args = (set_points, prior_means, information_increments, precision_increments, Qargs)
        
        # Instantiate the scipy.optimize wrapper
        initial_hyperparameters = self.hyperparameters.get_array()
        optimizer = Optimizer(initial_hyperparameters, cost_function, cost_function_args = cost_function_args, verbose = True)

        # Run the optimization
        optimizer.fit()

        # If sucessfully converged then update the hyperparameter values and store - > this conditional statement should be moved into the optimizer object module
        if optimizer.result.success == True:
            #raise NotImplementedError('Hyperparameter storage not implemented for optimization')
            self.hyperparameters.set_array(optimizer.result.x)
        else:
            #raise NotImplementedError('Hyperparameter storage not implemented for return to defaults')
            self.hyperparameters.set_array(initial_hyperparameters)
        
        print "optim_path:"
        print optimizer.optim_path