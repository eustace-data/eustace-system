"""Provide top-level analysis system for advanced standard method."""

import numpy
import os
import psutil
import sys

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import Regridder

from eustace.outputformats.definitions import GLOBAL_FIELD_OUTPUT_FLAGS
from eustace.outputformats.definitions import CLIMATOLOGY_FRACTION_FILL_VALUE

from eustace.timeutils import epoch

class AnalysisSystemInputLoader(object):
    """Abstract interface to provide data to an analysis system."""

    def load_observation_structure(self, observable, time_index, log=sys.stdout):
        """Load data for observable at specified time index.  Should return an ObservationStructure instance."""
        
        raise NotImplementedError

    def datetime_at_time_index(self, time_index):
        """The time index model that converts time indices to datetime objects."""

        raise NotImplementedError

class AnalysisSystem(object):
    """
    The analysis system is created from analysis components and allows creation of
    ComponentSolution instances which can be plied with data, and from that provide a solution.
    """

    def __init__(self, components, observable, log=sys.stdout):

        self.components = components
        self.observable = observable
        self.log = log
        self.field_flags = GLOBAL_FIELD_OUTPUT_FLAGS
        self.output_flags = ['POINTWISE', 'GRID_CELL_AREA_AVERAGE']

    def update(self, inputloaders, time_indices):
        """
        Update the whole analysis system using data from specified time indices,
        as loaded by the given list of AnalysisSystemInputLoader instances."""
        
        # Iterate over components to provide them with data
        for component_index, component in enumerate(self.components):

            self.update_component(inputloaders, component_index, time_indices)

    def update_component(self, inputloaders, component_index, time_indices):

        # Iterate over time indices
        for time_index in time_indices:

            self.update_component_time(inputloaders, component_index, time_index)

        # Update after data ingestion
        self.update_component_solution(component_index)

    def update_component_time(self, inputloaders, component_index, time_index):
        
        # Get current date for information
        current_date = inputloaders[0].datetime_at_time_index(time_index)

        # Report to user
        self.log.write('Component {0}: {1}: memory {2}MB\n'.format(
                component_index,
                current_date,
                memory_usage_megabytes()))
        self.log.flush()

        # The component being updated
        component = self.components[component_index]
        component_solution = component.component_solution()

        for inputloader in inputloaders:
            
            # Load input
            inputs = inputloader.load_observation_structure(self.observable, time_index, self.log)
    
            number_of_observations = inputs.number_of_observations()
            
            if number_of_observations > 0:
                
                # Sum the expected values provided by the other components
                subtract_offset = numpy.zeros((number_of_observations,))
                for other_index, other_component in enumerate(self.components):

                    # Only sum over other components (not the one about to be computed)
                    if other_index != component_index:                    

                        # Get component solution
                        other_solution = other_component.component_solution()

                        # Evaluate the expected value at given observation locations/times
                        expected_value = other_solution.solution_observation_expected_value(inputs)

                        # Ensure it's a flat array [some elements currently produce a matrix object here]
                        expected_value = numpy.array(expected_value).ravel()
                        
                        # self.log.write('Component {0} mean abs offset: {1}\n'.format(other_index, numpy.abs(expected_value).mean()))

                        # Add to offset
                        subtract_offset += expected_value

                # Process these observations
                component_solution.process_observations(inputs, subtract_offset)

        # Complete this time step
        component_solution.update_time_step()

    def update_component_times(self, inputloaders, component_index, time_indices):
        """Here each elements of inputloaders is a list of input loaders corresponding to observations at a time index in time_indices"""
        
        # TODO preload space time components if possible
        
        # TODO for each inputloader and time index compute a measurement update and append to the 'measurement' (the projection onto model for information vector and precision matrix)
        
        # TODO include similar to update_at_time for disk write for whole time period for space time model. Keep both in code with dummy methods that don't do anything in the component storage classes where not relevant.
        
        # The component being updated
        component = self.components[component_index]
        component_solution = component.component_solution()
        
        for time_index, inputloaders_at_time in zip(time_indices, inputloaders):
            # Get current date for information
            current_date = inputloaders_at_time[0].datetime_at_time_index(time_index)

            # Report to user
            self.log.write('Component {0}: {1}: memory {2}MB\n'.format(
                    component_index,
                    current_date,
                    memory_usage_megabytes()))
            self.log.flush()

            # Reset the expected date for observations in preparation for new time_index 
            component_solution.measurement_time_index = None

            for inputloader in inputloaders_at_time:
                # Load input
                inputs = inputloader.load_observation_structure(self.observable, time_index, self.log)
        
                number_of_observations = inputs.number_of_observations()
                
                if number_of_observations > 0:
                    
                    # Sum the expected values provided by the other components
                    subtract_offset = numpy.zeros((number_of_observations,))
                    for other_index, other_component in enumerate(self.components):

                        # Only sum over other components (not the one about to be computed)
                        if other_index != component_index:                    

                            # Get component solution
                            other_solution = other_component.component_solution() # TODO move this out of for loops so that the component solutions are only initialised once and disk reads can be minimised

                            # Evaluate the expected value at given observation locations/times
                            expected_value = other_solution.solution_observation_expected_value(inputs)

                            # Ensure it's a flat array [some elements currently produce a matrix object here]
                            expected_value = numpy.array(expected_value).ravel()
                            
                            # self.log.write('Component {0} mean abs offset: {1}\n'.format(other_index, numpy.abs(expected_value).mean()))

                            # Add to offset
                            subtract_offset += expected_value

                    # Process these observations
                    component_solution.process_observations(inputs, subtract_offset)

            if isinstance(component, SpatialComponent):
                # Spatial components update on each substep
                component_solution.update_time_step()
                print component_solution
                for key, value in component_solution.__dict__.iteritems():
                    print key, value
                
        if isinstance(component, SpaceTimeComponent):
            # Space time components update for all times after processing all observations
            component_solution.update_time_step()

    def update_component_solution(self, component_index):

        self.log.write('Component {0}: solving\n'.format(component_index))
        self.log.flush()
        self.components[component_index].component_solution().update()
         
    def evaluate_projected_sample(self, outputstructure):
        """
        Project sample onto outputgrid: only pointwise approach is used. Future development should incorporate it into methods that compute grid cell area averages
        """

        # Check all components have drawn samples of the same size
        list_of_sample_sizes = []
        for component in self.components:
	    list_of_sample_sizes.append(component.component_solution().sample_size)
	if len(set(list_of_sample_sizes)) > 1:
	    raise ValueError('Not all model components have the same sample size!')
	else:
	    sample_size = list_of_sample_sizes[0]
	    
        result_expected_value = numpy.zeros((outputstructure.number_of_observations(),sample_size))

        # Iterate over components and evaluate at
        for component_index, component in enumerate(self.components):

            # Get component solution
            component_solution = component.component_solution()

	    # Evaluate the expected value at given observation locations/times
            component_expected_value = component_solution.solution_observation_projected_sample(outputstructure)

            # Add to state
            result_expected_value += component_expected_value

        return result_expected_value

    def evaluate_expected_value(self, field, outputstructure, flag, grid_points_per_cell = [1,1], blocking = 100):
        """
        Evaluate expected value the locations int the given output structure.
        This can be used to produce gridded output by using a OutputRectilinearGridStructure instance
        to configure the required structure.	
        """

        if field not in self.field_flags:
            message = 'Field flag "{}", not recognized: available flags are {}'.format(field, self.field_flags)
            raise ValueError(message)

        if flag not in self.output_flags:
            message = 'Output flag "{}" not recognized: available flags are {}'.format(flag, self.output_flags)
            raise ValueError(message)
        
        if flag == self.output_flags[0]:
            return self.evaluate_pointwise_expected_value(field, outputstructure)
        else:
            return self.evaluate_grid_cell_average_expected_value(field, outputstructure, grid_points_per_cell, blocking)

    def evaluate_pointwise_expected_value(self, field, outputstructure):
        """
        Evaluate expected value the locations int the given output structure.
        This can be used to produce gridded output by using a OutputRectilinearGridStructure instance
        to configure the required structure.
        This method should be deprecated in future.
        """

        # Initialise result mean
        result_expected_value = numpy.zeros((outputstructure.number_of_observations(),))

        # Iterate over components and evaluate at
        for component_index, component in enumerate(self.components):

            # Get component solution
            component_solution = component.component_solution()

            # Evaluate the expected value at given observation locations/times
            if field == self.field_flags[0]:
                component_expected_value = component_solution.solution_observation_expected_value(outputstructure)
            elif field == self.field_flags[1]:
                component_expected_value = numpy.square(component_solution.solution_observation_expected_uncertainties(outputstructure))
            else:
                message = '{} should not be computed as a point like specific output, because it does not represent field uncertainties'.format(field)
                raise ValueError(message)

            # Ensure it's a flat array [some elements currently produce a matrix object here]
            component_expected_value = numpy.array(component_expected_value).ravel()

            # Add to state
            result_expected_value += component_expected_value
            
        if field == self.field_flags[0]:
            return result_expected_value
        else:
            return numpy.sqrt(result_expected_value)
    
    def evaluate_grid_cell_average_expected_value(self, field, outputstructure, grid_points_per_cell, blocking):
        """
        Evaluate grid-cell area averga expected value into the given output structure.
        """

        # Initialise result mean
        result_expected_value = numpy.zeros((outputstructure.number_of_observations(),))

        # Iterate over components and evaluate at 
        for component_index, component in enumerate(self.components):
            # Get component solution

            component_solution = component.component_solution()
            component_expected_value = self.evaluate_single_component_grid_cell_average_expected_value(field, component_solution, outputstructure, grid_points_per_cell, blocking)
	    
            # Add to state
            result_expected_value += component_expected_value

        if field == self.field_flags[0]:
            return result_expected_value
        elif field == self.field_flags[1]:
            return numpy.sqrt(result_expected_value)
        else:
            message = '{} should not be computed as a regridded like specific output, because it does not represent field uncertainties'.format(field)
            raise ValueError(message)

    def evaluate_single_component_grid_cell_average_expected_value(self, field, component_solution, outputstructure, grid_points_per_cell, blocking):
        """
        Evaluate grid-cell area average for a single model component
        """
        
        regridder =  Regridder(outputstructure, grid_points_per_cell, blocking)
        
        if field == self.field_flags[0]:
            return regridder.compute_gridded_expected_value(field, component_solution)
        else:
            return numpy.square(regridder.compute_gridded_expected_value(field, component_solution))

    def evaluate_climatology_fraction(self, outputstructure, grid_points_per_cell=[1,1], blocking=100):
        """Compute climatology fraction"""
        
        index = self.check_spatial_components()
        if index == None:
            return numpy.zeros((outputstructure.number_of_observations(),)) + CLIMATOLOGY_FRACTION_FILL_VALUE
        else:
            component_solution = self.components[index].component_solution()
        
        posterior_variance = self.evaluate_single_component_grid_cell_average_expected_value(GLOBAL_FIELD_OUTPUT_FLAGS[1], component_solution, outputstructure, grid_points_per_cell, blocking)
        prior_variance = self.evaluate_single_component_grid_cell_average_expected_value(GLOBAL_FIELD_OUTPUT_FLAGS[2], component_solution, outputstructure, grid_points_per_cell, blocking)
        
        return posterior_variance/prior_variance

    def check_spatial_components(self):
        """Check the model is using just a single spatial component"""

        list_of_spatial_components = []
        
        for index, component in enumerate(self.components):
            if isinstance(component, SpatialComponent):
                list_of_spatial_components.append(index)

        number_of_finds = len(list_of_spatial_components)
        if (number_of_finds < 1)  or (number_of_finds > 1):
            message = '{} Spatial components available: a unique one is needed.\n'.format(number_of_finds)
            sys.stderr.write('WARNING: ' + message)
            return None
        else:
            return list_of_spatial_components[0]
      
def memory_usage_megabytes():
    """Utility to find current memory usage in megabytes."""

    return psutil.Process(os.getpid()).memory_info().rss / (1024*1024)



