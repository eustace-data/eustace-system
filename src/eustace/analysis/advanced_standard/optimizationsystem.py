import sys
import numpy

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.analysissystem import memory_usage_megabytes
from eustace.analysis.advanced_standard.examples.inputloader_rawbinary import  AnalysisSystemInputLoaderRawBinary_OneDay
from eustace.analysis.advanced_standard.fileio import optimization_inputs
from eustace.timeutils import epoch
import dateutil.parser

class OptimizationSystem(AnalysisSystem):
    """
    
    Adds functionality to AnalysisSystem for reading observations from
    input_descriptor dictionaries (time step index dictionaries of input files).
    
    """
    
    def __init__(self, components, observable, log=sys.stdout):
        
        super(OptimizationSystem, self).__init__(components, observable, log=sys.stdout)
    
    def process_inputs(self, input_descriptor, component_index, time_keys):
        """Pre-process observations at specified times for a specified component
        
        Does not solve the system. To preprocess the observations and also
        solve the system run update_component.
        
        """
        
        for time_key in time_keys:
            
            # convert time_key string to days since epoch
            this_time = dateutil.parser.parse(time_key)
            time_index = int(epoch.days_since_epoch(this_time))
            
            # Build inputloaders from list of sources            
            inputloaders = [ AnalysisSystemInputLoaderRawBinary_OneDay(time_index=time_index, **source) for source in input_descriptor[time_key] ]
            
            # Build and store measurement systems for component
            self.update_component_time(inputloaders, component_index, time_index)

    def update_component(self, input_descriptor, component_index, time_keys):
        """Pre-process observations at specified times for a specified component and solve."""
        
        # Generate projection of observations onto model variables and store
        self.process_inputs(input_descriptor, component_index, time_keys)
        
        # Update after data ingestion
        self.update_component_solution(component_index)

    def optimize_component(self, component_index, hyperparameter_storage_file=None):
        """Run optimization for specified component"""
        
        self.components[component_index].component_solution().optimize()
        
        if hyperparameter_storage_file is not None:
            self.components[component_index].storage.hyperparameters.values_to_npy_savefile( hyperparameter_storage_file )
        
    def optimize_component_time(self, inputloaders, component_index, time_index):
        """Run optimization for a single component at a single time index
        
        Replaced by optimize_component and soon to be deprecated.
        
        """

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
                subtract_offset = numpy.zeros((inputs.number_of_observations(),))
                for other_index, other_component in enumerate(self.components):

                    # Only sum over other components (not the one about to be computed)
                    if other_index != component_index:                    

                        # Get component solution
                        other_solution = other_component.component_solution()

                        # Evaluate the expected value at given observation locations/times
                        expected_value = other_solution.solution_observation_expected_value(inputs)
                        print other_index, expected_value
                        # Ensure it's a flat array [some elements currently produce a matrix object here]
                        expected_value = numpy.array(expected_value).ravel()

                        # self.log.write('Component {0} mean abs offset: {1}\n'.format(other_index, numpy.abs(expected_value).mean()))

                        # Add to offset
                        subtract_offset += expected_value

                # Process these observations
                component_solution.process_observations(inputs, subtract_offset)
            
        # Run the optimization
        component_solution.optimize_time_step()

        # Complete this time step
        component_solution.update_time_step()
        
        return component_solution.hyperparameters
    
    def merge_local_parameterisations(local_spde_view_list, global_non_stationary_spde, local_hyperparameter_list):
        # move this into its own module
        
        sigma_design_accumulator = []
        rho_design_accumulator = []
        theta_accumulator = []
        
        sigma_contribution_counter = np.zeros( global_non_stationary_spde.n_latent_variables() )
        
        for (local_spde, hyperparameters) in zip(local_spde_view_list, local_hyperparameter_list):
            
            vertex_map = supertriangle.spde.get_vertex_mapping_matrix().T
            
            local_sigma_design = vertex_map.dot( local_spde.sigma_design_matrix() )
            sigma_contribution_counter += np.int64( local_sigma_design.getnnz(axis = 1)  > 0 )
            
            local_rho_design = vertex_map.dot( local_spde.rho_design_matrix() )
            rho_contribution_counter += np.int64( local_rho_design.getnnz(axis = 1)  > 0 )
            
            sigma_design_accumulator.append( local_sigma_design )
            rho_design_accumulator.append( local_rho_design )
            theta_accumulator.append( hyperparameters )
            
        sigma_design = scipy.sparse.hstack(sigma_design_accumulator)       
        sigma_normaliser = scipy.sparse.diags( 1.0 / sigma_contribution_counter, 0 )
        sigma_design = sigma_normaliser * sigma_design
        
        rho_design = scipy.sparse.hstack(rho_design_accumulator)    
        rho_normaliser = scipy.sparse.diags( 1.0 / rho_contribution_counter, 0 )
        rho_design = sigma_normaliser * rho_design
        
        merged_hyperparameter_values = numpy.concatenate( theta_accumulator )
        
        return merged_hyperparameter_values, sigma_design, rho_design
        
        