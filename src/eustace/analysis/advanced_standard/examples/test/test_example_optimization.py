"""Test component splitting and subregioning in example_optimization"""

from eustace.analysis.advanced_standard.examples import example_eustace
from eustace.analysis.advanced_standard.examples import example_optimization
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory, SpatialComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.spatialdelayed import DelayedSpatialComponent
from eustace.analysis.advanced_standard.elements.combination import CombinationElement, CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement, LocalHyperparameters
import unittest
import numpy

import tempfile
import os.path

class TestDirectory():
    
    def __init__(self, suffix=None, prefix=None, dir=None):
        self.test_directory = tempfile.mkdtemp()
        print self.test_directory
    
    def temporary_filename(self, filename):
        new_filename =  os.path.join(self.test_directory, filename)
        print new_filename
        return new_filename

class TestSystem():
    pass


class TestExampleOptimization(unittest.TestCase):
    
    def test_split_states_time(self):
        
        global_biases_group_list = ['global bias',]
        
        full_component_definition = example_eustace.LocalDefinition(bias_terms = True, global_biases_group_list = global_biases_group_list)
        full_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
        full_component = DelayedSpatialComponent(full_component_definition, full_component_storage_solution)
        
        local_sub_component_definition = example_optimization.PureLocalComponentDefinition()
        local_sub_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
        local_sub_component = DelayedSpatialComponent(local_sub_component_definition, local_sub_component_storage_solution)
        
        bias_sub_component_definition = example_optimization.PureBiasComponentDefinition(bias_terms = True, global_biases_group_list = global_biases_group_list)
        bias_sub_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
        bias_sub_component = DelayedSpatialComponent(bias_sub_component_definition, bias_sub_component_storage_solution)
        
        element_optimisation_flags = [True, False]
        time_key = 0
        
        state_length = full_component.storage.element_read().element_prior( full_component.storage.hyperparameters  ).prior_number_of_state_parameters()
        random_state = numpy.random.normal(0.0, 1.0, state_length)        
        full_component.solutionstorage.partial_state_write(random_state, time_key)
        
        example_optimization.split_states_time( full_component, local_sub_component, bias_sub_component, element_optimisation_flags, time_key )
        
        full_state = full_component.solutionstorage.partial_state_read( time_key )
        
        (local_target_state, bias_target_state) = full_component.storage.element_read().element_prior(full_component.storage.hyperparameters_read()).element_states(full_state)
        
        numpy.testing.assert_almost_equal( local_target_state, local_sub_component.solutionstorage.partial_state_read(time_key) )
        numpy.testing.assert_almost_equal( bias_target_state, bias_sub_component.solutionstorage.partial_state_read(time_key) )
        
        
    def test_extract_local_view_states_time(self):
        
        # Define model components for full global and local views
        global_sub_component_definition = example_optimization.PureLocalComponentDefinition()
        global_sub_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
        global_component = DelayedSpatialComponent(global_sub_component_definition, global_sub_component_storage_solution)
        
        # Assign a random state
        time_key = 0
        
        state_length = global_component.storage.element_read().element_prior( global_component.storage.hyperparameters  ).prior_number_of_state_parameters()
        random_state = numpy.random.normal(0.0, 1.0, state_length)        
        global_component.solutionstorage.partial_state_write(random_state, time_key)
        
        # Define subregions and extract their states
        neighbourhood_level = 0
        
        n_subregions = global_component.storage.element_read().combination[0].spde.n_triangles_at_level(neighbourhood_level)
        
        # split states for each local model
        setup = example_eustace.LocalSetup()
        
        subregion_component_list = []
        for region_index in range(n_subregions):
            print region_index
            view_flags = [True,]
            region_component_definition = example_optimization.LocalViewDefinition( neighbourhood_level, region_index )
            region_component_storage_solution = SpatialComponentSolutionStorage_InMemory()
            region_component = DelayedSpatialComponent(region_component_definition, region_component_storage_solution)
            
            print "extracting state"
            example_optimization.extract_local_view_states_time( global_component, region_component, view_flags, time_key )
            
            subregion_component_list.append(region_component)
        
        # reconstruct the full state by averaging over states for local splits
        state_accumulator = numpy.zeros( state_length )
        n_contributors_accumulator = numpy.zeros( state_length )
        
        for region_index in range(n_subregions):
            region_state_vector = subregion_component_list[region_index].solutionstorage.partial_state_read( time_key )
            full_state_indices = subregion_component_list[region_index].storage.element_read().combination[0].spde.active_vertex_indices
            state_accumulator[full_state_indices] += region_state_vector
            n_contributors_accumulator[full_state_indices] += 1.0
            
        reconstructed_state = state_accumulator / n_contributors_accumulator
        
        # assert that the reconstructed state matches the original full component state vector
        numpy.testing.assert_almost_equal( reconstructed_state, global_component.solutionstorage.partial_state_read(time_key) )

    def test_NonStationaryLocalDefinition(self):
        """Check that component can be constructed reading hyperparaters from disk"""
        
        # setup a hyperparameter file
        test_directory = TestDirectory()
        hyperparameter_file = test_directory.temporary_filename( 'test_hyperparameter_file.npy' )
        example_optimization.initialise_local_hyperparameter_file(hyperparameter_file)
        
        # setup a model component with hyperparameters read from the hyperparameter_file
        component_definition = example_eustace.NonStationaryLocalDefinition(bias_terms = False, global_biases_group_list = [], local_hyperparameter_file = hyperparameter_file)
        
        # assert that the number of hyperparameters is equal to twice the number of latent variables for the local spde
        n_hyperparameters = len(component_definition.hyperparameters.get_array())
        self.assertEqual(component_definition.element.combination[0].spde.n_latent_variables()*2, n_hyperparameters)
        
        my_hyperparameter_values = component_definition.hyperparameters.elementparameters[0].get_array()
        
        # assert that the hyperparameter array contains identical values for first half of array
        self.assertTrue( all(x == my_hyperparameter_values[0] for x in my_hyperparameter_values[:n_hyperparameters/2]) )
        
        # assert that the hyperparameter array contains identical values for second half of array 
        self.assertTrue( all(x == my_hyperparameter_values[n_hyperparameters/2] for x in my_hyperparameter_values[n_hyperparameters/2:]) )
        
    def test_merge_local_parameterisations(self):
        
        test_directory = TestDirectory()
        
        setup = example_eustace.LocalSetup()
        n_triangulation_divisions = setup.n_triangulation_divisions
        neighbourhood_level = 0
        
        n_regions = 20
        for level in range(neighbourhood_level):
            n_regions = n_regions * 4
        
        hyperparameter_filenames = []
        hyperparameter_store = []
        output_filename = test_directory.temporary_filename( "temp_hyperparameters.%d.%d.npy" % ( n_triangulation_divisions, neighbourhood_level ) )
        for region_index in range(n_regions):
            # save temporary hyperparmeter files to disk
            
            hyperparameters = numpy.log(numpy.random.uniform(1,2,2))
            hyperparameter_store.append(hyperparameters)
            local_filename = test_directory.temporary_filename( "temp_hyperparameters.%d.%d.%d.npy" % ( n_triangulation_divisions, neighbourhood_level, region_index ) )
            numpy.save(local_filename, hyperparameters)
            
            hyperparameter_filenames.append(local_filename)
        
        # merge parameterisations
        example_optimization.merge_local_parameterisations(n_triangulation_divisions, neighbourhood_level, hyperparameter_filenames, output_filename)
        
        merged_hyperparameters = numpy.load(output_filename)
        print merged_hyperparameters
        
        # check that that the number of hyperparameter entries equalling each regional (sigma, rho) equals the number of interior points in a region
        
        # find the number of interior points in a region
        level_difference = n_triangulation_divisions - neighbourhood_level
        n_edges_of_triangle = 3
        n_edge_points = 3
        for level in range(n_triangulation_divisions):
            n_edge_points += n_edges_of_triangle * 2**level
            
        from eustace.analysis.advanced_standard.stats.spde.spherical_view import SphereMeshViewSuperTriangle
        n_points_in_region = SphereMeshViewSuperTriangle(n_triangulation_divisions, neighbourhood_level, 0 ).n_latent_variables()
        interior_points_per_region = n_points_in_region - n_edge_points
        
        # count matching points and assert
        for region_index in range(n_regions):
            self.assertEqual( interior_points_per_region, numpy.sum( merged_hyperparameters == hyperparameter_store[region_index][0] ) )
            self.assertEqual( interior_points_per_region, numpy.sum( merged_hyperparameters == hyperparameter_store[region_index][1] ) )