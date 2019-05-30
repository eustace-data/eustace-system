"""Test space-time component."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory

class TestComponentStorage_InMemory(unittest.TestCase):
    
    def test_component_storage_in_memory(self):
      
      storage = ComponentStorage_InMemory('A', 'B')
      self.assertEqual('A', storage.element_read())
      self.assertEqual('B', storage.hyperparameters_read())
      
class TestSpatialComponentSolutionStorage_InMemory(unittest.TestCase):

    def test_init(self):

        solutionstorage=SpatialComponentSolutionStorage_InMemory()
        self.assertDictEqual(solutionstorage.state_at_time, { })
        self.assertDictEqual(solutionstorage.state_marginal_std_at_time, { })
        self.assertDictEqual(solutionstorage.state_sample_at_time, { })

    def test_partial_write_and_read(self):
        solutionstorage = SpatialComponentSolutionStorage_InMemory()
        
        input_arrays = [numpy.array([1, 2, 3, 4]), numpy.array([4, 3, 2]), numpy.array([1])]
        input_time_indices = [435, 86758, 96758]
        
        for array, time_index in zip(input_arrays, input_time_indices):
            solutionstorage.partial_state_write(array, time_index)
            solutionstorage.partial_state_marginal_std_write(array, time_index)
            solutionstorage.partial_state_sample_write(array, time_index)
	    
        for array, time_index in zip(input_arrays, input_time_indices):
            numpy.testing.assert_array_equal(solutionstorage.partial_state_read(time_index), array)
            numpy.testing.assert_array_equal(solutionstorage.partial_state_marginal_std_read(time_index), array)
            numpy.testing.assert_array_equal(solutionstorage.partial_state_sample_read(time_index), array)


class TestSpaceTimeComponentSolutionStorage_InMemory(unittest.TestCase):
    
    def test_init(self):
      solutionstorage = SpaceTimeComponentSolutionStorage_InMemory()
      self.assertEqual(None, solutionstorage.state)
      self.assertEqual(None, solutionstorage.state_marginal_std)
      self.assertEqual(None, solutionstorage.state_sample)
      self.assertEqual({}, solutionstorage.measurement)  
      
    def test_write_and_read(self):
        solutionstorage = SpaceTimeComponentSolutionStorage_InMemory()
        
        input_array = numpy.array([1, 2, 3, 4])
        solutionstorage.state_write(input_array)
        numpy.testing.assert_array_equal(input_array, solutionstorage.state_read())
        
        solutionstorage.state_marginal_std_write(input_array)
        numpy.testing.assert_array_equal(input_array, solutionstorage.state_marginal_std_read())

        solutionstorage.state_sample_write(input_array)
        numpy.testing.assert_array_equal(input_array, solutionstorage.state_sample_read())


    def test_measurement_write_and_read(self):
        solutionstorage = SpaceTimeComponentSolutionStorage_InMemory()
        
        input_arrays = [numpy.array([1, 2, 3, 4]), numpy.array([4, 3, 2]), numpy.array([1])]
        input_time_indices = [435, 86758, 96758]
            
        for array, time_index in zip(input_arrays, input_time_indices):
            solutionstorage.measurement_write(array, time_index)

        self.assertItemsEqual(solutionstorage.timeindices_read(), input_time_indices)
        for array, time_index in zip(input_arrays, input_time_indices):
            numpy.testing.assert_array_equal(solutionstorage.measurement_read(time_index), array)


