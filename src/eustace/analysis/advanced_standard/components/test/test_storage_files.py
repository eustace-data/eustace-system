"""Test space-time component."""

import unittest
import numpy
import tempfile
import scipy.sparse

from eustace.analysis.advanced_standard.components.storage_files import SpaceTimeComponentSolutionStorage_Files
from eustace.analysis.advanced_standard.components.storage_files import SpatialComponentSolutionStorage_Files
from eustace.analysis.advanced_standard.components.storage_files import DelayedSpatialComponentSolutionStorage_Files
from eustace.analysis.advanced_standard.components.storage_files import DelayedSpatialComponentSolutionStorageFlexible_Files
class TestSpaceTimeComponentSolutionStorage_Files(unittest.TestCase):
    
    def test_init(self):
      
        solutionstorage = SpaceTimeComponentSolutionStorage_Files()
        self.assertEqual(None, solutionstorage.statefilename_read)
        self.assertEqual(None, solutionstorage.marginal_std_filename_read)
        self.assertEqual(None, solutionstorage.sample_filename_read)
        self.assertEqual(None, solutionstorage.statefilename_write)
        self.assertEqual(None, solutionstorage.marginal_std_filename_write)
        self.assertEqual(None, solutionstorage.sample_filename_write)
        self.assertEqual(None, solutionstorage.measurement_time_index_write)
        self.assertEqual(None, solutionstorage.measurementfilename_write)
        self.assertEqual(None, solutionstorage.measurementfilelist_read)
      
    def test_state_write_and_read(self):
        solutionstorage = SpaceTimeComponentSolutionStorage_Files()
        self.assertRaises(ValueError, solutionstorage.state_write, 3)
        self.assertEqual(None, solutionstorage.state_read())
        self.assertRaises(ValueError, solutionstorage.state_marginal_std_write, 3)
        self.assertEqual(None, solutionstorage.state_marginal_std_read())
        self.assertRaises(ValueError, solutionstorage.state_sample_write, 3)
        self.assertEqual(None, solutionstorage.state_sample_read())

        pickle_file = tempfile.NamedTemporaryFile(suffix='.pickle')
       
        # Checking state
        input_array = numpy.array([1, 2, 3, 4, 5, 6])
        solutionstorage.statefilename_write = solutionstorage.statefilename_read = pickle_file.name
        solutionstorage.state_write(input_array)
        numpy.testing.assert_array_equal(input_array, solutionstorage.state_read())
  
        # Checking marginal variances
        input_array = numpy.array([[1, 2], [3, 4], [5, 6]])
        solutionstorage.marginal_std_filename_write = solutionstorage.marginal_std_filename_read = pickle_file.name
        solutionstorage.state_marginal_std_write(input_array)
        numpy.testing.assert_array_equal(input_array, solutionstorage.state_marginal_std_read())
  
        # Checking sample
        input_array = numpy.array([[1, 2, 3], [4, 5, 6]])
        solutionstorage.sample_filename_write = solutionstorage.sample_filename_read = pickle_file.name
        solutionstorage.state_sample_write(input_array)
        numpy.testing.assert_array_equal(input_array, solutionstorage.state_sample_read())

    def test_measurement_write_and_read(self):

        solutionstorage = SpaceTimeComponentSolutionStorage_Files()
        self.assertRaises(ValueError, solutionstorage.measurement_write, 1, 2)
        solutionstorage.measurementfilename_write = 'A'
        self.assertRaises(ValueError, solutionstorage.measurement_write, 1, 2)
        self.assertRaises(ValueError, solutionstorage.measurement_read, 1)
        solutionstorage.measurementfilelist_read = [1]
        self.assertRaises(ValueError, solutionstorage.measurement_read, -1)
        self.assertRaises(ValueError, solutionstorage.measurement_read, 3)

        names = ['A', 'B', 'C']
        pickle_files = [tempfile.NamedTemporaryFile(suffix=name+'.pickle') for name in names]
        list_of_inputs = [numpy.array([1, 2, 3, 4]), numpy.array([3, 4]), numpy.array([1])] 

        solutionstorage.measurementfilelist_read = [pickle_file.name for pickle_file in pickle_files]
        self.assertEqual(range(len(names)), solutionstorage.timeindices_read())

        for time_index, input_array in enumerate(list_of_inputs):
            solutionstorage.measurement_time_index_write = time_index
            solutionstorage.measurementfilename_write = solutionstorage.measurementfilelist_read[time_index]
            solutionstorage.measurement_write(input_array, time_index)

            numpy.testing.assert_array_equal(solutionstorage.measurement_read(time_index), input_array)

class TestSpatialComponentSolutionStorage_Files(unittest.TestCase):
    
    def test_init(self):
      
        solutionstorage = SpatialComponentSolutionStorage_Files()
        self.assertEqual(None, solutionstorage.time_index)
        self.assertEqual(None, solutionstorage.statefilename_read)
        self.assertEqual(None, solutionstorage.marginal_std_filename_read)
        self.assertEqual(None, solutionstorage.sample_filename_read)
        self.assertEqual(None, solutionstorage.statefilename_write)
        self.assertEqual(None, solutionstorage.marginal_std_filename_write)
        self.assertEqual(None, solutionstorage.sample_filename_write)
      
    def test_state_write_and_read(self):
        solutionstorage = SpatialComponentSolutionStorage_Files()
        self.assertRaises(ValueError, solutionstorage.partial_state_write, 1, 2)
        solutionstorage.statefilename_write = 'A'
        self.assertRaises(ValueError, solutionstorage.partial_state_write, 1, 2)
        self.assertEqual(None, solutionstorage.partial_state_read(1))
        solutionstorage.time_index = 3
        solutionstorage.statefilename_read = 'A'
        self.assertRaises(ValueError, solutionstorage.partial_state_read, 1)

        self.assertRaises(ValueError, solutionstorage.partial_state_marginal_std_write, 1, 2)
        solutionstorage.marginal_std_filename_write = 'A'
        self.assertRaises(ValueError, solutionstorage.partial_state_marginal_std_write, 1, 2)
        self.assertEqual(None, solutionstorage.partial_state_marginal_std_read(1))
        solutionstorage.time_index = 3
        solutionstorage.marginal_std_filename_read = 'A'
        self.assertRaises(ValueError, solutionstorage.partial_state_marginal_std_read, 1)

        self.assertRaises(ValueError, solutionstorage.partial_state_sample_write, 1, 2)
        solutionstorage.sample_filename_write = 'A'
        self.assertRaises(ValueError, solutionstorage.partial_state_sample_write, 1, 2)
        self.assertEqual(None, solutionstorage.partial_state_sample_read(1))
        solutionstorage.time_index = 3
        solutionstorage.sample_filename_read = 'A'
        self.assertRaises(ValueError, solutionstorage.partial_state_sample_read, 1)

        names = ['A', 'B', 'C']
        pickle_files = [tempfile.NamedTemporaryFile(suffix=name+'.pickle') for name in names]
        list_of_inputs = [numpy.array([1, 2, 3, 4]), numpy.array([3, 4]), numpy.array([1])] 
      
        # Checking state
        for time_index, input_array in enumerate(list_of_inputs):
            solutionstorage.time_index = time_index
            solutionstorage.statefilename_write = solutionstorage.statefilename_read = pickle_files[time_index].name
            solutionstorage.partial_state_write(input_array, time_index)

            numpy.testing.assert_array_equal(solutionstorage.partial_state_read(time_index), input_array)

        list_of_inputs = [numpy.array([[1, 2], [3, 4]]), numpy.array([3, 4,3]), numpy.array([10])] 
        # Checking marginal variances
        for time_index, input_array in enumerate(list_of_inputs):
            solutionstorage.time_index = time_index
            solutionstorage.marginal_std_filename_write = solutionstorage.marginal_std_filename_read = pickle_files[time_index].name
            solutionstorage.partial_state_marginal_std_write(input_array, time_index)
            
            numpy.testing.assert_array_equal(solutionstorage.partial_state_marginal_std_read(time_index), input_array)
    
        list_of_inputs = [numpy.array([[1, 2, 3, 4]]), numpy.array([[3, 3], [4,3]]), numpy.array([-10])] 
        # Checking marginal variances
        for time_index, input_array in enumerate(list_of_inputs):
            solutionstorage.time_index = time_index
            solutionstorage.sample_filename_write = solutionstorage.sample_filename_read = pickle_files[time_index].name
            solutionstorage.partial_state_sample_write(input_array, time_index)
            
            numpy.testing.assert_array_equal(solutionstorage.partial_state_sample_read(time_index), input_array)


class TestDelayedSpatialComponentSolutionStorage_Files(unittest.TestCase):
    
    def test_init(self):
      
        solutionstorage = DelayedSpatialComponentSolutionStorage_Files()
        
        self.assertEqual({}, solutionstorage.statefiledictionary_read)
        self.assertEqual({}, solutionstorage.marginal_std_filedictionary_read)
        self.assertEqual({}, solutionstorage.sample_filedictionary_read)
        
        self.assertEqual({}, solutionstorage.statefiledictionary_write)
        self.assertEqual({}, solutionstorage.marginal_std_filedictionary_write)
        self.assertEqual({}, solutionstorage.sample_filedictionary_write)
        
        self.assertEqual({}, solutionstorage.measurementfiledictionary_write)
        self.assertEqual({}, solutionstorage.measurementfiledictionary_read)
      
    def test_state_write_and_read(self):
        solutionstorage = DelayedSpatialComponentSolutionStorage_Files()
        
        self.assertRaises(ValueError, solutionstorage.state_write, 1, 2)
        solutionstorage.statefiledictionary_write = {0: 'A'}
        self.assertRaises(ValueError, solutionstorage.state_write, 1, 2)
        self.assertEqual(None, solutionstorage.state_read(1))

        self.assertRaises(ValueError, solutionstorage.state_marginal_std_write, 1, 2)
        solutionstorage.marginal_std_filedictionary_write = {0: 'A'}
        self.assertRaises(ValueError, solutionstorage.state_marginal_std_write, 1, 2)
        self.assertEqual(None, solutionstorage.state_marginal_std_read(1))
        
        names = ['A', 'B', 'C']
        pickle_files = [tempfile.NamedTemporaryFile(suffix=name+'.pickle') for name in names]
        list_of_inputs = [numpy.array([1, 2, 3, 4]), numpy.array([3, 4]), numpy.array([1])] 
      
        # Checking state
        for state_time_index, input_array in enumerate(list_of_inputs):
            solutionstorage.statefiledictionary_write[state_time_index] = pickle_files[state_time_index].name
            solutionstorage.statefiledictionary_read[state_time_index] = pickle_files[state_time_index].name
            solutionstorage.state_write(input_array, state_time_index)
    
            numpy.testing.assert_array_equal(solutionstorage.state_read(state_time_index), input_array)

        # Checking marginal variances
        for state_time_index, input_array in enumerate(list_of_inputs):
            solutionstorage.marginal_std_filedictionary_write[state_time_index] = pickle_files[state_time_index].name
            solutionstorage.marginal_std_filedictionary_read[state_time_index] = pickle_files[state_time_index].name
            solutionstorage.state_marginal_std_write(input_array, state_time_index)

            numpy.testing.assert_array_equal(solutionstorage.state_marginal_std_read(state_time_index), input_array)

    def test_measurement_write_and_read(self):

        solutionstorage = DelayedSpatialComponentSolutionStorage_Files()
        self.assertRaises(ValueError, solutionstorage.measurement_write, 1, 2)
        solutionstorage.measurementfilename_write = 'A'
        self.assertRaises(ValueError, solutionstorage.measurement_write, 1, 2)
        self.assertRaises(ValueError, solutionstorage.measurement_read, 1)
        solutionstorage.measurementfilelist_read = [1]
        self.assertRaises(ValueError, solutionstorage.measurement_read, -1)
        self.assertRaises(ValueError, solutionstorage.measurement_read, 3)

        names = ['A', 'B', 'C']
        pickle_files = [tempfile.NamedTemporaryFile(suffix=name+'.pickle') for name in names]
        list_of_inputs = [numpy.array([1, 2, 3, 4]), numpy.array([3, 4]), numpy.array([1])] 

        solutionstorage.measurementfiledictionary_read = {time_index: pickle_file.name for time_index, pickle_file in enumerate(pickle_files)}
        self.assertEqual(range(len(names)), solutionstorage.measurementfiledictionary_read.keys())
      
        for time_index, input_array in enumerate(list_of_inputs):
            
            solutionstorage.measurementfiledictionary_write[time_index] = solutionstorage.measurementfiledictionary_read[time_index]
            solutionstorage.measurement_write(input_array, time_index)

            numpy.testing.assert_array_equal(solutionstorage.measurement_read(time_index), input_array)
            
class TestDelayedSpatialComponentSolutionStorageFlexible_Files(unittest.TestCase):
    
    def test_init(self):
        
        solutionstorage = DelayedSpatialComponentSolutionStorageFlexible_Files( state_read_pattern="state_read_pattern.%Y.%m.%d.txt",
                                                                                marginal_std_read_pattern = "marginal_std_read_pattern.%Y.%m.%d.txt",
                                                                                sample_read_pattern="sample_read_pattern.%Y.%m.%d.txt",
                                                                                measurement_read_pattern="measurement_read_pattern.%Y.%m.%d.txt",
                                                                                state_write_pattern="state_write_pattern.%Y.%m.%d.txt",
                                                                                marginal_std_write_pattern = "marginal_std_write_pattern.%Y.%m.%d.txt",
                                                                                sample_write_pattern="sample_write_pattern.%Y.%m.%d.txt",
                                                                                measurement_write_pattern="measurement_write_pattern.%Y.%m.%d.txt")
        
        self.assertEqual("state_read_pattern.%Y.%m.%d.txt", solutionstorage.state_read_pattern)
        self.assertEqual("marginal_std_read_pattern.%Y.%m.%d.txt", solutionstorage.marginal_std_read_pattern)
        self.assertEqual("sample_read_pattern.%Y.%m.%d.txt", solutionstorage.sample_read_pattern)
        
        self.assertEqual("state_write_pattern.%Y.%m.%d.txt", solutionstorage.state_write_pattern)
        self.assertEqual("marginal_std_write_pattern.%Y.%m.%d.txt", solutionstorage.marginal_std_write_pattern)
        self.assertEqual("sample_write_pattern.%Y.%m.%d.txt", solutionstorage.sample_write_pattern)
        
        self.assertEqual("measurement_write_pattern.%Y.%m.%d.txt", solutionstorage.measurement_write_pattern)
        self.assertEqual("measurement_read_pattern.%Y.%m.%d.txt", solutionstorage.measurement_read_pattern)
    
    def test_filename_expansion(self):
        
        solutionstorage = DelayedSpatialComponentSolutionStorageFlexible_Files( state_read_pattern=["%Y", "state_read_pattern.%Y.%m.%d.txt"],
                                                                                marginal_std_read_pattern = ["%Y", "marginal_std_read_pattern.%Y.%m.%d.txt"],
                                                                                sample_read_pattern=["%Y", "sample_read_pattern.%Y.%m.%d.txt"],
                                                                                measurement_read_pattern=["%Y", "measurement_read_pattern.%Y.%m.%d.txt"],
                                                                                state_write_pattern=["%Y", "state_write_pattern.%Y.%m.%d.txt"],
                                                                                marginal_std_write_pattern = ["%Y", "marginal_std_write_pattern.%Y.%m.%d.txt"],
                                                                                sample_write_pattern=["%Y", "sample_write_pattern.%Y.%m.%d.txt"],
                                                                                measurement_write_pattern=["%Y", "measurement_write_pattern.%Y.%m.%d.txt"])
        
        self.assertEqual("1850/state_write_pattern.1850.01.01.txt", solutionstorage.filename_from_patterns(solutionstorage.state_write_pattern, 0))
        self.assertEqual("1850/state_write_pattern.1850.01.02.txt", solutionstorage.filename_from_patterns(solutionstorage.state_write_pattern, 1))
        self.assertEqual("1850/state_write_pattern.1850.01.03.txt", solutionstorage.filename_from_patterns(solutionstorage.state_write_pattern, 2))
        self.assertEqual("1850/state_write_pattern.1850.02.01.txt", solutionstorage.filename_from_patterns(solutionstorage.state_write_pattern, 31))
        self.assertEqual("1851/state_write_pattern.1851.01.01.txt", solutionstorage.filename_from_patterns(solutionstorage.state_write_pattern, 365))
    
        list_of_inputs = [numpy.array([1, 2, 3, 4]), numpy.array([3, 4]), numpy.array([1])] 
    
    def test_state_write_and_read(self):
        
        test_directory = tempfile.mkdtemp()

        solutionstorage = DelayedSpatialComponentSolutionStorageFlexible_Files( state_read_pattern=[test_directory, "%Y", "state.%Y.%m.%d.txt"],
                                                                                state_write_pattern=[test_directory, "%Y", "state.%Y.%m.%d.txt"])
        
        list_of_inputs = [numpy.array([1, 2, 3, 4]), numpy.array([3, 4]), numpy.array([1])] 
        for time_index, input_array in enumerate(list_of_inputs):
            
            solutionstorage.state_write(input_array, time_index)
            numpy.testing.assert_array_equal(solutionstorage.state_read(time_index), input_array)
            
    def test_state_marginal_std_write_and_read(self):
        
        test_directory = tempfile.mkdtemp()

        solutionstorage = DelayedSpatialComponentSolutionStorageFlexible_Files( marginal_std_read_pattern = [test_directory, "%Y", "marginal_std.%Y.%m.%d.txt"],
                                                                                marginal_std_write_pattern = [test_directory, "%Y", "marginal_std.%Y.%m.%d.txt"])
        
        list_of_inputs = [numpy.array([2, 3, 4, 5]), numpy.array([4, 5]), numpy.array([2])] 
        for time_index, input_array in enumerate(list_of_inputs):
            
            solutionstorage.state_marginal_std_write(input_array, time_index)
            numpy.testing.assert_array_equal(solutionstorage.state_marginal_std_read(time_index), input_array)
            
    def test_state_sample_write_and_read(self):
        
        test_directory = tempfile.mkdtemp()

        solutionstorage = DelayedSpatialComponentSolutionStorageFlexible_Files( sample_read_pattern=[test_directory, "%Y", "sample.%Y.%m.%d.txt"],
                                                                                sample_write_pattern=[test_directory, "%Y", "sample.%Y.%m.%d.txt"])
        
        list_of_inputs = [numpy.array([3, 4, 5, 6]), numpy.array([5, 6]), numpy.array([3])] 
        for time_index, input_array in enumerate(list_of_inputs):
            
            solutionstorage.state_sample_write(input_array, time_index)
            numpy.testing.assert_array_equal(solutionstorage.state_sample_read(time_index), input_array)
    
    def test_measurement_write_and_read(self):
        
        test_directory = tempfile.mkdtemp()

        solutionstorage = DelayedSpatialComponentSolutionStorageFlexible_Files( measurement_read_pattern=[test_directory, "%Y", "measurement.%Y.%m.%d.txt"],
                                                                                measurement_write_pattern=[test_directory, "%Y", "measurement.%Y.%m.%d.txt"])
        
        list_of_inputs = [numpy.array([4, 5, 6, 7]), numpy.array([6, 7]), numpy.array([4])] 
        for time_index, input_array in enumerate(list_of_inputs):
            
            solutionstorage.measurement_write(input_array, time_index)
            numpy.testing.assert_array_equal(solutionstorage.measurement_read(time_index), input_array)