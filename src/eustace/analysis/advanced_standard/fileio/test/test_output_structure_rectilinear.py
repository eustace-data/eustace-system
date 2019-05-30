"""Testing output rectilinear structure and regridder, for several different cases"""

import numpy
import tempfile
import unittest

from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import RegriddingRectilinearGridStructure
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import Regridder

from eustace.outputformats.definitions import GLOBAL_FIELD_OUTPUT_FLAGS

class TestOutputRectilinearGridStructure(unittest.TestCase):
  
    def test_init(self):
      A = OutputRectilinearGridStructure('A', 'B', 'C', 'D')
      
      self.assertEqual('A', A.time_index_number)
      self.assertEqual('B', A.corresponding_datetime)
      self.assertEqual('C', A.latitudes)
      self.assertEqual('D', A.longitudes)
      
    def test_number_of_observations(self):
      A = OutputRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3]), numpy.array([.1, .3]))
      
      self.assertEqual(6, A.number_of_observations())
      
    def test_location_polar_coordinates(self):
      A = OutputRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3]), numpy.array([.1, .3]))
      
      expected_array = numpy.array([[1, .1],[1, .3], [2, .1], [2, .3], [3 ,.1], [3, .3]])
      numpy.testing.assert_array_equal(expected_array, A.location_polar_coordinates())
      
      A = OutputRectilinearGridStructure('A', 'B', numpy.array([1, 2]), numpy.array([.1, .3, .2]))
      
      expected_array = numpy.array([[1, .1],[1, .3], [1, .2], [2, .1], [2 ,.3], [2, .2]])
      numpy.testing.assert_array_equal(expected_array, A.location_polar_coordinates())
      
class TestRegriddingRectilinearGridStructure(unittest.TestCase):
  
    def test_init(self):
      A = RegriddingRectilinearGridStructure('A', 'B', 'C', 'D', 'E', 'F')
      
      self.assertEqual('A', A.time_index_number)
      self.assertEqual('B', A.corresponding_datetime)
      self.assertEqual('C', A.latitudes)
      self.assertEqual('D', A.longitudes)
      self.assertEqual('E', A.Nlat)
      self.assertEqual('F', A.Nlon)
      
    def test_number_of_observations(self):
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3]), numpy.array([.1, .3]),2 ,3 )
      
      self.assertEqual(6, A.number_of_observations())
      
    def test_rearrange_latitude_and_longitudes(self):
      # First case: 6 grid cells, 2 latitude points for each cell
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3, 4, 5, 6]), numpy.array([.1, .2, .3]),3 ,1)
      
      expected_latitudes = numpy.array([[1, 2, 3], [1, 2, 3], [1, 2, 3], [4, 5, 6], [4, 5, 6], [4, 5, 6]])
      expected_longitudes = numpy.array([[0.1, 0.1, 0.1], [0.2, 0.2, 0.2], [0.3, 0.3, 0.3], [0.1, 0.1, 0.1], [0.2, 0.2, 0.2], [0.3, 0.3, 0.3]])
      numpy.testing.assert_array_equal(expected_latitudes, A.rearrange_latitudes())
      numpy.testing.assert_array_equal(expected_longitudes, A.rearrange_longitudes())

      # Second case: 9 grid cells, 3 latitude points for each cell
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3, 4, 5, 6]), numpy.array([.1, .2, .3]),2 ,1)
      
      expected_latitudes = numpy.array([[1, 2], [1, 2], [1, 2], [3, 4], [3, 4], [3, 4], [5, 6], [5, 6], [5, 6]])
      expected_longitudes = numpy.array([[0.1, 0.1], [0.2, 0.2], [0.3, 0.3], [0.1, 0.1], [0.2, 0.2], [0.3, 0.3], [0.1, 0.1], [0.2, 0.2], [0.3, 0.3]])
      numpy.testing.assert_array_equal(expected_latitudes, A.rearrange_latitudes())
      numpy.testing.assert_array_equal(expected_longitudes, A.rearrange_longitudes())

      #3d case: 3 grid cells, 2 latitude, 2 longitude points for each cell
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3, 4, 5, 6]), numpy.array([.1, .2]),2 ,2)
      
      expected_latitudes = numpy.array([[1, 2, 1, 2], [3, 4, 3, 4], [5, 6, 5, 6]])
      expected_longitudes = numpy.array([[.1, .1, .2, .2], [.1, .1, .2, .2], [.1, .1, .2, .2]])
      numpy.testing.assert_array_equal(expected_latitudes, A.rearrange_latitudes())
      numpy.testing.assert_array_equal(expected_longitudes, A.rearrange_longitudes())
      
      #4th case: 2 grid cells, 6 latitude points, 2 points for each cell
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3, 4, 5, 6]), numpy.array([.1, .2, .3, .4]),6 ,2)
      
      expected_latitudes = numpy.array([[1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6], [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6]])
      expected_longitudes = numpy.array([[.1, .1, .1, .1, .1, .1, .2, .2, .2, .2, .2, .2], [.3, .3, .3, .3, .3, .3, .4, .4, .4, .4, .4, .4]])
      numpy.testing.assert_array_equal(expected_latitudes, A.rearrange_latitudes())
      
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3, 4, 5, 6]), numpy.array([.1, .2, .3]),1 ,1)
      
      expected_latitudes = numpy.array([[1], [1], [1], [2], [2], [2], [3], [3], [3], [4], [4], [4], [5], [5], [5], [6], [6], [6]])
      numpy.testing.assert_array_equal(expected_latitudes, A.rearrange_latitudes())
      
    def test_location_polar_coordinates(self):
      # First case: 6 grid cells, 2 latitude points for each cell
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3, 4, 5, 6]), numpy.array([.1, .2, .3]),3 ,1)

      expected_coordinates = numpy.array([[1, 0.1], [2, 0.1], [3, 0.1], 
				          [1, 0.2], [2, 0.2], [3, 0.2], 
				          [1, 0.3], [2, 0.3], [3, 0.3], 
				          [4, 0.1], [5, 0.1], [6, 0.1],
				          [4, 0.2], [5, 0.2], [6, 0.2],
				          [4, 0.3], [5, 0.3], [6, 0.3]])
      numpy.testing.assert_array_equal(expected_coordinates, A.location_polar_coordinates())
     
      #3d case: 3 grid cells, 2 latitude, 2 longitude points for each cell
      A = RegriddingRectilinearGridStructure('A', 'B', numpy.array([1, 2, 3, 4, 5, 6]), numpy.array([.1, .2]),2 ,2)
      
      expected_coordinates = numpy.array([[1, 0.1], [2, 0.1],
					  [1, 0.2], [2, 0.2],
					  [3, 0.1], [4, 0.1],
					  [3, 0.2], [4, 0.2],
					  [5, 0.1], [6, 0.1],
					  [5, 0.2], [6, 0.2]])
					  
class TestRegridder(unittest.TestCase):
  
    class DummyComponentSolution(object):
	
	def __init__(self):
	  pass
	
	def solution_observation_expected_value(self, observations):
	    return numpy.zeros([observations.number_of_observations(),1])+1.
		    
	def solution_observation_expected_uncertainties(self, observations):
	    return numpy.zeros([observations.number_of_observations(),1])+.5
	    
	def solution_observation_prior_uncertainties(self, observations):
	    return numpy.zeros([observations.number_of_observations(),1])+.7
	    
    class LessDummyComponentSolution(object):
	
	def __init__(self):
	  pass
	
	def solution_observation_expected_value(self, observations):
	    vector = numpy.zeros([observations.number_of_observations(),1])
	    vector[0] = 1.
	    vector[1] = .5
	    vector[-1] = -.5
	    return vector

	def solution_observation_expected_uncertainties(self, observations):
	    vector = numpy.zeros([observations.number_of_observations(),1])
	    vector[2] = 1.
	    vector[3] = .5
	    vector[-2] = -.5
	    return vector

	def solution_observation_prior_uncertainties(self, observations):
	    vector = numpy.zeros([observations.number_of_observations(),1])
	    vector[1] = 1.
	    vector[2] = .5
	    vector[-1] = -.5
	    return vector

    class HarmonicComponentSolution(object):
	
	def __init__(self):
	  pass
	
	def solution_observation_expected_value(self, observations):
	    harmonic_vector = numpy.cos(numpy.radians(observations.location_polar_coordinates()[:,0]))
	    return harmonic_vector

	def solution_observation_expected_uncertainties(self, observations):
	    harmonic_vector = numpy.sin(numpy.radians(observations.location_polar_coordinates()[:,0]))
	    return harmonic_vector

	def solution_observation_prior_uncertainties(self, observations):
	    harmonic_vector = numpy.square(numpy.sin(numpy.radians(observations.location_polar_coordinates()[:,0])))
	    return harmonic_vector

    def setUp(self):
      """Set up mock objects to be used for testing the Regridder class functionalities"""
      
      self.structure = OutputRectilinearGridStructure(1, None, numpy.array([1, 2, 3, 4]), numpy.array([.2, .3, .4]))
      
      # Example of latitude and longitude values for different sets of points
      self.points = [numpy.array([[-180.],[-90.],[0.]]), 
		     numpy.array([[-180., -90.], [0., 90.], [-45, 45]]),
		     numpy.array([[-180., -90., 0., 90., -45, 45]])]
		     
      # Expected harmonic factors
      self.expected_arrays = [numpy.array([[-1],[ 0.],[ 1.]]), 
		              numpy.array([[-1, 0.],[ 1., 0.], [numpy.sqrt(2.)/2., numpy.sqrt(2.)/2.]]),
		              numpy.array([[-1, 0., 1., 0., numpy.sqrt(2.)/2., numpy.sqrt(2.)/2.]])]

      self.expected_normalizations = [numpy.array([[-1], [0.], [ 1.]]),
                                      numpy.array([[-1], [1.], [numpy.sqrt(2.)]]),
                                      numpy.array([[numpy.sqrt(2.)]])]

      # We cannot divide by zero, we discard the first normalization factor
      self.expected_weithing_factors = []
      for index in range(1, len(self.expected_arrays)):
	self.expected_weithing_factors.append(self.expected_arrays[index]/self.expected_normalizations[index])
      
      self.cell_dimensions = [[1,1], [2,1], [6,1]]
    
    def test_init(self):
      A = Regridder(self.structure, [1,1])
      
      self.assertEqual(1, A.time_index_number)
      self.assertEqual(None, A.corresponding_datetime)
      numpy.testing.assert_array_equal(numpy.array([1, 2, 3, 4]), A.latitudes)
      numpy.testing.assert_array_equal(numpy.array([.2, .3, .4]), A.longitudes)
      self.assertEqual(1, A.Nlat)
      self.assertEqual(1, A.Nlon)
      self.assertEqual(100, A.blocking)
      self.assertEqual(1, A.latitude_spacing)
      self.assertEqual(.1, A.longitude_spacing)  
      self.assertEqual(12, A.number_of_observations)
      self.assertListEqual(GLOBAL_FIELD_OUTPUT_FLAGS, A.field_flags)
     
    def test_harmonic_factor(self):
      
      for cell_dimension, points, expected_array in zip(self.cell_dimensions, self.points, self.expected_arrays):
	  A = Regridder(self.structure, cell_dimension)
	  numpy.testing.assert_array_almost_equal(expected_array, A.harmonic_factor(points))
      points = numpy.array([[-30., -90.], [0., 90.], [-45, 45]])
      self.assertTrue((A.harmonic_factor(points)>=0.).all())

    def test_normalization_factor(self):
      
      for cell_dimension, points, expected_normalization in zip(self.cell_dimensions, self.points, self.expected_normalizations):
	  A = Regridder(self.structure, cell_dimension)
	  numpy.testing.assert_array_almost_equal(expected_normalization, A.normalization_factor(points))
      
    def test_weighting_factors(self):
      
      for index in range(len(self.expected_weithing_factors)):
	  shifted_index = index+1
	  A = Regridder(self.structure,  self.cell_dimensions[shifted_index])
  	  numpy.testing.assert_array_almost_equal(self.expected_weithing_factors[index], A.weighting_factors(self.points[shifted_index]))
	  numpy.testing.assert_array_almost_equal(1., A.weighting_factors(self.points[shifted_index]).sum(axis=1))
	  
    def test_create_sub_cell_grids(self):

      input_points = [self.structure.latitudes, self.structure.longitudes]

      A = Regridder(self.structure, [1,1])
      structure = A.create_sub_cell_grids(input_points)
      self.assertTrue(isinstance(structure, RegriddingRectilinearGridStructure))
      numpy.testing.assert_array_equal(input_points[0], structure.latitudes)
      numpy.testing.assert_array_equal(input_points[1], structure.longitudes)
      
      A = Regridder(self.structure, [2,1])
      structure = A.create_sub_cell_grids(input_points)
      expected_points = numpy.array([[0.75, 0.2], [1.25, 0.2], 
				     [0.75, 0.3], [1.25, 0.3],
				     [0.75, 0.4], [1.25, 0.4],
				     [1.75, 0.2], [2.25, 0.2], 
				     [1.75, 0.3], [2.25, 0.3],
				     [1.75, 0.4], [2.25, 0.4],
				     [2.75, 0.2], [3.25, 0.2], 
				     [2.75, 0.3], [3.25, 0.3],
				     [2.75, 0.4], [3.25, 0.4],
				     [3.75, 0.2], [4.25, 0.2], 
				     [3.75, 0.3], [4.25, 0.3],
				     [3.75, 0.4], [4.25, 0.4]])
      numpy.testing.assert_array_equal(expected_points, structure.location_polar_coordinates())
      
      A = Regridder(self.structure, [1,2])
      structure = A.create_sub_cell_grids(input_points)
      expected_points = numpy.array([[1, 0.175], [1, 0.225], 
                                     [1, 0.275], [1, 0.325], 
                                     [1, 0.375], [1, 0.425],
                                     [2, 0.175], [2, 0.225], 
                                     [2, 0.275], [2, 0.325], 
                                     [2, 0.375], [2, 0.425],
                                     [3, 0.175], [3, 0.225], 
                                     [3, 0.275], [3, 0.325], 
                                     [3, 0.375], [3, 0.425],
                                     [4, 0.175], [4, 0.225], 
                                     [4, 0.275], [4, 0.325], 
                                     [4, 0.375], [4, 0.425]])
      numpy.testing.assert_array_almost_equal(expected_points, structure.location_polar_coordinates())
      
      A = Regridder(self.structure, [3,2])
      structure = A.create_sub_cell_grids(input_points)
      expected_points = numpy.array([[2./3., 0.175], [1., 0.175], [4./3., 0.175], [2./3., 0.225], [1., 0.225], [4./3., 0.225],
                                     [2./3., 0.275], [1., 0.275], [4./3., 0.275], [2./3., 0.325], [1., 0.325], [4./3., 0.325],
                                     [2./3., 0.375], [1., 0.375], [4./3., 0.375], [2./3., 0.425], [1., 0.425], [4./3., 0.425],
				     [5./3., 0.175], [2., 0.175], [7./3., 0.175], [5./3., 0.225], [2., 0.225], [7./3., 0.225],
				     [5./3., 0.275], [2., 0.275], [7./3., 0.275], [5./3., 0.325], [2., 0.325], [7./3., 0.325],
				     [5./3., 0.375], [2., 0.375], [7./3., 0.375], [5./3., 0.425], [2., 0.425], [7./3., 0.425],
				     [8./3., 0.175], [3., 0.175], [10./3., 0.175], [8./3., 0.225], [3., 0.225], [10./3., 0.225],
				     [8./3., 0.275], [3., 0.275], [10./3., 0.275], [8./3., 0.325], [3., 0.325], [10./3., 0.325],
				     [8./3., 0.375], [3., 0.375], [10./3., 0.375], [8./3., 0.425], [3., 0.425], [10./3., 0.425],
				     [11./3., 0.175], [4., 0.175], [13./3., 0.175], [11./3., 0.225], [4., 0.225], [13./3., 0.225],
				     [11./3., 0.275], [4., 0.275], [13./3., 0.275], [11./3., 0.325], [4., 0.325], [13./3., 0.325],
				     [11./3., 0.375], [4., 0.375], [13./3., 0.375], [11./3., 0.425], [4., 0.425], [13./3., 0.425]])
      numpy.testing.assert_array_almost_equal(expected_points, structure.location_polar_coordinates())
    
    def test_compute_blocked_gridded_expected_value_dummy_projection(self):

	blocked_latitudes, longitudes = self.structure.latitudes, self.structure.longitudes

	# Different integral estimation for a trivial case

	component_solution = TestRegridder.DummyComponentSolution()

	for field, value in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [1., .5, .7]):
	    for shapes in [[1,1], [1,2], [2,1]]:
		A = Regridder(self.structure, shapes)
		numpy.testing.assert_array_equal(value, A.compute_blocked_gridded_expected_value(field, blocked_latitudes, longitudes, component_solution))
	    	
	# Less trivial case: non trivial projection
	component_solution = TestRegridder.LessDummyComponentSolution()

	# Point like estimation for less trivial case
	A = Regridder(self.structure, [1,1])
	numpy.testing.assert_array_equal(component_solution.solution_observation_expected_value(self.structure).ravel(), A.compute_blocked_gridded_expected_value('MAP', blocked_latitudes, longitudes, component_solution))
	numpy.testing.assert_array_equal(component_solution.solution_observation_expected_uncertainties(self.structure).ravel(), A.compute_blocked_gridded_expected_value('post_STD', blocked_latitudes, longitudes, component_solution))
	numpy.testing.assert_array_equal(component_solution.solution_observation_prior_uncertainties(self.structure).ravel(), A.compute_blocked_gridded_expected_value('prior_STD', blocked_latitudes, longitudes, component_solution))
	
	#1 latitude point, 2 longitude points per grid cell
	A = Regridder(self.structure, [1,2])

	expected_MAP_result = numpy.array([1.5/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -.25])
	expected_post_STD_result = numpy.array([0., 1.5/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., -.25])
	expected_prior_STD_result = numpy.array([0.5, .25, 0., 0., 0., 0., 0., 0., 0., 0., 0., -.25])
	
	for field, array in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [expected_MAP_result, expected_post_STD_result, expected_prior_STD_result]):
	    numpy.testing.assert_array_equal(array, A.compute_blocked_gridded_expected_value(field, blocked_latitudes, longitudes, component_solution))

	#2 latitude points, 3 longitude points per grid cell
	A = Regridder(self.structure, [2,3])
	expected_weight_latitudes = numpy.array([[0.75, 1.25, 0.75, 1.25, 0.75, 1.25], [0.75, 1.25, 0.75, 1.25, 0.75, 1.25],[0.75, 1.25, 0.75, 1.25, 0.75, 1.25],
					         [1.75, 2.25, 1.75, 2.25, 1.75, 2.25], [1.75, 2.25, 1.75, 2.25, 1.75, 2.25],[1.75, 2.25, 1.75, 2.25, 1.75, 2.25],
					         [2.75, 3.25, 2.75, 3.25, 2.75, 3.25], [2.75, 3.25, 2.75, 3.25, 2.75, 3.25],[2.75, 3.25, 2.75, 3.25, 2.75, 3.25],
					         [3.75, 4.25, 3.75, 4.25, 3.75, 4.25], [3.75, 4.25, 3.75, 4.25, 3.75, 4.25],[3.75, 4.25, 3.75, 4.25, 3.75, 4.25]])
					                                        
	expected_MAP_array = numpy.array([[1., .5, 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					  [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					  [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])
	expected_post_STD_array = numpy.array([[0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					       [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					       [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., -.5, 0.]])
	expected_prior_STD_array = numpy.array([[0., 1., .5, 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					        [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					        [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])


	expected_MAP_result = (A.weighting_factors(expected_weight_latitudes)*expected_MAP_array).sum(axis=1)
	expected_post_STD_result = (A.weighting_factors(expected_weight_latitudes)*expected_post_STD_array).sum(axis=1)
	expected_prior_STD_result = (A.weighting_factors(expected_weight_latitudes)*expected_prior_STD_array).sum(axis=1)
	
	for field, array in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [expected_MAP_result, expected_post_STD_result, expected_prior_STD_result]):
	    numpy.testing.assert_array_equal(array, A.compute_blocked_gridded_expected_value(field, blocked_latitudes, longitudes, component_solution))

    def test_compute_gridded_expected_value_dummy_projection(self):
      
	# Point like estimation for a trivial case
	
	component_solution = TestRegridder.DummyComponentSolution()
	
	for field, value in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [1., .5, .7]):
	    for blocking in [1, 2, 3, 6, 10]:
		A = Regridder(self.structure, [1,1], blocking)
		numpy.testing.assert_array_equal(value, A.compute_gridded_expected_value(field, component_solution))

	# Less trivial component solution
	
	component_solution = TestRegridder.LessDummyComponentSolution()

    	A = Regridder(self.structure, [1,2], blocking=4)
	#1 latitude point, 2 longitude points per grid cell	
	expected_MAP_result = numpy.array([1.5/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., -.25])
	expected_post_STD_result = numpy.array([0., 1.5/2., 0., 0., 0., 0., 0., 0., 0., 0., 0., -.25])
	expected_prior_STD_result = numpy.array([0.5, .25, 0., 0., 0., 0., 0., 0., 0., 0., 0., -.25])

	for field, array in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [expected_MAP_result, expected_post_STD_result, expected_prior_STD_result]):
	    numpy.testing.assert_array_equal(array, A.compute_gridded_expected_value(field, component_solution))
    
	A = Regridder(self.structure, [2,3], blocking=5)
	expected_weight_latitudes = numpy.array([[0.75, 1.25, 0.75, 1.25, 0.75, 1.25], [0.75, 1.25, 0.75, 1.25, 0.75, 1.25],[0.75, 1.25, 0.75, 1.25, 0.75, 1.25],
					         [1.75, 2.25, 1.75, 2.25, 1.75, 2.25], [1.75, 2.25, 1.75, 2.25, 1.75, 2.25],[1.75, 2.25, 1.75, 2.25, 1.75, 2.25],
					         [2.75, 3.25, 2.75, 3.25, 2.75, 3.25], [2.75, 3.25, 2.75, 3.25, 2.75, 3.25],[2.75, 3.25, 2.75, 3.25, 2.75, 3.25],
					         [3.75, 4.25, 3.75, 4.25, 3.75, 4.25], [3.75, 4.25, 3.75, 4.25, 3.75, 4.25],[3.75, 4.25, 3.75, 4.25, 3.75, 4.25]])
					                                        
	expected_MAP_array = numpy.array([[1., .5, 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					  [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					  [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])
	expected_post_STD_array = numpy.array([[0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
			    		       [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					       [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., -.5, 0.]])
	expected_prior_STD_array = numpy.array([[0., 1., .5, 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					        [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					        [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])

	expected_MAP_result = (A.weighting_factors(expected_weight_latitudes)*expected_MAP_array).sum(axis=1)
	expected_post_STD_result = (A.weighting_factors(expected_weight_latitudes)*expected_post_STD_array).sum(axis=1)
	expected_prior_STD_result = (A.weighting_factors(expected_weight_latitudes)*expected_prior_STD_array).sum(axis=1)

	for field, array in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [expected_MAP_result, expected_post_STD_result, expected_prior_STD_result]):
	    numpy.testing.assert_array_equal(array, A.compute_gridded_expected_value(field, component_solution))
    
    def test_compute_gridded_expected_value_harmonic_projection(self):
	#Compute cell grid area averages, with projection equalling the cosine of the latitude of each point.
	
	new_structure = OutputRectilinearGridStructure(1, None, numpy.array([-60., 0., 60.]), numpy.array([-90., -60., 60., 90.]))
	A = Regridder(new_structure, [2, 3], blocking=2)
	component_solution = TestRegridder.HarmonicComponentSolution()

	# Check it raises an exception if wrong field flags are given
	self.assertRaises(ValueError, A.compute_gridded_expected_value, 'MMP', component_solution)

	my_expected_grid = numpy.array([[-75., -45., -75., -45.,-75., -45.],[-75., -45., -75., -45.,-75., -45.], [-75., -45., -75., -45.,-75., -45.], [-75., -45., -75., -45.,-75., -45.],
					[-15., 15., -15., 15.,-15., 15.], [-15., 15., -15., 15.,-15., 15.], [-15., 15., -15., 15.,-15., 15.], [-15., 15., -15., 15.,-15., 15.],
					[45., 75., 45., 75., 45., 75.], [45., 75., 45., 75., 45., 75.], [45., 75., 45., 75., 45., 75.], [45., 75., 45., 75., 45., 75.]])

	expected_MAP_result = (A.weighting_factors(my_expected_grid)*numpy.cos(numpy.radians(my_expected_grid))).sum(axis=1)
	expected_post_STD_result = (A.weighting_factors(my_expected_grid)*numpy.sin(numpy.radians(my_expected_grid))).sum(axis=1)
	expected_prior_STD_result = (A.weighting_factors(my_expected_grid)*numpy.square(numpy.sin(numpy.radians(my_expected_grid)))).sum(axis=1)

	for field, array in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [expected_MAP_result, expected_post_STD_result, expected_prior_STD_result]):
	    numpy.testing.assert_array_equal(array, A.compute_gridded_expected_value(field, component_solution))
    