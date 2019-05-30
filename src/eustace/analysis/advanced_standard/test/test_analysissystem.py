"""Test of example using simulated data, low-resolution mesh, and climatology model."""


import numpy
import scipy.sparse
import unittest

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem

from StringIO import StringIO

from eustace.analysis.mesh.mesh import MeshIcosahedronSubdivision
from eustace.analysis.mesh.geometry import cartesian_to_polar2d
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystemInputLoader
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters
from eustace.analysis.advanced_standard.elements.seasonal import datetime_to_decimal_year
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import Regridder

from eustace.outputformats.definitions import GLOBAL_FIELD_OUTPUT_FLAGS

from eustace.timeutils.epoch import epoch_plus_days

class TestAnalysisSystem(unittest.TestCase):
  class DummyComponentSolution(object):
      
      def __init__(self, sample_size=1):
	self.sample_size=sample_size
      
      def solution_observation_expected_value(self, observations):
	  return numpy.zeros([observations.number_of_observations(),1])+1.
		  
      def solution_observation_expected_uncertainties(self, observations):
	  return numpy.zeros([observations.number_of_observations(),1])+.5
	  
      def solution_observation_prior_uncertainties(self, observations):
	  return numpy.zeros([observations.number_of_observations(),1])+.7

      def solution_observation_projected_sample(self, observations):
	  vector = numpy.zeros([observations.number_of_observations(),self.sample_size])
	  for index in range(self.sample_size):
	      vector[:, index] = -index
	  return vector
	  
  class LessDummyComponentSolution(object):
      
      def __init__(self, sample_size=1):
	self.sample_size=sample_size
      
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

      def solution_observation_projected_sample(self, observations):
	  vector = numpy.zeros([observations.number_of_observations(),self.sample_size])
	  for index in range(self.sample_size):
	      vector[:, index] = index
	  return vector
	  
  class FirstDummyComponent(object):
    
    def __init__(self, sample_size=1):
      self.sample_size = sample_size
    
    def component_solution(self):
      return TestAnalysisSystem.DummyComponentSolution(self.sample_size)
      
  class SecondDummyComponent(object):
    
    def __init__(self, sample_size=1):
      self.sample_size = sample_size
    
    def component_solution(self):
      return TestAnalysisSystem.LessDummyComponentSolution(self.sample_size)
      
  class DummySpatialComponent(SpatialComponent):
    """We need an instance of the SpatialComponent class to check the correct computation of climatology fraction""" 
    
    def component_solution(self):
      return TestAnalysisSystem.LessDummyComponentSolution(self.sample_size)
       
  def setUp(self):
      """Set up mock objects to be used for testing the Regridder class functionalities"""

      self.structure = OutputRectilinearGridStructure(1, None, numpy.array([1, 2, 3, 4]), numpy.array([.2, .3, .4]))

  def test_init(self):
      A = AnalysisSystem('A', 'B', 'C')
      
      self.assertEqual('A', A.components)
      self.assertEqual('B', A.observable)   
      self.assertEqual('C', A.log)

      self.assertListEqual(['POINTWISE', 'GRID_CELL_AREA_AVERAGE'], A.output_flags)
      self.assertListEqual(GLOBAL_FIELD_OUTPUT_FLAGS, A.field_flags)
      
  def test_evaluate_single_component_grid_cell_average_expected_value(self):
      A = AnalysisSystem('A', 'B', 'C')
      
      # Dummy component solution
      component_solution = TestAnalysisSystem.DummyComponentSolution()
      
      for field, value in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [1., .25, .49]):
	  numpy.testing.assert_array_almost_equal(value, A.evaluate_single_component_grid_cell_average_expected_value(field, component_solution, self.structure, [1, 1], 3))

      component_solution = TestAnalysisSystem.LessDummyComponentSolution()
      B =  Regridder(self.structure, [2,3], blocking=3)
 
      # Less dummy component solution
      expected_weight_latitudes = numpy.array([[0.75, 1.25, 0.75, 1.25, 0.75, 1.25], [0.75, 1.25, 0.75, 1.25, 0.75, 1.25],[0.75, 1.25, 0.75, 1.25, 0.75, 1.25],
					       [1.75, 2.25, 1.75, 2.25, 1.75, 2.25], [1.75, 2.25, 1.75, 2.25, 1.75, 2.25],[1.75, 2.25, 1.75, 2.25, 1.75, 2.25],
					       [2.75, 3.25, 2.75, 3.25, 2.75, 3.25], [2.75, 3.25, 2.75, 3.25, 2.75, 3.25],[2.75, 3.25, 2.75, 3.25, 2.75, 3.25],
					       [3.75, 4.25, 3.75, 4.25, 3.75, 4.25], [3.75, 4.25, 3.75, 4.25, 3.75, 4.25],[3.75, 4.25, 3.75, 4.25, 3.75, 4.25]])
					                                        
      expected_MAP_array = numpy.array([[1., .5, 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., 0., -.5], [1., .5, 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])
      expected_post_STD_array = numpy.array([[0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., -.5, 0.], [0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., -.5, 0.]])
      expected_prior_STD_array = numpy.array([[0., 1., .5, 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					      [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					      [0., 0., 0., 0., 0., -.5], [0., 1., .5, 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])


      expected_MAP_result = (B.weighting_factors(expected_weight_latitudes)*expected_MAP_array).sum(axis=1)
      expected_post_STD_result = numpy.square((B.weighting_factors(expected_weight_latitudes)*expected_post_STD_array).sum(axis=1))
      expected_prior_STD_result = numpy.square((B.weighting_factors(expected_weight_latitudes)*expected_prior_STD_array).sum(axis=1))
      for field, array in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [expected_MAP_result, expected_post_STD_result, expected_prior_STD_result]):
	  numpy.testing.assert_array_equal(array, A.evaluate_single_component_grid_cell_average_expected_value(field, component_solution, self.structure, [2, 3], 3))
 
 
  def test_evaluate_grid_cell_average_expected_value(self):
      # Combined components solutions
      A = AnalysisSystem([TestAnalysisSystem.FirstDummyComponent(), TestAnalysisSystem.SecondDummyComponent()], 'B', 'C')

      B =  Regridder(self.structure, [2,3], blocking=3)
      expected_weight_latitudes = numpy.array([[0.75, 1.25, 0.75, 1.25, 0.75, 1.25], [0.75, 1.25, 0.75, 1.25, 0.75, 1.25],[0.75, 1.25, 0.75, 1.25, 0.75, 1.25],
					       [1.75, 2.25, 1.75, 2.25, 1.75, 2.25], [1.75, 2.25, 1.75, 2.25, 1.75, 2.25],[1.75, 2.25, 1.75, 2.25, 1.75, 2.25],
					       [2.75, 3.25, 2.75, 3.25, 2.75, 3.25], [2.75, 3.25, 2.75, 3.25, 2.75, 3.25],[2.75, 3.25, 2.75, 3.25, 2.75, 3.25],
					       [3.75, 4.25, 3.75, 4.25, 3.75, 4.25], [3.75, 4.25, 3.75, 4.25, 3.75, 4.25],[3.75, 4.25, 3.75, 4.25, 3.75, 4.25]])
					                                        
      expected_MAP_array = numpy.array([[1., .5, 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., 0., -.5], [1., .5, 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])
      expected_post_STD_array = numpy.array([[0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., -.5, 0.], [0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., -.5, 0.]])

      expected_MAP_result = (B.weighting_factors(expected_weight_latitudes)*(expected_MAP_array+1.)).sum(axis=1)
      
      first_post_STD = (B.weighting_factors(expected_weight_latitudes)*expected_post_STD_array).sum(axis=1)
      second_post_STD = (B.weighting_factors(expected_weight_latitudes)*.5).sum(axis=1)
      expected_post_STD_result = numpy.sqrt(numpy.square(first_post_STD)+ numpy.square(second_post_STD))

      for field, array in zip(GLOBAL_FIELD_OUTPUT_FLAGS, [expected_MAP_result, expected_post_STD_result]):
	  numpy.testing.assert_array_almost_equal(array, A.evaluate_grid_cell_average_expected_value(field, self.structure, [2, 3], 3))

      self.assertRaises(ValueError, A.evaluate_grid_cell_average_expected_value, GLOBAL_FIELD_OUTPUT_FLAGS[2], self.structure, [2, 3], 3)
      
  def test_evaluate_pointlike_limit_expected_value(self):
      A = AnalysisSystem([TestAnalysisSystem.FirstDummyComponent(), TestAnalysisSystem.SecondDummyComponent()], 'B', 'C')
      
      # Checking only valid flags gets used
      self.assertRaises(ValueError, A.evaluate_expected_value, 'A', 'MAP', 'C')
      self.assertRaises(ValueError, A.evaluate_expected_value, 'POINTWISE', 'B', 'C')
      
      # Checking computations only for MAP and posterior marginal STD:
      for field in GLOBAL_FIELD_OUTPUT_FLAGS[:2]:
	  pointwise_result = A.evaluate_expected_value(field, self.structure, 'POINTWISE')
	  pointwise_limit_result = A.evaluate_expected_value(field, self.structure,  'GRID_CELL_AREA_AVERAGE', [1, 1], 10)
	  numpy.testing.assert_array_equal(pointwise_result, pointwise_limit_result)

      self.assertRaises(ValueError, A.evaluate_expected_value, GLOBAL_FIELD_OUTPUT_FLAGS[2], self.structure, 'POINTWISE')
      self.assertRaises(ValueError, A.evaluate_expected_value, GLOBAL_FIELD_OUTPUT_FLAGS[2], self.structure, 'POINTWISE', [2,3], 3)

  def test_check_spatial_components(self):
    
      A = AnalysisSystem([TestAnalysisSystem.FirstDummyComponent(), TestAnalysisSystem.SecondDummyComponent()], 'B', 'C')      
      self.assertEqual(None, A.check_spatial_components())
      
      A = AnalysisSystem([TestAnalysisSystem.SecondDummyComponent(), SpatialComponent('A', 'B'), SpatialComponent('C', 'D')], 'E', 'F')
      self.assertEqual(None, A.check_spatial_components())
      
      A = AnalysisSystem([TestAnalysisSystem.FirstDummyComponent(), TestAnalysisSystem.SecondDummyComponent(), SpatialComponent('A', 'B')], 'C', 'D')
      self.assertEqual(2, A.check_spatial_components())

      A = AnalysisSystem([SpatialComponent('A', 'B'), TestAnalysisSystem.FirstDummyComponent(), TestAnalysisSystem.SecondDummyComponent()], 'C', 'D')
      self.assertEqual(0, A.check_spatial_components())
     
  def test_evaluate_climatology_fraction(self):
      A = AnalysisSystem([TestAnalysisSystem.FirstDummyComponent(), TestAnalysisSystem.DummySpatialComponent('A', 'B')], 'E', 'F')
    
      B =  Regridder(self.structure, [2,3], blocking=3)
      expected_weight_latitudes = numpy.array([[0.75, 1.25, 0.75, 1.25, 0.75, 1.25], [0.75, 1.25, 0.75, 1.25, 0.75, 1.25],[0.75, 1.25, 0.75, 1.25, 0.75, 1.25],
					       [1.75, 2.25, 1.75, 2.25, 1.75, 2.25], [1.75, 2.25, 1.75, 2.25, 1.75, 2.25],[1.75, 2.25, 1.75, 2.25, 1.75, 2.25],
					       [2.75, 3.25, 2.75, 3.25, 2.75, 3.25], [2.75, 3.25, 2.75, 3.25, 2.75, 3.25],[2.75, 3.25, 2.75, 3.25, 2.75, 3.25],
					       [3.75, 4.25, 3.75, 4.25, 3.75, 4.25], [3.75, 4.25, 3.75, 4.25, 3.75, 4.25],[3.75, 4.25, 3.75, 4.25, 3.75, 4.25]])
					                                        
      expected_post_STD_array = numpy.array([[0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					[0., 0., 0., 0., -.5, 0.], [0., 0., 1., .5, 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., -.5, 0.]])
      expected_prior_STD_array = numpy.array([[0., 1., .5, 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
			  		      [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0.], 
					      [0., 0., 0., 0., 0., -.5], [0., 1., .5, 0., 0., 0.], [0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., -.5]])

      expected_post_STD_result = numpy.square((B.weighting_factors(expected_weight_latitudes)*expected_post_STD_array).sum(axis=1))
      expected_prior_STD_result = numpy.square((B.weighting_factors(expected_weight_latitudes)*expected_prior_STD_array).sum(axis=1))
      # numpy is raising a warning due to 0/0 division, but is perfectly fine, as some of the prior std are set to zero in this mock test
      expected_climatology_fraction = expected_post_STD_result/expected_prior_STD_result
      
      numpy.testing.assert_array_almost_equal(expected_climatology_fraction, A.evaluate_climatology_fraction(self.structure, [2, 3], 3))
    
      # Checking pointlike limit
      B =  Regridder(self.structure, [1,1], blocking=20)
					                                        
      expected_post_STD_array = numpy.array([0., 0., 1., .5, 0., 0., 0., 0., 0., 0., -.5, 0.])
      expected_prior_STD_array = numpy.array([0., 1., .5, 0., 0., 0., 0., 0., 0., 0., 0., -.5])
      expected_post_STD_result = numpy.square(expected_post_STD_array)
      expected_prior_STD_result = numpy.square(expected_prior_STD_array)
      expected_climatology_fraction = expected_post_STD_result/expected_prior_STD_result
      
      numpy.testing.assert_array_almost_equal(expected_climatology_fraction, A.evaluate_climatology_fraction(self.structure, [1, 1], 20))
      
  def test_evaluate_sample(self):
      # Testing the system is checking all model components have drawn sample with the same sample size
      A = AnalysisSystem([TestAnalysisSystem.FirstDummyComponent(3), TestAnalysisSystem.DummySpatialComponent('A', 'B', sample_size=2)], 'E', 'F')
      self.assertRaises(ValueError, A.evaluate_projected_sample, self.structure)

      # The first model component generates a sample of the form [[0, 0, ..], [-1, -1, ...], [-(sample_size-1), -(sample_size-1),...]]
      # The second model component generates a sample of the form [[0, 0, ..], [1, 1, ...], [(sample_size-1), (sample_size-1),...]]
      # By summing the projection we should obtained a vector of zeros.
      A = AnalysisSystem([TestAnalysisSystem.FirstDummyComponent(5), TestAnalysisSystem.DummySpatialComponent('A', 'B', sample_size=5)], 'E', 'F')
      numpy.testing.assert_array_equal(0., A.evaluate_projected_sample(self.structure))