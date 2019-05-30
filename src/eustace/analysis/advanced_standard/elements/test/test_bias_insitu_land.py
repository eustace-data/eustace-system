"""Tests for bias insitu land element."""

import datetime
import unittest
import numpy
import scipy.sparse
import tempfile

from eumopps.timeutils.timebase import TimeBaseDays
from eustace.analysis.advanced_standard.elements.bias_insitu_land import InsituLandBiasElement
from eustace.analysis.advanced_standard.elements.bias_insitu_land import InsituLandBiasDesign
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT
from eustace.analysis.advanced_standard.fileio.observation_structure_source_connector import ObservationStructureSourceConnector
from eustace.analysis.fileio.observationsource_rawbinary import ObservableFileSpec
from eustace.analysis.observationsource import Observations
from eustace.analysis.observationsource import ObservationsBreakPoints
from eustace.analysis.observationsource import ObservationSource
from eustace.outputformats import definitions 
from netCDF4 import Dataset
                   
class TestInsituLandBiasElement(unittest.TestCase):

    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 23.0, 3.2 ], [ 15.0, -7.0 ], [  5.0, 100.0 ] ])

        def covariate_effect(self, groupname, breakpoints):

            if groupname == InsituLandBiasElement.GROUPNAME:

                # obs 1 --> bias 0
                # obs 2 --> bias 4
                return numpy.array( [ [ 1, 0 ],
                                      [ 2, 4 ] ], numpy.int64 )

            elif groupname == 'NotBob':

                # obs 0 --> bias 2
                return numpy.array( [ [ 0, 2 ] ], numpy.int64 )

            else:

                # (empty map)
                return numpy.array( [ [ ] ], numpy.int64)

        def number_of_observations(self):

            return self.location_polar_coordinates().shape[0]

    def setUp(self):
	self.breakpoints_file = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.elements.test_bias_insitu_land', suffix='.nc',delete=True)

	self.break_stations = numpy.array([4, 4, 4, 4, 6, 6, 7, 9 , 9 , 9], dtype=numpy.int32)
	self.break_times = numpy.array([47116, 51499, 53325, 54056, 38715, 39811, 47116, 23740, 24836, 27027], dtype=numpy.int32)
	self.break_likelihood = numpy.array([1, 8, 7, 10, 11, 2, 23, 4, 7, 10], dtype=numpy.int8)
	self.detection_feasibility = numpy.array(range(1,21), dtype=numpy.int8)
	dataset=Dataset(self.breakpoints_file.name,'w','NETCDF4')
	dataset.createDimension('merged_break',size=10)
	dataset.createDimension('station',size=len(self.detection_feasibility))
	merged_break_time=dataset.createVariable('merged_break_time',numpy.int32,dimensions=('merged_break'))
	merged_break_station=dataset.createVariable('merged_break_station',numpy.int32,dimensions=('merged_break'))
	merged_break_likelihood=dataset.createVariable('merged_break_likelihood',numpy.int8,dimensions=('merged_break'))
	detection_feasibility=dataset.createVariable('detection_feasibility',numpy.int8,dimensions=('station'))
	merged_break_time[:]=self.break_times
	merged_break_station[:]=self.break_stations
	merged_break_likelihood[:]=self.break_likelihood
	detection_feasibility[:]=self.detection_feasibility
	dataset.close()

    def test_init(self):

	# Basic check: is the InsituBiasElement class initialized with the right input file?
        self.assertRaises(IOError,InsituLandBiasElement,'Bob')
	
	# Using the right input file
	bias = InsituLandBiasElement(self.breakpoints_file.name)
	self.assertEqual(bias.groupname, InsituLandBiasElement.GROUPNAME)
        self.assertEqual(bias.number_of_biases, len(self.break_times))
	self.assertTrue(isinstance(bias.observed_breakpoints,ObservationsBreakPoints))
	
	# Using the right input file, enhancing reduction
	bias = InsituLandBiasElement(self.breakpoints_file.name, apply_policy=True, cut_value = 10)
	self.assertEqual(bias.groupname, InsituLandBiasElement.GROUPNAME)
        self.assertEqual(bias.number_of_biases, 4)
	self.assertTrue(isinstance(bias.observed_breakpoints,ObservationsBreakPoints))


    def test_element_design(self):

        design = InsituLandBiasElement(self.breakpoints_file.name).element_design(TestInsituLandBiasElement.SimulatedObservationStructure())

        # Check type
        self.assertTrue(isinstance(design, InsituLandBiasDesign))

        # Covariate info should be copied from constructor
        self.assertEqual(InsituLandBiasElement.GROUPNAME, design.groupname)
        self.assertEqual(10, design.number_of_biases)

        # And effect info comes from observation structure
        numpy.testing.assert_equal(design.effect, [ [ 1, 0 ], [ 2, 4 ] ])


    def test_element_prior(self):

        prior = InsituLandBiasElement(self.breakpoints_file.name).element_prior(CovariateHyperparameters(9.9))

        self.assertTrue(isinstance(prior, CovariatePrior))
        self.assertEqual(9.9, prior.hyperparameters.value)

        # The number of state parameters should be the total number of biases (irrespective of number observed)
        self.assertEqual(10, prior.number_of_state_parameters)

class TestInsituLandBiasDesign(unittest.TestCase):

    def test_init(self):

        design_insitu_land = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), InsituLandBiasElement.GROUPNAME, 5, None)
        self.assertEqual(InsituLandBiasElement.GROUPNAME, design_insitu_land.groupname)
        self.assertEqual(5, design_insitu_land.number_of_biases)
        self.assertEqual(3, design_insitu_land.number_of_observations)
        numpy.testing.assert_equal(design_insitu_land.effect, [ [ 1, 0 ], [ 2, 4 ] ])

        design_notbob = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), 'NotBob', 10, None)
        self.assertEqual('NotBob', design_notbob.groupname)
        self.assertEqual(10, design_notbob.number_of_biases)
        self.assertEqual(3, design_notbob.number_of_observations)
        numpy.testing.assert_equal(design_notbob.effect, [ [ 0, 2 ] ])

        design_nobody = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), 'Nobody', 7, None)
        self.assertEqual('Nobody', design_nobody.groupname)
        self.assertEqual(7, design_nobody.number_of_biases)
        self.assertEqual(3, design_nobody.number_of_observations)
        numpy.testing.assert_equal(design_nobody.effect, [ [ ] ])

    def test_design_number_of_state_parameters(self):

        design = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), 'Nobody', 7, None)
        self.assertEqual(7, design.design_number_of_state_parameters())


    def test_design_matrix(self):

        A = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), InsituLandBiasElement.GROUPNAME, 5, None).design_matrix()
        self.assertEqual(SPARSEFORMAT, A.getformat())
        self.assertEqual((3, 5), A.shape)
        self.assertEqual(2, A.nnz)

        # This mapping comes from observation definition
        # in TestBiasElement.SimulatedObservationStructure:
        # obs 1 --> bias 0, obs 2 --> bias 4
        numpy.testing.assert_equal(A.todense(), 
                                   [ [ 0.0, 0.0, 0.0, 0.0, 0.0 ],
                                     [ 1.0, 0.0, 0.0, 0.0, 0.0 ],
                                     [ 0.0, 0.0, 0.0, 0.0, 1.0 ] ])

    def test_design_function(self):

        design = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), InsituLandBiasElement.GROUPNAME, 5, None)
        y = design.design_function(numpy.array([ [ 22.2 ], [ 33.3 ], [ 44.4 ], [ 55.5 ], [ 66.6 ] ]))

        # Again this comes from observation definition
        # in TestBiasElement.SimulatedObservationStructure:
        # obs 0 has no associated bias (so is zero)
        # obs 1 --> bias 0 (which is 22.2)
        # obs 2 --> bias 4 (which is 66.6)
        numpy.testing.assert_almost_equal(y, numpy.array([ [ 0.0 ], [ 22.2 ], [ 66.6 ] ]))

    def test_isnonlinear(self):
        
        self.assertFalse(InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), InsituLandBiasElement.GROUPNAME, 5, None).isnonlinear())

    def test_design_jacobian(self):

        A = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), InsituLandBiasElement.GROUPNAME, 5, None).design_matrix()
        J = InsituLandBiasDesign(TestInsituLandBiasElement.SimulatedObservationStructure(), InsituLandBiasElement.GROUPNAME, 5, None).design_jacobian(numpy.array([ ]))
        self.assertEqual(SPARSEFORMAT, J.getformat())
        numpy.testing.assert_almost_equal(J.todense(), A.todense())
        
class TestInsituLandBias(unittest.TestCase):
  
    class SimulatedObservationSource(ObservationSource):
      """Mock observation source class: used to probe the correctness of the ObservationStructureConnector"""

      OBSERVATIONMAPS = { ObservationSource.TMIN: definitions.TASMIN, ObservationSource.TMEAN: definitions.TAS, ObservationSource.TMAX: definitions.TASMAX }

      def __init__(self, times, look_up, masks, observations_collection, daynumber = None):
	  self.times = times
	  self.look_up = look_up
	  self.masks = masks
	  self.observations_collection = observations_collection                                                     
	  self.daynumber = daynumber
          self.filespecs = {'Tmin': ObservableFileSpec('insitu_land_Tmin_20140202.bin',None), 'Tmax':ObservableFileSpec('insitu_land_Tmax_20140202.bin',None)}
          
      def observables(self):
	  """The names of variables estimated from this source."""
	  return SimulatedObservationSource.OBSERVATIONMAPS.keys()

      def number_of_observations(self):
	  """Total number of observations."""

	  if self.daynumber is None:
	      return len(self.times)*self.look_up.shape[1]
	  else:
	      return self.look_up.shape[1]

      def observation_location_lookup(self):

	  return self.look_up

      def observations(self, observable):

	  variable = self.observations_collection[observable]
	  mask = self.masks[observable]

	  # Restrict by time if requested
	  if self.daynumber is None:

	      measurement = variable
	      time = self.times
	  else:

	      dayindex = numpy.argwhere(self.times == self.daynumber)

	      if len(dayindex) == 1:
		  measurement = variable[dayindex.ravel(), :]
		  mask = self.masks[observable][dayindex.ravel(), :]
		  time = self.times[dayindex]
	      else:
		  return None                        

	  # Formulate output
	  all_time = numpy.tile(time, (measurement.shape[1], 1)).transpose()
	  all_locations = numpy.tile(numpy.array(range(measurement.shape[1]), numpy.uint64), (measurement.shape[0], 1))
	  uncorrelatederror = numpy.zeros(measurement.shape)+1.
	  locallycorrelatederror = [ ]
	  return Observations(
	      mask.ravel(),
	      all_time.ravel(),
	      all_locations.ravel(),
	      measurement.ravel().astype(numpy.float32),
	      uncorrelatederror.ravel(),
	      locallycorrelatederror)

      def local_correlation_length_scale(self, observable):
	  """Length scale for locally correlated component of uncertainty (empty array in this case)."""
	  pass

    def setUp(self):
      """We will have the following setup:
	  
	  9 stations (locations)
	  10 time measurements
	  3 different temperatures series (MIN, MAX, MEAN)
	  3 different masks
      """
      self.times = numpy.array([23740, 24836, 27027, 38715, 39811, 47116, 47117, 51499, 53325, 54056])
      self.observations = {'TMIN': numpy.array([[ 41,  44,  56, 234, 118, 123,  47,  26, 243],        
						[125, 160, 115,  37, 150, 201,  98, 160,  66],
						[185, 118, 271, 232, 160, 198, 257,  30,  34],
						[292,  12,  86, 215,  23, 268,  21, 170, 152],
						[178, 185,  68, 133,  11, 153, 187,  46, 211],
						[183, 290, 298,  70, 295, 172, 268, 217, 122],
						[ 75, 239, 236,  56, 222,  26,  79, 293,  26],
						[268, 268,  76, 192, 113, 120,  71, 210,  88],
						[ 91, 257, 169,   6, 206,   1,  53, 265, 161],
						[148,  30, 251,  56, 188,  91, 276,  98, 190]]),
			    'TMEAN': numpy.array([[197, 298, 287,  55, 276, 238, 138, 131,   3],        
						  [232, 262, 141,  23,  68, 190,  50,  85, 209],
						  [ 62,  41, 119,  95,  80, 141, 117, 183,  51],
						  [148, 137, 298,  76, 284, 220, 151, 198, 224],
						  [148, 129, 214, 191, 230,  70, 288, 227, 178],
						  [ 39, 144, 273, 155, 247,  85,  57, 132, 187],
						  [281, 271, 274, 131,  54, 299,   9,  33, 111],
						  [  1, 216, 235,  61, 131,  50,  57, 247, 288],
						  [125, 223,  71, 252, 180,  14,   1,  60, 269],
						  [ 85,  74, 260, 229, 144,  32, 185,  28,  46]]),
			      'TMAX': numpy.array([[224, 101, 110, 276,  49,  95, 200,  90,  23],        
						  [ 93, 265, 162, 259,   4,  40,  61, 252,  18],
						  [ 85,  22, 209, 170, 121,  52, 299,   9,  61],
						  [  7, 227, 164, 194, 160, 140, 116,  95,  89],
						  [270,  66, 219, 211, 214, 133, 250, 119, 117],
						  [  1,  53, 134, 128, 193, 164,  13, 117, 257],
						  [235, 208, 181, 220,  68,  83,  90,  53,  28],
						  [258, 171,  48, 271, 134, 111,  40, 198,  99],
						  [232,  36, 238, 259,  79, 163, 238,  11,  87],
						  [ 48,  40,  93,  46,  98, 269,  69, 166,  35]])}

      self.masks = {'TMIN': numpy.array([[False, False, False,  True, False, False,  True,  True,  True], 
					  [False, False, False,  True,  True,  True, False, False,  True],
					  [ True, False, False,  True, False,  True,  True,  True,  True],
					  [ True,  True,  True, False, False,  True, False,  True,  True],
					  [ True, False,  True,  True, False,  True, False, False, False],
					  [ True, False,  True,  True,  True, False, False, False,  True],
					  [ True, False, False, False,  True, False, False,  True, False],
					  [ True,  True,  True, False, False, False,  True,  True, False],
					  [ True, False,  True,  True, False, False, False, False, False],
					  [False,  True,  True, False, False, False,  True,  True, False]]),
		    'TMEAN': numpy.array([[ True,  True, False,  True,  True, False, False,  True, False], 
					  [False,  True, False, False, False, False, False,  True,  True],
					  [ True, False,  True,  True, False, False, False, False,  True],
					  [False, False,  True,  True, False,  True,  True, False,  True],
					  [False,  True,  True, False, False, False,  True, False,  True],
					  [False, False,  True, False,  True, False,  True, False,  True],
					  [ True, False, False,  True, False,  True,  True, False,  True],
					  [False,  True,  True, False, False,  True, False, False, False],
					  [False,  True, False,  True,  True, False, False,  True, False],
					  [False, False,  True, False,  True,  True,  True,  True,  True]]),
		      'TMAX': numpy.array([[ True, False, False, False, False, False,  False, False,  False], 
					  [False, False,  True, False,  True, False,  True,  True, False],
					  [False, False, False,  True,  True, False,  True, False, False],
					  [ True,  True, False,  True,  True, False,  True, False, False],
					  [False, False,  True, False, False, False, False, False,  True],
					  [ True,  True,  True,  True,  True, False,  True, False, False],
					  [ True, False,  True, False, False, False, False, False,  True],
					  [False,  True,  True,  True,  True,  True,  True,  True, False],
					  [False, False,  True,  True,  True, False, False,  True, False],
					  [False, False,  True, False, False, False, False, False,  True]])}

      self.look_up = numpy.array([[-132.18908084,  -66.80417394,  -58.83865167,   88.72957003, 38.73584399, -144.16397268,   11.25863152,  150.2133835 , 107.26345795],
				  [-65.33291349,  70.30566194,  69.27239724, -44.01245973, -30.80704767,  64.12609839,  59.4841331 ,  62.79291718, -25.13326873]])

      self.breakpoints_file = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.elements.test_bias_insitu_land', suffix='.nc',delete=True)

      self.break_stations = numpy.array([4, 4, 4, 4, 6, 6, 7, 8 , 8 , 8], dtype=numpy.int32)
      self.break_times = numpy.array([47116, 51499, 53325, 54056, 38715, 39811, 47116, 23740, 24836, 27027], dtype=numpy.int32)
      self.break_likelihood = numpy.array([1, 8, 7, 10, 11, 2, 23, 4, 7, 10], dtype=numpy.int8)
      self.detection_feasibility = numpy.array(range(1,21), dtype=numpy.int8)
      dataset=Dataset(self.breakpoints_file.name,'w','NETCDF4')
      dataset.createDimension('merged_break',size=10)
      dataset.createDimension('station',size=len(self.detection_feasibility))
      merged_break_time=dataset.createVariable('merged_break_time',numpy.int32,dimensions=('merged_break'))
      merged_break_station=dataset.createVariable('merged_break_station',numpy.int32,dimensions=('merged_break'))
      merged_break_likelihood=dataset.createVariable('merged_break_likelihood',numpy.int8,dimensions=('merged_break'))
      detection_feasibility=dataset.createVariable('detection_feasibility',numpy.int8,dimensions=('station'))
      merged_break_time[:]=self.break_times
      merged_break_station[:]=self.break_stations
      merged_break_likelihood[:]=self.break_likelihood
      detection_feasibility[:]=self.detection_feasibility
      dataset.close()

    def test_bias_insitu_land(self):

      bias = InsituLandBiasElement(self.breakpoints_file.name)
      
      index = 9
      obs = 'TMIN'
      # valid observation indices = [0, 3, 4, 5, 8]
      # t = 54056 -> [[2, 3]]
      # resulting design matrix
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
	    
      test_source = TestInsituLandBias.SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
      test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))

      bias_design = bias.element_design(test_connector)
      self.assertEqual(bias_design.groupname, 'insitu_land')
      self.assertEqual(bias_design.number_of_biases, len(self.break_times))
      self.assertEqual(bias_design.number_of_observations, (~self.masks[obs][index,:]).sum())
      numpy.testing.assert_array_equal(bias_design.effect, numpy.array([[2, 3]]))
      numpy.testing.assert_array_equal(bias_design.design_matrix().todense(), 
				       numpy.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]]))
      
      index = 4
      obs = 'TMEAN'
      # valid observation indices = [0, 3, 4, 5, 7]
      # t = 39811 -> [[2, 0], [4, 6]]
      # resulting design matrix
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[1. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
      
      test_source = TestInsituLandBias.SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
      test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))

      bias_design = bias.element_design(test_connector)
      self.assertEqual(bias_design.groupname, 'insitu_land')
      self.assertEqual(bias_design.number_of_biases, len(self.break_times))
      self.assertEqual(bias_design.number_of_observations, (~self.masks[obs][index,:]).sum())
      numpy.testing.assert_array_equal(bias_design.effect, numpy.array([[2, 0], [4, 6]]))
      numpy.testing.assert_array_equal(bias_design.design_matrix().todense(), 
				       numpy.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 1., 0., 0., 0.]]))
						    
						    
						    
      index = 9
      obs = 'TMAX'
      # valid observation indices = [0, 1, 3, 4, 5, 6, 7]
      # t = 54056 -> [[3, 3]]
      # resulting design matrix
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 1. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      
      test_source = TestInsituLandBias.SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
      test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))

      bias_design = bias.element_design(test_connector)
      self.assertEqual(bias_design.groupname, 'insitu_land')
      self.assertEqual(bias_design.number_of_biases, len(self.break_times))
      self.assertEqual(bias_design.number_of_observations, (~self.masks[obs][index,:]).sum())
      numpy.testing.assert_array_equal(bias_design.effect, numpy.array([[3, 3]]))
      numpy.testing.assert_array_equal(bias_design.design_matrix().todense(), 
				       numpy.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 1., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]]))
						    

      index = 0
      obs = 'TMAX'
      # valid observation indices = [ 1, 2, 3, 4, 5, 6, 7, 8]
      # t = 23740 -> [[3, 0], [5, 4], [6, 6], [7, 7]]
      # resulting design matrix
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[1. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 1. 0. 0. 0.]
      #[0. 0. 0. 0. 0. 0. 0. 1. 0. 0.]
      
      test_source = TestInsituLandBias.SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
      test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))

      bias_design = bias.element_design(test_connector)
      self.assertEqual(bias_design.groupname, 'insitu_land')
      self.assertEqual(bias_design.number_of_biases, len(self.break_times))
      self.assertEqual(bias_design.number_of_observations, (~self.masks[obs][index,:]).sum())
      numpy.testing.assert_array_equal(bias_design.effect, numpy.array([[3, 0], [5, 4], [6, 6], [7, 7]]))
      numpy.testing.assert_array_equal(bias_design.design_matrix().todense(), 
				       numpy.array([[0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [1., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 1., 0., 0., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 1., 0., 0., 0.],
						    [0., 0., 0., 0., 0., 0., 0., 1., 0., 0.]]))

