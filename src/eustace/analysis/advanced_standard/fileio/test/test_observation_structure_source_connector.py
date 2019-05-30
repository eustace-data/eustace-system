"""Testing the connector, for several different cases"""

import datetime
import numpy
import scipy.sparse
import tempfile
import unittest

from eumopps.timeutils.timebase import TimeBaseDays
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.fileio.observation_structure_source_connector import ObservationStructureSourceConnector
from eustace.analysis.fileio.observationsource_rawbinary import ObservableFileSpec
from eustace.analysis.observationsource import ObservationSource 
from eustace.analysis.observationsource import Observations
from eustace.outputformats import definitions 
from eustace.preprocess.fileio.insitu_land_breakpoints import ObservationBreakPointSourceInsituLand
from netCDF4 import Dataset

class SimulatedObservationSource(ObservationSource):
    """Mock observation source class: used to probe the correctness of the ObservationStructureConnector"""

    OBSERVATIONMAPS = { ObservationSource.TMIN: definitions.TASMIN, ObservationSource.TMEAN: definitions.TAS, ObservationSource.TMAX: definitions.TASMAX }

    def __init__(self, times, look_up, masks, observations_collection, daynumber = None):
	self.times = times
	self.look_up = look_up
	self.masks = masks
	self.observations_collection = observations_collection
	self.daynumber = daynumber
	self.filespecs = {SimulatedObservationSource.OBSERVATIONMAPS[ObservationSource.TMIN]: ObservableFileSpec('some_file',None)}

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

class TestObservationStructureSourceConnector(unittest.TestCase):

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
		       'TMAX': numpy.array([[ True, False, False, False, False, False,  True, False,  True], 
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

     def test_simulated_observation_source(self):  
	"""Just check out simulated version of ObservationSource class is correctly working"""

	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations)
	numpy.testing.assert_array_equal(self.times, test_source.times)
	numpy.testing.assert_array_equal(self.look_up, test_source.look_up)
	numpy.testing.assert_array_equal(self.masks, test_source.masks)
	numpy.testing.assert_array_equal(self.observations, test_source.observations_collection)

	self.assertItemsEqual(SimulatedObservationSource.OBSERVATIONMAPS.keys(), test_source.observables())
	self.assertEqual(test_source.number_of_observations(), 90)

	numpy.testing.assert_array_equal(self.look_up, test_source.observation_location_lookup())

	obs = 'TMAX'
	observation = test_source.observations(obs)
	numpy.testing.assert_array_equal(self.times[1], observation.time[9:18])
	numpy.testing.assert_array_equal(self.masks[obs][1,:], observation.mask[9:18])
	numpy.testing.assert_array_equal(range(0,9),observation.location[18:27])
	numpy.testing.assert_array_equal(self.observations[obs][2,:], observation.measurement[18:27])

	# select a specific day = 27027
	index = 6
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	numpy.testing.assert_array_equal(self.times, test_source.times)
	numpy.testing.assert_array_equal(self.look_up, test_source.look_up)
	numpy.testing.assert_array_equal(self.masks, test_source.masks)
	numpy.testing.assert_array_equal(self.observations, test_source.observations_collection)

	self.assertItemsEqual(SimulatedObservationSource.OBSERVATIONMAPS.keys(), test_source.observables())
	self.assertEqual(test_source.number_of_observations(), 9)

	numpy.testing.assert_array_equal(self.look_up, test_source.observation_location_lookup())

	obs = 'TMEAN'
	observation = test_source.observations(obs)
	numpy.testing.assert_array_equal(self.times[index], observation.time)
	numpy.testing.assert_array_equal(self.masks[obs][index,:], observation.mask)
	numpy.testing.assert_array_equal(range(0,9),observation.location)
	numpy.testing.assert_array_equal(self.observations[obs][index,:], observation.measurement)

     def test_init(self):
	index = 7
	obs='TMIN'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))

	numpy.testing.assert_array_equal(self.look_up, test_connector.location_lookup)
      	numpy.testing.assert_array_equal(numpy.nonzero(~self.masks[obs][index,:]), test_connector.valid_indices)
	self.assertEqual(TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]), test_connector.corresponding_datetime)
	self.assertEqual(['some_file'], test_connector.observable_filenames)
	self.assertTrue(isinstance(test_connector.observations,Observations))

     def test_time_index(self):
	index = 8
	obs='TMEAN'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	numpy.testing.assert_array_equal(self.times[index], test_connector.time_index())

     def test_time_datetime(self):
	index = 3
	obs='TMAX'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	self.assertEqual(TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]), test_connector.time_datetime())

     def test_number_of_observations(self):
	index = 4
	obs = 'TMAX'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	self.assertEqual((~self.masks[obs][index,:]).sum(), test_connector.number_of_observations()) 

     def test_location_polar_coordinates(self):
	index = 4
	obs = 'TMIN'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	numpy.testing.assert_array_equal(self.look_up[:,[1, 4, 6, 7, 8]].T , test_connector.location_polar_coordinates())

     def test_observation_vector(self):
	index = 6
	obs = 'TMEAN'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	numpy.testing.assert_array_equal(self.observations['TMEAN'][index,[1, 2, 4, 7]], test_connector.observation_vector())

     def test_observation_precision(self):
	index = 6
	obs = 'TMEAN'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	test_array = numpy.zeros(9)+1.
	test_matrix = scipy.sparse.diags(1.0 / (test_array[[1,2,4,7]]**2), format='csc')
	numpy.testing.assert_array_equal(test_matrix.toarray(),test_connector.observation_precision().toarray())

     def test_check_all_candidates(self):
	index = 6
	obs = 'TMEAN'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	
	pattern = '%a%C%B.b1.%Y%m%d.000000.cdf'
	self.assertFalse(test_connector.check_all_candidates(pattern))
	
	pattern = 'some_file'
	self.assertTrue(test_connector.check_all_candidates(pattern))

	test_connector.observable_filenames = ['inZitu_land_Tmin_20140202.bin', 'insitu_land_Tmax_20140202.bin']
	pattern = 'insitu_land_%A_%Y%m%d.bin'
	self.assertFalse(test_connector.check_all_candidates(pattern))


	test_connector.observable_filenames = ['insitu_land_Tmin_20140202.bin', 'insitu_land_Tmax_20140202.bin']
	pattern = 'insitu_land_%A_%Y%m%d.bin'
	self.assertTrue(test_connector.check_all_candidates(pattern))
	
     def test_covariate_effect_exception(self):
	index = 6
	obs = 'TMEAN'
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	self.assertRaises(ValueError, test_connector.covariate_effect,'gababob')
	self.assertEqual(None, test_connector.covariate_effect('insitu_land'))

     def test_insitu_land_covariate_effect(self):
       	reader = ObservationBreakPointSourceInsituLand(self.breakpoints_file.name)
	observed_breakpoints = reader.observations('merged_break')
	
	# Testing covariate effect for all the available observables
	# break points =  [47116, 51499, 53325, 54056, 38715, 39811, 47116, 23740, 24836, 27027]
	# stations     =  [4,     4,     4,     4,     6,     6,     7,     8 ,    8 ,    8]

	index = 9
	obs = 'TMIN'
	# valid observation indices = [0, 3, 4, 5, 8]
	# t = 54056 -> [[2, 3]]
	
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	test_connector.observable_filenames = ['insitu_land_Tmin_20140202.bin', 'insitu_land_Tmax_20140202.bin']
	numpy.testing.assert_array_equal(test_connector.covariate_effect('insitu_land',**{'breakpoints': observed_breakpoints}), numpy.array([[2,3]]))


	index = 6
	obs = 'TMEAN'
	# valid observation indices = [1, 2, 4, 7]
	# t = 47117 -> [[2,1]]
	
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	test_connector.observable_filenames = ['insitu_land_Tmin_20140202.bin', 'insitu_land_Tmax_20140202.bin']
	numpy.testing.assert_array_equal(test_connector.covariate_effect('insitu_land',**{'breakpoints': observed_breakpoints}), numpy.array([[2,1]]))

	index = 2
	obs = 'TMAX'
	# valid observation indices = [0, 1, 2, 5, 7, 8]
	# t = 27027 -> [[4, 6], [5, 9]]
	
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	test_connector.observable_filenames = ['insitu_land_Tmin_20140202.bin', 'insitu_land_Tmax_20140202.bin']
	numpy.testing.assert_array_equal(test_connector.covariate_effect('insitu_land',**{'breakpoints': observed_breakpoints}), numpy.array([[4, 6], [5, 9]]))

     def test_global_satellite_effect(self):

	index = 9
	obs = 'TMIN'
	# valid observation indices = [0, 3, 4, 5, 8]
	
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	test_connector.observable_filenames = ['surfaceairmodel_land_Tmin_20140202.bin', 'surfaceairmodel_land_Tmax_20140202.bin']
	numpy.testing.assert_array_equal(test_connector.covariate_effect('surfaceairmodel_land_global'), numpy.array([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0]]))


	index = 6
	obs = 'TMEAN'
	# valid observation indices = [1, 2, 4, 7]
	
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	test_connector.observable_filenames = ['surfaceairmodel_ocean_Tmean_20140202.bin']
	numpy.testing.assert_array_equal(test_connector.covariate_effect('surfaceairmodel_ocean_global'), numpy.array([[0, 0], [1, 0], [2, 0], [3, 0]]))

	index = 2
	obs = 'TMAX'
	# valid observation indices = [0, 1, 2, 5, 7, 8]
	
	test_source = SimulatedObservationSource(self.times, self.look_up, self.masks, self.observations, daynumber=self.times[index])
	test_connector = ObservationStructureSourceConnector(test_source, obs, TimeBaseDays(datetime.datetime(1850,1,1)).number_to_datetime(self.times[index]))
	test_connector.observable_filenames = ['surfaceairmodel_ice_Tmin_20140202.bin', 'surfaceairmodel_ice_Tmean_20140202.bin', 'surfaceairmodel_ice_Tmax_20140202.bin']
	numpy.testing.assert_array_equal(test_connector.covariate_effect('surfaceairmodel_ice_global'), numpy.array([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [5, 0]]))
