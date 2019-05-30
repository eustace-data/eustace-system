"""Tests for input observation sources."""

import unittest
import numpy
#from ..observationsource import Observations
#from ..observationsource import ObservationsBreakPoints
#from ..observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from eustace.analysis.observationsource import ObservationsBreakPoints
from eustace.analysis.observationsource import ObservationSource
from datetime import datetime

class TestObservations(unittest.TestCase):

    def test_init_and_count(self):

        obs = Observations(           
           numpy.array( [ False, True, False, False ], numpy.bool ),
           numpy.float32( 51000.0 ),
           numpy.array( [  1   ,   3  ,    0   ,   2   ], numpy.uint64),
           numpy.array( [ 272.3, 278.4, 271.111, 288.0 ], numpy.float32 ),
           numpy.array( [ 8.0, 9.33, 0.1, 2.2222 ], numpy.float32 ),
           [ numpy.array( [ 1.0, 3.3, 0.2, 7.0 ], numpy.float32 ),
             numpy.array( [ 2.0, 0.1, 0.0, 1.0 ], numpy.float32 ) ] )

        numpy.testing.assert_equal([ False, True, False, False ], obs.mask)
        self.assertEqual(numpy.float32( 51000.0 ), obs.time)
        numpy.testing.assert_equal([  1   ,   3  ,    0   ,   2   ], obs.location)
        numpy.testing.assert_almost_equal([ 272.3, 278.4, 271.111, 288.0 ], obs.measurement, decimal=4)
        numpy.testing.assert_almost_equal([ 8.0, 9.33, 0.1, 2.2222 ], obs.uncorrelatederror, decimal=4)
        self.assertEqual(2, len(obs.locallycorrelatederror))
        numpy.testing.assert_almost_equal([ 1.0, 3.3, 0.2, 7.0 ], obs.locallycorrelatederror[0], decimal=4)
        numpy.testing.assert_almost_equal([ 2.0, 0.1, 0.0, 1.0 ], obs.locallycorrelatederror[1], decimal=4)
        self.assertEqual(4, obs.number_of_observations())

    def test_mean(self):
        # masked data at common locations but different masks
        a = Observations(
           numpy.array( [ False, True, False, False ], numpy.bool ),
           numpy.float32( 52000.0 ),
           numpy.array( [  0   ,   1  ,    2   ,   3   ], numpy.uint64),
           numpy.array( [ 270, 280, 290, 300 ], numpy.float32 ),
           numpy.array( [ 3, 6, 4, 2 ], numpy.float32 ),
           [ numpy.array( [ 4, 5, 3, 1 ], numpy.float32 ),
             numpy.array( [ 2, 3, 5, 6 ], numpy.float32 ) ] )

        b = Observations(
           numpy.array( [ False, False, False, True ], numpy.bool ),
           numpy.float32( 52000.0 ),
           numpy.array( [  0   ,   1  ,    2   ,   3   ], numpy.uint64),
           numpy.array( [ 280, 280, 300, 310 ], numpy.float32 ),
           numpy.array( [ 4, 6, 3, 2 ], numpy.float32 ),
           [ numpy.array( [ 3, 5, 4, 1 ], numpy.float32 ),
             numpy.array( [ 8, 3, 10, 3 ], numpy.float32 ) ] )

        c = Observations.mean(a, b)

        numpy.testing.assert_equal([False, True, False, True], c.mask)
        numpy.testing.assert_equal(numpy.float32(52000.0), c.time)
        numpy.testing.assert_equal([  0   ,   1  ,    2   ,   3   ], c.location)
        numpy.testing.assert_almost_equal(275.0, c.measurement[0], decimal=4)
        numpy.testing.assert_almost_equal(295.0, c.measurement[2], decimal=4)
        numpy.testing.assert_almost_equal(3.5355, c.uncorrelatederror[0], decimal=4)
        numpy.testing.assert_almost_equal(3.5355, c.uncorrelatederror[2], decimal=4)
        numpy.testing.assert_almost_equal(3.5355, c.locallycorrelatederror[0][0], decimal=4)
        numpy.testing.assert_almost_equal(3.5355, c.locallycorrelatederror[0][2], decimal=4)
        numpy.testing.assert_almost_equal(5.8310, c.locallycorrelatederror[1][0], decimal=4)
        numpy.testing.assert_almost_equal(7.9057, c.locallycorrelatederror[1][2], decimal=4)

    def test_mean_unstructured_locations(self):
        # data at different number of locations with some values masked
        a = Observations(
           numpy.array( [ False, True, False, False ], numpy.bool ),
           numpy.float32( 52000.0 ),
           numpy.array( [  0   ,   1  ,  3,  5   ], numpy.uint64),
           numpy.array( [ 270, 280, 290, 300, ], numpy.float32 ),
           numpy.array( [ 3, 6, 4, 2 ], numpy.float32 ),
           [ numpy.array( [ 4, 5, 3, 1 ], numpy.float32 ),
             numpy.array( [ 2, 3, 5, 6 ], numpy.float32 ) ] )

        b = Observations(
           numpy.array( [ False, False, False, False, False ], numpy.bool ),
           numpy.float32( 52000.0 ),
           numpy.array( [   0,  1,    2,   3,   4 ], numpy.uint64),
           numpy.array( [ 280, 280, 300, 310, 320 ], numpy.float32 ),
           numpy.array( [   4,   6,   3,   2,   4 ], numpy.float32 ),
           [ numpy.array( [ 3,   5,   4,   1,   3 ], numpy.float32 ),
             numpy.array( [ 8,   3,  10,   3,   1 ], numpy.float32 ) ] )

        target = Observations(
           numpy.array( [ False, True, False ], numpy.bool ),
           numpy.float32( 52000.0 ),
           numpy.array( [  0   ,   1  ,  3   ], numpy.uint64),
           numpy.array( [ 275, 280, 300 ], numpy.float32 ),
           numpy.array(   [ numpy.sqrt((3.**2+4.**2)/2.), numpy.sqrt((6.**2+6.**2)/2.), numpy.sqrt((4.**2+2.**2)/2.) ], numpy.float32 ),
           [ numpy.array( [ numpy.sqrt((4.**2+3.**2)/2.), numpy.sqrt((5.**2+5.**2)/2.), numpy.sqrt((3.**2+1.**2)/2.) ], numpy.float32 ),
             numpy.array( [ numpy.sqrt((2.**2+8.**2)/2.), numpy.sqrt((3.**2+3.**2)/2.), numpy.sqrt((5.**2+3.**2)/2.) ], numpy.float32 ) ] )

        c = Observations.mean(a, b)
        print c.location
        #numpy.testing.assert_equal(target.mask, c.mask)
        numpy.testing.assert_equal(target.time, c.time)
        numpy.testing.assert_equal(target.location, c.location)
        numpy.testing.assert_almost_equal(target.measurement, c.measurement)
        numpy.testing.assert_almost_equal(target.uncorrelatederror, c.uncorrelatederror)
        numpy.testing.assert_almost_equal(target.locallycorrelatederror, c.locallycorrelatederror)

    def test_append(self):

        obs = Observations(           
           numpy.array( [ False, True, False, False ], numpy.bool ),
           numpy.float32( 51000.0 ),
           numpy.array( [  1   ,   3  ,    0   ,   2   ], numpy.uint64),
           numpy.array( [ 272.3, 278.4, 271.111, 288.0 ], numpy.float32 ),
           numpy.array( [ 8.0, 9.33, 0.1, 2.2222 ], numpy.float32 ),
           [ numpy.array( [ 1.0, 3.3, 0.2, 7.0 ], numpy.float32 ),
             numpy.array( [ 2.0, 0.1, 0.0, 1.0 ], numpy.float32 ) ] )

        second = Observations(           
           numpy.array( [ True, False ], numpy.bool ),
           numpy.float32( 51000.0 ),
           numpy.array( [  0, 1   ], numpy.uint64),
           numpy.array( [ 209.0, 210.0 ], numpy.float32 ),
           numpy.array( [ 1.3, 2.4 ], numpy.float32 ),
           [ numpy.array( [ 1.11, 2.22 ], numpy.float32 ),
             numpy.array( [ 3.33, 4.44 ], numpy.float32 ) ] )

        obs.append(second, locationoffset=1000)

        numpy.testing.assert_equal([ False, True, False, False, True, False ], obs.mask)
        self.assertEqual(numpy.float32( 51000.0 ), obs.time)
        numpy.testing.assert_equal([  1   ,   3  ,    0   ,   2   , 1000   , 1001   ], obs.location)
        numpy.testing.assert_almost_equal([ 272.3, 278.4, 271.111, 288.0, 209.0, 210.0 ], obs.measurement, decimal=4)
        numpy.testing.assert_almost_equal([ 8.0, 9.33, 0.1, 2.2222, 1.3, 2.4 ], obs.uncorrelatederror, decimal=4)
        self.assertEqual(2, len(obs.locallycorrelatederror))
        numpy.testing.assert_almost_equal([ 1.0, 3.3, 0.2, 7.0, 1.11, 2.22 ], obs.locallycorrelatederror[0], decimal=4)
        numpy.testing.assert_almost_equal([ 2.0, 0.1, 0.0, 1.0, 3.33, 4.44 ], obs.locallycorrelatederror[1], decimal=4)

class TestObservationsBreakPoints(unittest.TestCase):

    def test_init_and_count(self):
	
	obs = ObservationsBreakPoints(numpy.array([54786, 56978, 57343, 57708, 58074], dtype=numpy.int32), 
					       numpy.array([4, 4, 4, 4, 4], dtype=numpy.int32),
					       numpy.array([1, 1, 1, 1, 1], dtype=numpy.int8),)

	numpy.testing.assert_array_equal([54786, 56978, 57343, 57708, 58074], obs.break_time,err_msg='Testing initialization \"break_time\" attribute')
	numpy.testing.assert_array_equal([4, 4, 4, 4, 4], obs.break_station,err_msg='Testing initialization \"break_station\" attribute')
	numpy.testing.assert_array_equal([1, 1, 1, 1, 1], obs.break_likelihood,err_msg='Testing initialization of \"break_likelihood\" attribute')
	self.assertEqual(None, obs.break_type,'Testing initialization of \"break_type\" attribute')
	self.assertEqual(None, obs.detection_feasibility,'Testing initialization of \"detection_feasibility\" attribute')
	self.assertEqual(None, obs.detection_score,'Testing initialization of \"detection_score\" attribute')
	self.assertEqual(5, obs.number_of_observations())
	self.assertItemsEqual(obs.POLICY_CODES,['HARD_CUTOFF', 'LAPLACE_KERNEL'])

    def test_station_related_methods(self):

	obs = ObservationsBreakPoints(numpy.array([54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074], dtype=numpy.int32), 
					       numpy.array([4, 3, 3, 3, 2, 0, 0, 0], dtype=numpy.int32),
					       numpy.array([1, 1, 1, 1, 1, 1, 1, 1], dtype=numpy.int8),)
	
	self.assertEqual(4,obs.total_number_of_stations())
	numpy.testing.assert_array_equal(numpy.array([0,2,3,4]), obs.stations_collection())
	numpy.testing.assert_array_equal(numpy.array([3,1,3,1]), obs.stations_count())
      
	test_arrays = [numpy.array([5, 6, 7]), numpy.array([4]), numpy.array([1, 2, 3]), numpy.array([0])]
	test_values = [numpy.array([57343, 57708, 58074]), numpy.array([58074]), numpy.array([56978, 57343, 57708]), numpy.array([54786])]
	for index, item in enumerate(zip(test_arrays, test_values)):
	    numpy.testing.assert_array_equal(item[0], obs.stations_break_point_indices([obs.stations_collection()[index]]))
	    numpy.testing.assert_array_equal(item[1], obs.filtered_break_points([obs.stations_collection()[index]]))

    def test_breaking_points_conversion(self):

	obs = ObservationsBreakPoints(numpy.array([54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074], dtype=numpy.int32), 
					       numpy.array([4, 3, 3, 3, 2, 0, 0, 0], dtype=numpy.int32),
					       numpy.array([1, 1, 1, 1, 1, 1, 1, 1], dtype=numpy.int8),)
	
	test_array =  numpy.array(['2000-01-01', '2006-01-01', '2007-01-01', '2008-01-01', '2009-01-01', '2007-01-01', '2008-01-01', '2009-01-01'])
	numpy.testing.assert_array_equal(obs.break_points_timestamp(datetime(1850,1,1)), test_array, err_msg = 'Testing breaking points convertion into datetime array')

    def test_apply_policy(self):
	
	# Test value error raising and HARD CUT OFF policy
	obs = ObservationsBreakPoints(numpy.array([54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074], dtype=numpy.int32), 
				      numpy.array([4, 3, 3, 3, 2, 0, 0, 0], dtype=numpy.int32),
				      numpy.array([1, 100, 112, 13, 10, 21, 2, 3], dtype=numpy.int8),
				      detection_feasibility = numpy.array([0, 1, 2, 0, 1, 0]))
	self.assertRaises(ValueError, obs.apply_policy,'brombobom')

	obs.apply_policy('HARD_CUTOFF', threshold = 13)
	numpy.testing.assert_array_equal([56978, 57343, 57708, 57343], obs.break_time,err_msg='Testing \"HARD_CUTOFF\" policy application on \"break_time\" attribute')
	numpy.testing.assert_array_equal([3, 3, 3, 0], obs.break_station,err_msg='Testing \"HARD_CUTOFF\" policy application on \"break_station\" attribute')
	numpy.testing.assert_array_equal([100, 112, 13, 21], obs.break_likelihood,err_msg='Testing \"HARD_CUTOFF\" policy application on  \"break_likelihood\" attribute')

	# Test LAPLACE KERNEL policy
	obs = ObservationsBreakPoints(numpy.array([54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074], dtype=numpy.int32), 
				      numpy.array([4, 3, 3, 3, 2, 0, 0, 0], dtype=numpy.int32),
				      numpy.array([1, 100, 112, 13, 10, 21, 2, 3], dtype=numpy.int8),
				      detection_feasibility = numpy.array([0, 1, 2, 1, 1, 1]))

	obs.apply_policy('LAPLACE_KERNEL',decay_constant=.01,threshold=.5)
	numpy.testing.assert_array_equal([56978, 57343], obs.break_time,err_msg='Testing \"LAPLACE_KERNEL\" policy application (standard threshold) on \"break_time\" attribute')
	numpy.testing.assert_array_equal([3, 3], obs.break_station,err_msg='Testing \"LAPLACE_KERNEL\" policy application (standard threshold) on \"break_station\" attribute')
	numpy.testing.assert_array_equal([100, 112], obs.break_likelihood,err_msg='Testing \"LAPLACE_KERNEL\" policy application (standard threshold) on  \"break_likelihood\" attribute')

	obs = ObservationsBreakPoints(numpy.array([54786, 56978, 57343, 57708, 58074, 57343, 57708, 58074], dtype=numpy.int32), 
				      numpy.array([4, 3, 3, 3, 2, 0, 0, 0], dtype=numpy.int32),
				      numpy.array([115, 100, 112, 13, 10, 21, 125, 3], dtype=numpy.int8),
				      detection_feasibility = numpy.array([0, 1, 2, 1, 1, 1]))

	obs.apply_policy('LAPLACE_KERNEL',decay_constant=.02,threshold=.4)
	numpy.testing.assert_array_equal([54786, 56978, 57343, 57708], obs.break_time,err_msg='Testing \"LAPLACE_KERNEL\" policy application (user defined threshold) on \"break_time\" attribute')
	numpy.testing.assert_array_equal([4, 3, 3, 0], obs.break_station,err_msg='Testing \"LAPLACE_KERNEL\" policy application (user defined threshold) on \"break_station\" attribute')
	numpy.testing.assert_array_equal([115, 100, 112, 125], obs.break_likelihood,err_msg='Testing \"LAPLACE_KERNEL\" policy application (user defined threshold) on  \"break_likelihood\" attribute')

class TestObservationSource(unittest.TestCase):

    def test_init(self):

        s = ObservationSource()
        self.assertRaises(NotImplementedError, ObservationSource.observation_location_lookup, s)
        self.assertRaises(NotImplementedError, ObservationSource.observations, s, ObservationSource.TMEAN)
        self.assertRaises(NotImplementedError, ObservationSource.observables, s)
        self.assertRaises(NotImplementedError, ObservationSource.local_correlation_length_scale, s, ObservationSource.TMEAN)
