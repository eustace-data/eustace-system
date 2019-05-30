"""Test interface to insitu land breaking points preprocessing."""

import unittest
import eustaceconfig
import numpy
import os
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import ObservationsBreakPoints
from eustace.outputformats import definitions
from ..insitu_land_breakpoints import ObservationBreakPointSourceInsituLand

class TestObservationBreakPointSourceInsituLand(unittest.TestCase):

    def test_load_example_D_1_7(self):

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001127_daily_status.nc')

        # attempt load: no station specified
        result = ObservationBreakPointSourceInsituLand(FILENAME_EXAMPLE)

        # check keys correctness
        self.assertItemsEqual(['merged_break'  ], result.observables())

        sub_keys = ['dimension','break_time', 'break_station', 'break_likelihood', 'detection_feasibility']
        self.assertItemsEqual(sub_keys, result.OBSERVATIONMAPS['merged_break'].keys())

        # get number of observations 
        keys = ['merged_break']

        expected_count = [130950]
        for observable, value in zip(keys,expected_count):
            self.assertEqual(result.number_of_observations(observable),value,'Testing number of observed breaking points for observable '+observable)

        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 2 rows with a column per station location
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, 35375), location_lookup.shape)

        # check coordinates
        numpy.testing.assert_almost_equal([ 48.4000, -123.4833 ], location_lookup[:,0].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([ 48.8167, -124.1333 ], location_lookup[:,10].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([-22.2170,   30.0000 ], location_lookup[:,35374].transpose(), decimal=4)

        # check shape and type of observations (merged_break)
        merged_break = result.observations('merged_break')
        self.assertTrue(isinstance(merged_break, ObservationsBreakPoints))

	observables = [merged_break.break_time, merged_break.break_station, merged_break.break_likelihood]
	array_types = [numpy.int32, numpy.int32, numpy.int8]
	array_shapes = [(130950,), (130950,), (130950,)]
	array_indices = [24, 119]
	array_values = [[47116, 32871],[9, 43],[4, 4]]

	for obs, array_type, array_shape, array_value in zip(observables, array_types, array_shapes, array_values):
            self.assertTrue(isinstance(obs, numpy.ndarray))
	    self.assertTrue(obs.dtype == array_type)
            self.assertEqual(array_shape, obs.shape)
	    self.assertEqual(array_value[0], obs[array_indices[0]])
	    self.assertEqual(array_value[1], obs[array_indices[1]])

	obs = merged_break.detection_feasibility
	self.assertTrue(isinstance(obs, numpy.ndarray))
	self.assertTrue(obs.dtype == numpy.int8)
	self.assertEqual((35375,), obs.shape)
	self.assertEqual(1, obs[24])
	self.assertEqual(1, obs[19])

    def test_load_example_D_1_7_single_station(self):

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001127_daily_status.nc')

        # attempt load: we specify a station.
        result = ObservationBreakPointSourceInsituLand(FILENAME_EXAMPLE, 15)

        # check keys correctness
        self.assertItemsEqual(['merged_break'  ], result.observables())

        # check sub-keys correctness
        sub_keys = ['dimension','break_time', 'break_station', 'break_likelihood', 'detection_feasibility']
        self.assertItemsEqual(sub_keys, result.OBSERVATIONMAPS['merged_break'].keys())

        # get number of observations 
        keys = ['merged_break']

        expected_count = [4]
        for observable, value in zip(keys, expected_count):
            self.assertEqual(result.number_of_observations(observable),value,'Testing number of observed breaking points for observable '+observable)

        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 2 rows with a column per station location
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, 35375), location_lookup.shape)

        # check coordinates
        numpy.testing.assert_almost_equal([ 48.4000, -123.4833 ], location_lookup[:,0].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([ 48.8167, -124.1333 ], location_lookup[:,10].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([-22.2170,   30.0000 ], location_lookup[:,35374].transpose(), decimal=4)
        
        # check shape and type of observations (merged_break)
        merged_break = result.observations('merged_break')
        self.assertTrue(isinstance(merged_break, ObservationsBreakPoints))

	observables = [merged_break.break_time, merged_break.break_station, merged_break.break_likelihood]
	array_types = [numpy.int32, numpy.int32, numpy.int8]
	array_shapes = [(4,), (4,), (4,)]
	array_indices = [0, 3]
	array_values = [[30315, 37620],[15, 15],[14, 18]]

	for obs, array_type, array_shape, array_value in zip(observables, array_types, array_shapes, array_values):
            self.assertTrue(isinstance(obs, numpy.ndarray))
	    self.assertTrue(obs.dtype == array_type)
            self.assertEqual(array_shape, obs.shape)
	    self.assertEqual(array_value[0], obs[array_indices[0]])
	    self.assertEqual(array_value[1], obs[array_indices[1]])

	obs = merged_break.detection_feasibility
	self.assertTrue(isinstance(obs, numpy.ndarray))
	self.assertTrue(obs.dtype == numpy.int8)
	self.assertEqual((1,), obs.shape)
	self.assertEqual(0, obs[0])
	      