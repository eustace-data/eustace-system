"""Test interface to insitu land preprocessing."""

import unittest
import numpy
import eustaceconfig
import os
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from ..insitu_land import ObservationSourceInsituLand

class MockAdjustment():
    """Mockup of AltitudeAdjustment class for known adjustment values at all locations"""
    
    def __init__(self, dummy_adjustment, dummy_adjustment_uncertainty):
        
        self.dummy_adjustment = dummy_adjustment
        self.dummy_adjustment_uncertainty = dummy_adjustment_uncertainty
    
    def load_reference_map(self):
        pass
    
    def adjustment(self, altitudes, locations):
        return numpy.ones(locations.shape[0]) * self.dummy_adjustment
    
    def adjustment_uncertainty(self, altitudes, locations):
        return numpy.ones(locations.shape[0]) * self.dummy_adjustment_uncertainty

class TestObservationSourceInsituLand(unittest.TestCase):
    
    def atest_load_example_D_1_7(self):

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')

        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE)
	self.assertFalse(result.hold_out_dataset)
        # only supports tmean
        self.assertEqual([ ObservationSource.TMAX, ObservationSource.TMIN  ], result.observables())

        # get number of observations [ = total stations x 365 days ]
        # (some will be masked out)
        self.assertEqual(12911875, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 2 rows with a column per station location
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, 35375), location_lookup.shape)

        # check coordinates
        numpy.testing.assert_almost_equal([ 48.4000, -123.4833 ], location_lookup[:,0].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([ 48.5000, -124.0000 ], location_lookup[:,1].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([-22.2170,   30.0000 ], location_lookup[:,35374].transpose(), decimal=4)

        # check shape and type of observations (TMIN)
        tmin = result.observations(ObservationSource.TMIN, False)
        self.assertTrue(isinstance(tmin, Observations))
        self.assertTrue(isinstance(tmin.mask, numpy.ndarray))
        self.assertTrue(tmin.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmin.measurement, numpy.ndarray))
        self.assertTrue(tmin.measurement.dtype == numpy.float32)
        self.assertEqual((12911875,), tmin.measurement.shape)
        self.assertTrue(isinstance(tmin.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmin.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((12911875,), tmin.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmin.locallycorrelatederror, list))
        self.assertEqual(0, len(tmin.locallycorrelatederror))

        # check first few values
        self.assertTrue(tmin.mask[0])
        self.assertTrue(tmin.mask[1])
        self.assertTrue(tmin.mask[2])
        self.assertFalse(tmin.mask[3])
        self.assertAlmostEqual(277.15, tmin.measurement[3], places=4)
        self.assertTrue(tmin.mask[4])
        self.assertTrue(tmin.mask[5])
        self.assertTrue(tmin.mask[6])
        self.assertTrue(tmin.mask[7])
        self.assertTrue(tmin.mask[8])
        self.assertTrue(tmin.mask[9])
        self.assertFalse(tmin.mask[10])
        self.assertAlmostEqual(276.15, tmin.measurement[10], places=4)

        # at present the uncorrelated error is set to 5 everywhere
        self.assertTrue(numpy.all(tmin.uncorrelatederror == numpy.float32(0.3)))

        # check shape and type of observations (TMAX)
        tmax = result.observations(ObservationSource.TMAX, False)
        self.assertTrue(isinstance(tmax, Observations))
        self.assertTrue(isinstance(tmax.mask, numpy.ndarray))
        self.assertTrue(tmax.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmax.measurement, numpy.ndarray))
        self.assertTrue(tmax.measurement.dtype == numpy.float32)
        self.assertEqual((12911875,), tmax.measurement.shape)
        self.assertTrue(isinstance(tmax.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmax.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((12911875,), tmax.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmax.locallycorrelatederror, list))
        self.assertEqual(0, len(tmax.locallycorrelatederror))

        # check first few values
        self.assertTrue(tmax.mask[0])
        self.assertTrue(tmax.mask[1])
        self.assertTrue(tmax.mask[2])
        self.assertFalse(tmax.mask[3])
        self.assertAlmostEqual(280.65, tmax.measurement[3], places=4)
        self.assertTrue(tmin.mask[4])
        self.assertTrue(tmax.mask[5])
        self.assertTrue(tmax.mask[6])
        self.assertTrue(tmax.mask[7])
        self.assertTrue(tmax.mask[8])
        self.assertTrue(tmax.mask[9])
        self.assertFalse(tmax.mask[10])
        self.assertAlmostEqual(283.65, tmax.measurement[10], places=4)
        self.assertTrue(tmax.mask[11])
        self.assertFalse(tmax.mask[12])
        self.assertAlmostEqual(283.15, tmax.measurement[12], places=4)

        # at present the uncorrelated error is set to 5 everywhere
        self.assertTrue(numpy.all(tmax.uncorrelatederror == numpy.float32(0.3)))

        # check empty array of correlation length scales
        tmin_length_scale = result.local_correlation_length_scale(ObservationSource.TMIN)
        self.assertTrue(isinstance(tmin_length_scale, numpy.ndarray))
        self.assertEqual(0, len(tmin_length_scale))
        tmax_length_scale = result.local_correlation_length_scale(ObservationSource.TMAX)
        self.assertTrue(isinstance(tmax_length_scale, numpy.ndarray))
        self.assertEqual(0, len(tmax_length_scale))


    def atest_load_example_D_1_7_single_day(self):

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')

        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE, daynumber=57344)
	self.assertFalse(result.hold_out_dataset)
        # only supports tmean
        self.assertEqual([ ObservationSource.TMAX, ObservationSource.TMIN  ], result.observables())

        # get number of observations [ = total stations x 1 day ]
        # (some will be masked out)
        self.assertEqual(35375, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 2 rows with a column per station location
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, 35375), location_lookup.shape)

        # check coordinates
        numpy.testing.assert_almost_equal([ 48.4000, -123.4833 ], location_lookup[:,0].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([ 48.5000, -124.0000 ], location_lookup[:,1].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([-22.2170,   30.0000 ], location_lookup[:,35374].transpose(), decimal=4)

        # check shape and type of observations (TMIN)
        tmin = result.observations(ObservationSource.TMIN, False)
        self.assertTrue(isinstance(tmin, Observations))
        self.assertTrue(isinstance(tmin.mask, numpy.ndarray))
        self.assertTrue(tmin.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmin.measurement, numpy.ndarray))
        self.assertTrue(tmin.measurement.dtype == numpy.float32)
        self.assertEqual((35375,), tmin.measurement.shape)
        self.assertTrue(isinstance(tmin.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmin.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((35375,), tmin.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmin.locallycorrelatederror, list))
        self.assertEqual(0, len(tmin.locallycorrelatederror))

        # check first few values
        self.assertTrue(tmin.mask[0])
        self.assertTrue(tmin.mask[1])
        self.assertTrue(tmin.mask[2])
        self.assertFalse(tmin.mask[3])
        self.assertAlmostEqual(279.65, tmin.measurement[3], places=4)
        self.assertTrue(tmin.mask[4])

        # at present the uncorrelated error is set to 5 everywhere
        self.assertTrue(numpy.all(tmin.uncorrelatederror == numpy.float32(0.3)))

        # check shape and type of observations (TMAX)
        tmax = result.observations(ObservationSource.TMAX, False)
        self.assertTrue(isinstance(tmax, Observations))
        self.assertTrue(isinstance(tmax.mask, numpy.ndarray))
        self.assertTrue(tmax.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmax.measurement, numpy.ndarray))
        self.assertTrue(tmax.measurement.dtype == numpy.float32)
        self.assertEqual((35375,), tmax.measurement.shape)
        self.assertTrue(isinstance(tmax.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmax.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((35375,), tmax.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmax.locallycorrelatederror, list))
        self.assertEqual(0, len(tmax.locallycorrelatederror))

        # check first few values
        self.assertTrue(tmax.mask[0])
        self.assertTrue(tmax.mask[1])
        self.assertTrue(tmax.mask[2])
        self.assertFalse(tmax.mask[3])
        self.assertAlmostEqual(285.65, tmax.measurement[3], places=4)
        self.assertTrue(tmin.mask[4])

        # at present the uncorrelated error is set to 5 everywhere
        self.assertTrue(numpy.all(tmax.uncorrelatederror == numpy.float32(0.3)))

        # check empty array of correlation length scales
        tmin_length_scale = result.local_correlation_length_scale(ObservationSource.TMIN)
        self.assertTrue(isinstance(tmin_length_scale, numpy.ndarray))
        self.assertEqual(0, len(tmin_length_scale))
        tmax_length_scale = result.local_correlation_length_scale(ObservationSource.TMAX)
        self.assertTrue(isinstance(tmax_length_scale, numpy.ndarray))
        self.assertEqual(0, len(tmax_length_scale))
            

    def atest_load_example_D_1_7_hold_out_data(self):

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')
        FILENAME_HOLD_OUT = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/M27/station_selection/Land/Bern_WP3_reserved_stations_land_matched2R001127.nc')

        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE, hold_out_filename=FILENAME_HOLD_OUT)

        # only supports tmean
        self.assertEqual([ ObservationSource.TMAX, ObservationSource.TMIN  ], result.observables())

        # get number of observations [ = total stations x 365 days ]
        # (some will be masked out)
        self.assertEqual(12911875, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 2 rows with a column per station location
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, 35375), location_lookup.shape)

        # check coordinates
        numpy.testing.assert_almost_equal([ 48.4000, -123.4833 ], location_lookup[:,0].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([ 48.5000, -124.0000 ], location_lookup[:,1].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([-22.2170,   30.0000 ], location_lookup[:,35374].transpose(), decimal=4)

        # check shape and type of observations (TMIN)
        tmin = result.observations(ObservationSource.TMIN, False)
        self.assertTrue(isinstance(tmin, Observations))
        self.assertTrue(isinstance(tmin.mask, numpy.ndarray))
        self.assertTrue(tmin.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmin.measurement, numpy.ndarray))
        self.assertTrue(tmin.measurement.dtype == numpy.float32)
        self.assertEqual((12911875,), tmin.measurement.shape)
        self.assertTrue(isinstance(tmin.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmin.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((12911875,), tmin.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmin.locallycorrelatederror, list))
        self.assertEqual(0, len(tmin.locallycorrelatederror))

        # check first few values
        self.assertTrue(tmin.mask[13])
        self.assertTrue(tmin.mask[23])
        self.assertTrue(tmin.mask[24])
        self.assertFalse(tmin.mask[3])
        self.assertAlmostEqual(277.15, tmin.measurement[3], places=4)
        self.assertTrue(tmin.mask[35351])
        self.assertTrue(tmin.mask[35352])
        self.assertTrue(tmin.mask[35354])
        self.assertTrue(tmin.mask[30])
        self.assertTrue(tmin.mask[35])
        self.assertTrue(tmin.mask[99])
        self.assertFalse(tmin.mask[10])
        self.assertAlmostEqual(276.15, tmin.measurement[10], places=4)

        # at present the uncorrelated error is set to 5 everywhere
        self.assertTrue(numpy.all(tmin.uncorrelatederror == numpy.float32(0.3)))

        # check shape and type of observations (TMAX)
        tmax = result.observations(ObservationSource.TMAX, False)
        self.assertTrue(isinstance(tmax, Observations))
        self.assertTrue(isinstance(tmax.mask, numpy.ndarray))
        self.assertTrue(tmax.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmax.measurement, numpy.ndarray))
        self.assertTrue(tmax.measurement.dtype == numpy.float32)
        self.assertEqual((12911875,), tmax.measurement.shape)
        self.assertTrue(isinstance(tmax.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmax.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((12911875,), tmax.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmax.locallycorrelatederror, list))
        self.assertEqual(0, len(tmax.locallycorrelatederror))

        # check first few values
        self.assertTrue(tmax.mask[157])
        self.assertTrue(tmax.mask[172])
        self.assertTrue(tmax.mask[225])
        self.assertFalse(tmax.mask[3])
        self.assertAlmostEqual(280.65, tmax.measurement[3], places=4)
        self.assertTrue(tmin.mask[400])
        self.assertTrue(tmax.mask[505])
        self.assertTrue(tmax.mask[535])
        self.assertTrue(tmax.mask[622])
        self.assertTrue(tmax.mask[832])
        self.assertTrue(tmax.mask[883])
        self.assertFalse(tmax.mask[10])
        self.assertAlmostEqual(283.65, tmax.measurement[10], places=4)
        self.assertTrue(tmax.mask[11])
        self.assertFalse(tmax.mask[12])
        self.assertAlmostEqual(283.15, tmax.measurement[12], places=4)

        # at present the uncorrelated error is set to 5 everywhere
        self.assertTrue(numpy.all(tmax.uncorrelatederror == numpy.float32(0.3)))

        # check empty array of correlation length scales
        tmin_length_scale = result.local_correlation_length_scale(ObservationSource.TMIN)
        self.assertTrue(isinstance(tmin_length_scale, numpy.ndarray))
        self.assertEqual(0, len(tmin_length_scale))
        tmax_length_scale = result.local_correlation_length_scale(ObservationSource.TMAX)
        self.assertTrue(isinstance(tmax_length_scale, numpy.ndarray))
        self.assertEqual(0, len(tmax_length_scale))
       

    def atest_D_1_7_qc_data_content(self):
        """ This test compares the number of observation before and after qc flags
            are applied, to verify the number of observations reduces 
            after the qc flags have been applied
        """

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')
        FILENAME_HOLD_OUT = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/M27/station_selection/Land/Bern_WP3_reserved_stations_land_matched2R001127.nc')
        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE, hold_out_filename=FILENAME_HOLD_OUT)

        # get number of observations [ = total stations x 365 days ]
        # (some will be masked out)
        self.assertEqual(12911875, result.number_of_observations())

        # access the object - not invoking qc flag use
        tmin = result.observations(ObservationSource.TMIN, False)
        tmax = result.observations(ObservationSource.TMAX, False)

        # get the number of days and total number of stations, total observations
        days = result.number_of_days()
        stations = result.number_of_stations()
        total_obs = result.number_of_observations()

        # calculate total number of actual observations by ignoring fill values
        total_actual_obs_tmin = (numpy.count_nonzero(numpy.where(tmin.measurement > -32768)))
        total_actual_obs_tmax = (numpy.count_nonzero(numpy.where(tmax.measurement > -32768)))

        # total number of valid values - non masked
        tmin_values = numpy.sum(tmin.mask == 0)
        tmax_values = numpy.sum(tmax.mask == 0)

        # calculate valid values as % of the total number of actual observations 
        pcnt_tmin = (tmin_values / float(total_actual_obs_tmin))*100 
        pcnt_tmax = (tmax_values / float(total_actual_obs_tmax))*100

        print(" --- ")
        print(" --- Observations - validation holdout only --- ")   
        print(" --- ")
        print("Observations - days x stations:  {} {}".format(days, stations))
        print("Observations - measurement shape: {}".format(tmin.measurement.shape))
        print("Observations - total potential obs:  {}".format(total_obs))
        print("Observations - total actual obs: tmin  {}".format(total_actual_obs_tmin))
        print("Observations - total actual obs: tmax  {}".format(total_actual_obs_tmax))
        print("Observations - Tmin used values: {} {}% of total actual obs".format(tmin_values, pcnt_tmin))
        print("Observations - Tmax used values: {} {}% of total actual obs".format(tmax_values, pcnt_tmax))
        print(" --- ")

        # access the object - invoking qc flag use - so more values are masked out 
        tmin_qc = result.observations(ObservationSource.TMIN, True)
        tmax_qc = result.observations(ObservationSource.TMAX, True)

        # total number of valid values - non masked
        tmin_values_qc = numpy.sum(tmin_qc.mask == 0)
        tmax_values_qc = numpy.sum(tmax_qc.mask == 0)

        # calculate valid values as % of the total number of actual observations 
        pcnt_tmin_qc = (tmin_values_qc / float(total_actual_obs_tmin))*100
        pcnt_tmax_qc = (tmax_values_qc / float(total_actual_obs_tmax))*100

        print(" --- ")
        print(" --- Observations - validation holdout and qc flags --- ")
        print(" --- ")
        print("Observations - days x stations:  {} {}".format(days, stations))
        print("Observations - measurement shape: {}".format(tmin.measurement.shape))
        print("Observations - total potential obs:  {}".format(total_obs))
        print("Observations - total actual obs: tmin  {}".format(total_actual_obs_tmin))
        print("Observations - total actual obs: tmax  {}".format(total_actual_obs_tmax))
        print("Observations - Tmin used values qc: {} {}% of total actual obs".format(tmin_values_qc, pcnt_tmin_qc))
        print("Observations - Tmax used values qc: {} {}% of total actual obs".format(tmax_values_qc, pcnt_tmax_qc))
        print(" --- ")

        # test that there are fewer valid observations once qc flags have been applied
        self.assertLess(tmin_values_qc, total_actual_obs_tmin)
        self.assertLess(tmax_values_qc, total_actual_obs_tmax)

        # test the number of valid, available observations after holdout and qc flag masks
        self.assertEqual(4508706, numpy.sum(tmin_qc.mask ==0))
        self.assertEqual(4475502, numpy.sum(tmax_qc.mask ==0))

        # test there are more masked values (True) after applying the qc flag masks
        self.assertGreater(numpy.sum(tmin_qc.mask), numpy.sum(tmin.mask))
        self.assertGreater(numpy.sum(tmax_qc.mask), numpy.sum(tmax.mask))


    def atest_load_example_D_1_7_single_day_with_qc_flags(self):
        """ This test compares the number of observation before and after qc flags
            are applied, to verify the number of observations reduces  
            after the qc flags have been applied - for a single day observation value
        """

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')

        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE, daynumber=57344)
        self.assertFalse(result.hold_out_dataset)

        # get number of observations [ = total stations x 1 day ]
        # (some will be masked out)
        # the value will be 35375 stations
        self.assertEqual(35375, result.number_of_observations())
        
        # access the object - not invoking qc flag use
        tmin = result.observations(ObservationSource.TMIN, False)
        tmax = result.observations(ObservationSource.TMAX, False)

        # get the number of days and total number of stations, total observations
        days = result.number_of_days()
        stations = result.number_of_stations()
        total_obs = result.number_of_observations()

        # calculate total number of actual observations by ignoring fill values
        total_actual_obs_tmin = (numpy.count_nonzero(numpy.where(tmin.measurement > -32768)))
        total_actual_obs_tmax = (numpy.count_nonzero(numpy.where(tmax.measurement > -32768)))

        # total number of valid values - non masked
        tmin_values = numpy.sum(tmin.mask == 0)
        tmax_values = numpy.sum(tmax.mask == 0)

        # calculate valid values as % of the total number of actual observations 
        pcnt_tmin = (tmin_values / float(total_actual_obs_tmin))*100
        pcnt_tmax = (tmax_values / float(total_actual_obs_tmax))*100


        print(" --- ")
        print(" --- Observations - single day - no qc flag use --- ")
        print(" --- ")
        print("Observations - days x stations:  {} {}".format(days, stations))
        print("Observations - measurement shape: {}".format(tmin.measurement.shape))
        print("Observations - total potential obs:  {}".format(total_obs))
        print("Observations - total actual obs: tmin  {}".format(total_actual_obs_tmin))
        print("Observations - total actual obs: tmax  {}".format(total_actual_obs_tmax))
        print("Observations - Tmin used values: {} {}% of total actual obs".format(tmin_values, pcnt_tmin))
        print("Observations - Tmax used values: {} {}% of total actual obs".format(tmax_values, pcnt_tmax))
        print(" --- ")

        # access the object - invoking qc flag use - so more values are masked out 
        tmin_qc = result.observations(ObservationSource.TMIN, True)
        tmax_qc = result.observations(ObservationSource.TMAX, True)

        # total number of valid values - non masked
        tmin_values_qc = numpy.sum(tmin_qc.mask == 0)
        tmax_values_qc = numpy.sum(tmax_qc.mask == 0)

        # calculate valid values as % of the total number of actual observations 
        pcnt_tmin_qc = (tmin_values_qc / float(total_actual_obs_tmin))*100
        pcnt_tmax_qc = (tmax_values_qc / float(total_actual_obs_tmax))*100

        print(" --- ")
        print(" --- Observations - single day - with qc flags --- ")
        print(" --- ")
        print("Observations - days x stations:  {} {}".format(days, stations))
        print("Observations - measurement shape: {}".format(tmin.measurement.shape))
        print("Observations - total potential obs:  {}".format(total_obs))
        print("Observations - total actual obs: tmin  {}".format(total_actual_obs_tmin))
        print("Observations - total actual obs: tmax  {}".format(total_actual_obs_tmax))
        print("Observations - Tmin used values qc: {} {}% of total actual obs".format(tmin_values_qc, pcnt_tmin_qc))
        print("Observations - Tmax used values qc: {} {}% of total actual obs".format(tmax_values_qc, pcnt_tmax_qc))
        print(" --- ")

       
        # test there are fewer valid observations once qc flags have been applied
        self.assertLess(tmin_values_qc, total_actual_obs_tmin)
        self.assertLess(tmax_values_qc, total_actual_obs_tmax)

        # test the number of valid, available observations after holdout and qc flag masks
        self.assertEqual(13692, numpy.sum(tmin_qc.mask ==0))
        self.assertEqual(13625, numpy.sum(tmax_qc.mask ==0))

        # test there are more masked values (True) after applying the qc flag masks
        self.assertGreater(numpy.sum(tmin_qc.mask), numpy.sum(tmin.mask))
        self.assertGreater(numpy.sum(tmax_qc.mask), numpy.sum(tmax.mask))

   
    def atest_load_example_D_1_7_single_day_with_qc_flags_and_holdout_file(self):
        """ This test compares the number of observation before and after qc flags
            and holdout files are applied, to verify the number of observations reduces, 
            after the qc flags have been applied - for a single day observation value
        """

        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')

        FILENAME_HOLD_OUT = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/M27/station_selection/Land/Bern_WP3_reserved_stations_land_matched2R001127.nc')

        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE, daynumber=57344)
        self.assertFalse(result.hold_out_dataset)

        # get number of observations [ = total stations x 1 day ]
        # (some will be masked out)
        # the value will be 35375 stations
        self.assertEqual(35375, result.number_of_observations())

        # access the object - not invoking qc flag use
        tmin = result.observations(ObservationSource.TMIN, False)
        tmax = result.observations(ObservationSource.TMAX, False)

        # get the number of days and total number of stations, total observations
        days = result.number_of_days()
        stations = result.number_of_stations()
        total_obs = result.number_of_observations()

        # calculate total number of actual observations by ignoring fill values
        total_actual_obs_tmin = (numpy.count_nonzero(numpy.where(tmin.measurement > -32768)))
        total_actual_obs_tmax = (numpy.count_nonzero(numpy.where(tmax.measurement > -32768)))

        # total number of valid values - non masked
        tmin_values = numpy.sum(tmin.mask == 0)
        tmax_values = numpy.sum(tmax.mask == 0)

        # calculate valid values as % of the total number of actual observations 
        pcnt_tmin = (tmin_values / float(total_actual_obs_tmin))*100
        pcnt_tmax = (tmax_values / float(total_actual_obs_tmax))*100

        print(" --- ")
        print(" --- Observations - single day - no qc flag use - no holdout file --- ")
        print(" --- ")
        print("Observations - days x stations:  {} {}".format(days, stations))
        print("Observations - measurement shape: {}".format(tmin.measurement.shape))
        print("Observations - total potential obs:  {}".format(total_obs))
        print("Observations - total actual obs: tmin  {}".format(total_actual_obs_tmin))
        print("Observations - total actual obs: tmax  {}".format(total_actual_obs_tmax))
        print("Observations - Tmin used values: {} {}% of total actual obs".format(tmin_values, pcnt_tmin))
        print("Observations - Tmax used values: {} {}% of total actual obs".format(tmax_values, pcnt_tmax))
        print(" --- ")


        # ---  repeat for single day, no qc flag, with holdout file --- 

        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE, hold_out_filename=FILENAME_HOLD_OUT, daynumber=57344)

        # access the object - not invoking qc flag use
        tmin_hold = result.observations(ObservationSource.TMIN, False)
        tmax_hold = result.observations(ObservationSource.TMAX, False)

        total_obs = result.number_of_observations()

        # calculate total number of actual observations by ignoring fill values
        total_actual_obs_tmin_hold = (numpy.count_nonzero(numpy.where(tmin_hold.measurement > -32768)))
        total_actual_obs_tmax_hold = (numpy.count_nonzero(numpy.where(tmax_hold.measurement > -32768)))

        # total number of valid values - non masked
        tmin_hold_values = numpy.sum(tmin_hold.mask == 0)
        tmax_hold_values = numpy.sum(tmax_hold.mask == 0)

        # calculate valid values as % of the total number of actual observations 
        pcnt_tmin_hold = (tmin_hold_values / float(total_actual_obs_tmin_hold))*100
        pcnt_tmax_hold = (tmax_hold_values / float(total_actual_obs_tmax_hold))*100

        print(" --- ")
        print(" --- Observations - single day - no qc flag use - WITH holdout file --- ")
        print(" --- ")
        print("Observations - days x stations:  {} {}".format(days, stations))
        print("Observations - measurement shape: {}".format(tmin.measurement.shape))
        print("Observations - total potential obs:  {}".format(total_obs))
        print("Observations - total actual obs: tmin  {}".format(total_actual_obs_tmin_hold))
        print("Observations - total actual obs: tmax  {}".format(total_actual_obs_tmax_hold))
        print("Observations - Tmin used values: {} {}% of total actual obs".format(tmin_hold_values, pcnt_tmin_hold))
        print("Observations - Tmax used values: {} {}% of total actual obs".format(tmax_hold_values, pcnt_tmax_hold))
        print(" --- ")


        # ---  repeat for single day, WITH qc flag, WITH holdout file --- 

        # attempt load
        result = ObservationSourceInsituLand(FILENAME_EXAMPLE, hold_out_filename=FILENAME_HOLD_OUT, daynumber=57344)

        # access the object - invoking qc flag use - so more values are masked out 
        tmin_qc_hold = result.observations(ObservationSource.TMIN, True)
        tmax_qc_hold = result.observations(ObservationSource.TMAX, True)

        total_obs = result.number_of_observations()

        # calculate total number of actual observations by ignoring fill values
        total_actual_obs_tmin_qc_hold = (numpy.count_nonzero(numpy.where(tmin_qc_hold.measurement > -32768)))
        total_actual_obs_tmax_qc_hold = (numpy.count_nonzero(numpy.where(tmax_qc_hold.measurement > -32768)))

        # total number of valid values - non masked
        tmin_values_qc_hold = numpy.sum(tmin_qc_hold.mask == 0)
        tmax_values_qc_hold = numpy.sum(tmax_qc_hold.mask == 0)

        # calculate valid values as % of the total number of actual observations 
        pcnt_tmin_qc_hold = (tmin_values_qc_hold / float(total_actual_obs_tmin_qc_hold))*100
        pcnt_tmax_qc_hold = (tmax_values_qc_hold / float(total_actual_obs_tmax_qc_hold))*100


        print(" --- ")
        print(" --- Observations - single day - with qc flags - with holdout file --- ")
        print(" --- ")
        print("Observations - days x stations:  {} {}".format(days, stations))
        print("Observations - measurement shape: {}".format(tmin_qc_hold.measurement.shape))
        print("Observations - total potential obs:  {}".format(total_obs))
        print("Observations - total actual obs: tmin  {}".format(total_actual_obs_tmin_qc_hold))
        print("Observations - total actual obs: tmax  {}".format(total_actual_obs_tmax_qc_hold))
        print("Observations - Tmin used values qc: {} {}% of total actual obs".format(tmin_values_qc_hold, pcnt_tmin_qc_hold))
        print("Observations - Tmax used values qc: {} {}% of total actual obs".format(tmax_values_qc_hold, pcnt_tmax_qc_hold))
        print(" --- ")


        # test there are fewer valid observations once qc flags have been applied
        self.assertLess(tmin_values_qc_hold, total_actual_obs_tmin)
        self.assertLess(tmax_values_qc_hold, total_actual_obs_tmax)

        # test the number of valid, available observations after holdout and qc flag masks
        self.assertEqual(11775, numpy.sum(tmin_qc_hold.mask ==0))
        self.assertEqual(11776, numpy.sum(tmax_qc_hold.mask ==0))

        # test there are more masked values (True) after applying the qc flag masks
        self.assertGreater(numpy.sum(tmin_qc_hold.mask), numpy.sum(tmin.mask))
        self.assertGreater(numpy.sum(tmax_qc_hold.mask), numpy.sum(tmax.mask))

    def test_adjustment(self):
        """Test a location based altitude adjustment using a mockup adjustment object"""
        
        expected_adjustment = 2.0
        expected_adjustment_uncertainty = 3.0
        altitude_adjustment = MockAdjustment(expected_adjustment, expected_adjustment_uncertainty)
        
        # example filename to use
        FILENAME_EXAMPLE = os.path.join(
            eustaceconfig.WORKSPACE_PATH,
            'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')
        
        
        # test adjusted multi day load
        unadjusted_result = ObservationSourceInsituLand(FILENAME_EXAMPLE)
        adjusted_result = ObservationSourceInsituLand(FILENAME_EXAMPLE, altitude_adjustment = altitude_adjustment)
        
        unadjusted_observations = unadjusted_result.observations(ObservationSource.TMAX, True).measurement 
        adjusted_observations = adjusted_result.observations(ObservationSource.TMAX, True).measurement 
        
        mask =  unadjusted_result.observations(ObservationSource.TMAX, True).mask         
        unadjusted_unc = unadjusted_result.observations(ObservationSource.TMAX, True).uncorrelatederror 
        adjusted_unc = adjusted_result.observations(ObservationSource.TMAX, True).uncorrelatederror 
        
        self.assertTrue( numpy.array_equal(adjusted_observations[~mask], unadjusted_observations[~mask]+expected_adjustment) )
        self.assertTrue( numpy.array_equal(adjusted_unc[~mask], numpy.sqrt(unadjusted_unc[~mask]**2 + expected_adjustment_uncertainty**2)) )
        
        # test adjusted single day load
        unadjusted_result = ObservationSourceInsituLand(FILENAME_EXAMPLE, daynumber=57344)
        adjusted_result = ObservationSourceInsituLand(FILENAME_EXAMPLE, daynumber=57344, altitude_adjustment = altitude_adjustment)
        
        unadjusted_observations = unadjusted_result.observations(ObservationSource.TMAX, True).measurement 
        adjusted_observations = adjusted_result.observations(ObservationSource.TMAX, True).measurement 
        
        mask =  unadjusted_result.observations(ObservationSource.TMAX, True).mask         
        unadjusted_unc = unadjusted_result.observations(ObservationSource.TMAX, True).uncorrelatederror 
        adjusted_unc = adjusted_result.observations(ObservationSource.TMAX, True).uncorrelatederror 
        
        self.assertTrue( numpy.array_equal(adjusted_observations[~mask], unadjusted_observations[~mask]+expected_adjustment) )
        self.assertTrue( numpy.array_equal(adjusted_unc[~mask], numpy.sqrt(unadjusted_unc[~mask]**2 + expected_adjustment_uncertainty**2)) )
        
        
        
def example_adjustment():
    """Example use of altitude adjustement for insitu land"""
    
    import matplotlib.pyplot as plt
    from ..height_adjustment import AltitudeAdjustment
    
    altitude_filename = os.path.join(
        eustaceconfig.WORKSPACE_PATH,
        'data/internal/climatology_covariates/DEM_global_0.25_0.25.nc')
    
    # setup using the analysis' DEM covariate file
    filename = os.path.join(eustaceconfig.WORKSPACE_PATH,
                            'data/internal/climatology_covariates/DEM_global_0.25_0.25.nc')
    latitude_label = "lat"
    longitude_label = "lon" 
    covariate_label = "dem"
    rescale_factor = 0.001 # rescale altitude map to to K km^-1

    # initialise the model
    altitude_adjustment = AltitudeAdjustment(filename, latitude_label, longitude_label, covariate_label, rescale_factor)
    
    # example filename to use
    FILENAME_EXAMPLE = os.path.join(
        eustaceconfig.WORKSPACE_PATH,
        'data/internal/D1.7/daily/eustace_stations_global_R001124_2007_daily_temperature.nc')

    # attempt load
    unadjusted_result = ObservationSourceInsituLand(FILENAME_EXAMPLE, daynumber=57344)
    adjusted_result = ObservationSourceInsituLand(FILENAME_EXAMPLE, daynumber=57344, altitude_adjustment = altitude_adjustment)
    
    unadjusted_observations = unadjusted_result.observations(ObservationSource.TMAX, True).measurement 
    adjusted_observations = adjusted_result.observations(ObservationSource.TMAX, True).measurement 
    
    mask =  unadjusted_result.observations(ObservationSource.TMAX, True).mask 
    
    unadjusted_unc = unadjusted_result.observations(ObservationSource.TMAX, True).uncorrelatederror 
    adjusted_unc = adjusted_result.observations(ObservationSource.TMAX, True).uncorrelatederror 
    
    locations = unadjusted_result.observation_location_lookup()
    mask = ~adjusted_result.observations(ObservationSource.TMAX, True).mask[:,0]
    
    # plot output
    plt.figure()
    plt.scatter(locations[1,mask], locations[0,mask], c = adjusted_observations[mask,0] - unadjusted_observations[mask,0], s=1.5, vmin = -3.25, vmax = 3.25)
    plt.colorbar()
    
    plt.figure()
    plt.scatter(locations[1,mask], locations[0,mask], c =  numpy.sqrt(adjusted_unc[mask,0]**2 - unadjusted_unc[mask,0]**2), s=1.5, vmin = 0.0, vmax = 0.55)
    plt.colorbar()

    plt.show()
    
if __name__ == '__main__':
    
    example_adjustment()