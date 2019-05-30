"""Tests for ice surface model file I/O."""

import unittest
import numpy
import os
import eustaceconfig
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from ..ice import IceSurfaceTemperatureQualityControlNetCDF
from ..ice import ObservationSourceIceDMI
from ..ice import ObservationSourceIceDMICombined

class TestIceSurfaceTemperatureQualityControlNetCDF(unittest.TestCase):

    def test_init(self):

        # Test data
        fileA = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/AASTI_v1_L3_solartime/2007/01/2007-01-01_sh.nc')
        fileB = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/AASTI_v1_L3_solartime/2003/08/2003-08-03_nh.nc')

        # Build
        resultA = IceSurfaceTemperatureQualityControlNetCDF(fileA)
        resultB = IceSurfaceTemperatureQualityControlNetCDF(fileB)

        # Check the auto-detected fields
        self.assertEqual(57343, resultA.daynumber)
        self.assertEqual('sh', resultA.hemisphere())
        self.assertEqual(56096, resultB.daynumber)
        self.assertEqual('nh', resultB.hemisphere())

        # Quick check of some data (retrieved using ncdump)
        numpy.testing.assert_equal(resultA.surface_temperature.mask[0,0:5], [ True, False, False, True, False ])
        self.assertAlmostEqual(resultA.surface_temperature.data[0,1], 250.2197265625)
        self.assertAlmostEqual(resultA.surface_temperature.data[0,2], 247.509765625)
        self.assertAlmostEqual(resultA.surface_temperature.data[0,4], 244.3701171875)


class TestObservationSourceIceDMI(unittest.TestCase):

    def test_load_example(self):

        # example filename to use
        # FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/MS8_version2/Ice/data_v2.0/2007/eustace_satellite_test2ice_nh_20070201.nc'
        # FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_nh_20070201.nc'
        FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_sh_20070201.nc'

        # attempt load
        result = ObservationSourceIceDMI(FILENAME_EXAMPLE)

        # supports all default observations
        self.assertEqual(['Tmean', 'Tmin', 'Tmax'], result.observables())

        # get number of observations
        self.assertEqual(230400, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 2 rows of 230400 columns
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, 230400), location_lookup.shape)

        # check ranges
        ranges = result.local_correlation_length_scale(ObservationSource.TMEAN)
        self.assertTrue(isinstance(ranges, list))

        # check shape and type of observations
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertTrue(isinstance(tmean, Observations))
        self.assertTrue(isinstance(tmean.mask, numpy.ndarray))
        self.assertTrue(tmean.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmean.measurement, numpy.ndarray))
        self.assertTrue(tmean.measurement.dtype == numpy.float32)
        self.assertEqual((230400,), tmean.measurement.shape)
        self.assertTrue(isinstance(tmean.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmean.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((230400,), tmean.uncorrelatederror.shape)
        self.assertEqual(1, len(tmean.locallycorrelatederror))
        self.assertTrue(isinstance(tmean.locallycorrelatederror[0], numpy.ndarray))
        self.assertTrue(tmean.locallycorrelatederror[0].dtype == numpy.float32)
        self.assertEqual((230400,), tmean.locallycorrelatederror[0].shape)

        tmin = result.observations(ObservationSource.TMIN)
        tmax = result.observations(ObservationSource.TMAX)

        # check values
        # (from ncks)
        #
        # -- Observations
        #
        # time[0]=57374 latitude[1]=-89.625 longitude[82]=-159.375 tas[1522]=246.337 K
        # time[0]=57374 latitude[104]=-63.875 longitude[484]=-58.875 tas[150244]=269.353 K
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[26]=-173.375 tasmin[26]=248.647 K
        # time[0]=57374 latitude[115]=-61.125 longitude[500]=-54.875 tasmin[166100]=268.285 K
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[1]=-179.625 tasmax[1]=245.438 K
        # time[0]=57374 latitude[110]=-62.375 longitude[478]=-60.375 tasmax[158878]=270.439 K
        #
        # -- Uncorrelated uncertainty
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[1]=-179.625 RU[1]=1.67472 K
        # time[0]=57374 latitude[1]=-89.625 longitude[82]=-159.375 RU[1522]=1.67362 K
        # time[0]=57374 latitude[115]=-61.125 longitude[500]=-54.875 RU[166100]=1.67154 K
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[1]=-179.625 RUmin[1]=1.52424 K
        # time[0]=57374 latitude[0]=-89.875 longitude[26]=-173.375 RUmin[26]=1.52333 K
        # time[0]=57374 latitude[115]=-61.125 longitude[500]=-54.875 RUmin[166100]=1.52095 K
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[1]=-179.625 RUmax[1]=3.8489 K
        # time[0]=57374 latitude[110]=-62.375 longitude[478]=-60.375 RUmax[158878]=3.84652 K
        # time[0]=57374 latitude[115]=-61.125 longitude[500]=-54.875 RUmax[166100]=3.84759 K
        #
        # -- Locally correlated uncertainty
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[1]=-179.625 SSU[1]=1.52248
        # time[0]=57374 latitude[1]=-89.625 longitude[82]=-159.375 SSU[1522]=1.51589 K
        # time[0]=57374 latitude[104]=-63.875 longitude[484]=-58.875 SSU[150244]=1.51363 K
        # time[0]=57374 latitude[115]=-61.125 longitude[500]=-54.875 SSU[166100]=1.51166 K
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[1]=-179.625 SSUmin[1]=1.81769 K
        # time[0]=57374 latitude[0]=-89.875 longitude[26]=-173.375 SSUmin[26]=1.81634 K
        # time[0]=57374 latitude[115]=-61.125 longitude[500]=-54.875 SSUmin[166100]=1.80916 K
        #
        # time[0]=57374 latitude[0]=-89.875 longitude[1]=-179.625 SSUmax[1]=2.01594 K
        # time[0]=57374 latitude[110]=-62.375 longitude[478]=-60.375 SSUmax[158878]=2.00499 K
        # time[0]=57374 latitude[115]=-61.125 longitude[500]=-54.875 SSUmax[166100]=2.00825 K

        # Check masks bracket the valid period
        self.assertEqual( True, tmean.mask[1521])
        self.assertEqual(False, tmean.mask[1522])
        self.assertEqual(False, tmean.mask[150244])
        self.assertEqual( True, tmean.mask[150245])
        self.assertEqual( True, tmin.mask[25])
        self.assertEqual(False, tmin.mask[26])
        self.assertEqual(False, tmin.mask[166100])
        self.assertEqual( True, tmin.mask[166101])
        self.assertEqual( True, tmax.mask[0])
        self.assertEqual(False, tmax.mask[1])
        self.assertEqual(False, tmax.mask[158878])
        self.assertEqual( True, tmax.mask[158879])

        # Values
        self.assertAlmostEqual(246.337, tmean.measurement[1522], places=3)
        self.assertAlmostEqual(269.353, tmean.measurement[150244], places=3)
        self.assertAlmostEqual(248.647, tmin.measurement[26], places=3)
        self.assertAlmostEqual(268.285, tmin.measurement[166100], places=3)
        self.assertAlmostEqual(245.438, tmax.measurement[1], places=3)
        self.assertAlmostEqual(270.439, tmax.measurement[158878], places=3)
        self.assertAlmostEqual(1.67362, tmean.uncorrelatederror[1522], places=5)
        self.assertAlmostEqual(1.66910, tmean.uncorrelatederror[150244], places=5)
        self.assertAlmostEqual(1.52333, tmin.uncorrelatederror[26], places=5)
        self.assertAlmostEqual(1.52095, tmin.uncorrelatederror[166100], places=5)
        self.assertAlmostEqual(3.84890, tmax.uncorrelatederror[1], places=5)
        self.assertAlmostEqual(3.84652, tmax.uncorrelatederror[158878], places=5)
        self.assertAlmostEqual(1.51589, tmean.locallycorrelatederror[0][1522], places=5)
        self.assertAlmostEqual(1.51363, tmean.locallycorrelatederror[0][150244], places=5)
        self.assertAlmostEqual(1.81634, tmin.locallycorrelatederror[0][26], places=5)
        self.assertAlmostEqual(1.80916, tmin.locallycorrelatederror[0][166100], places=5)
        self.assertAlmostEqual(2.01594, tmax.locallycorrelatederror[0][1], places=5)
        self.assertAlmostEqual(2.00499, tmax.locallycorrelatederror[0][158878], places=5)

class TestObservationSourceIceDMICombined(unittest.TestCase):

    def test_files_wrong_way_around(self):

        # example filenames to use
        FILENAME_EXAMPLE_NH = '/gws/nopw/j04/eustace/data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_nh_20070201.nc'
        FILENAME_EXAMPLE_SH = '/gws/nopw/j04/eustace/data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_sh_20070201.nc'

        # attempt load the wrong way around
        with self.assertRaises(ValueError):
            ObservationSourceIceDMICombined([FILENAME_EXAMPLE_NH, FILENAME_EXAMPLE_SH])

    def test_load_example_combined(self):
        
        # example filenames to use
        FILENAME_EXAMPLE_NH = '/gws/nopw/j04/eustace/data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_nh_20070201.nc'
        FILENAME_EXAMPLE_SH = '/gws/nopw/j04/eustace/data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_sh_20070201.nc'

        # attempt load
        result = ObservationSourceIceDMICombined([FILENAME_EXAMPLE_SH, FILENAME_EXAMPLE_NH])

        # supports all default observations
        self.assertEqual(['Tmean', 'Tmin', 'Tmax'], result.observables())

        # get number of observations
        self.assertEqual(460800, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 2 rows at quarter degree resolution
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, (4 * 180) * (4 * 360)), location_lookup.shape)

        # Expected values of locations are first 40 degrees from south pole, then gap, then last 40 degrees to north pole
        numpy.testing.assert_almost_equal(location_lookup[:,0:2], [ [ -89.875, -89.875 ], [ -179.875, -179.625 ] ])

        # check ranges
        ranges = result.local_correlation_length_scale(ObservationSource.TMEAN)
        self.assertTrue(isinstance(ranges, list))

        # check shape and type of observations
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertTrue(isinstance(tmean, Observations))
        self.assertTrue(isinstance(tmean.mask, numpy.ndarray))
        self.assertTrue(tmean.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmean.measurement, numpy.ndarray))
        self.assertTrue(tmean.measurement.dtype == numpy.float32)
        self.assertEqual((460800,), tmean.measurement.shape)
        self.assertTrue(isinstance(tmean.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmean.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((460800,), tmean.uncorrelatederror.shape)
        self.assertEqual(1, len(tmean.locallycorrelatederror))
        self.assertTrue(isinstance(tmean.locallycorrelatederror[0], numpy.ndarray))
        self.assertTrue(tmean.locallycorrelatederror[0].dtype == numpy.float32)
        self.assertEqual((460800,), tmean.locallycorrelatederror[0].shape)

        tmin = result.observations(ObservationSource.TMIN)
        tmax = result.observations(ObservationSource.TMAX)

        # Check masks bracket the valid period for southern hemisphere
        self.assertEqual( True, tmean.mask[1521])
        self.assertEqual(False, tmean.mask[1522])
        self.assertEqual(False, tmean.mask[150244])
        self.assertEqual( True, tmean.mask[150245])
        self.assertEqual( True, tmin.mask[25])
        self.assertEqual(False, tmin.mask[26])
        self.assertEqual(False, tmin.mask[166100])
        self.assertEqual( True, tmin.mask[166101])
        self.assertEqual( True, tmax.mask[0])
        self.assertEqual(False, tmax.mask[1])
        self.assertEqual(False, tmax.mask[158878])
        self.assertEqual( True, tmax.mask[158879])

        # Values from southern hemisphere (as in test_load_example above)
        self.assertAlmostEqual(246.337, tmean.measurement[1522], places=3)
        self.assertAlmostEqual(269.353, tmean.measurement[150244], places=3)
        self.assertAlmostEqual(248.647, tmin.measurement[26], places=3)
        self.assertAlmostEqual(268.285, tmin.measurement[166100], places=3)
        self.assertAlmostEqual(245.438, tmax.measurement[1], places=3)
        self.assertAlmostEqual(270.439, tmax.measurement[158878], places=3)
        self.assertAlmostEqual(1.67362, tmean.uncorrelatederror[1522], places=5)
        self.assertAlmostEqual(1.66910, tmean.uncorrelatederror[150244], places=5)
        self.assertAlmostEqual(1.52333, tmin.uncorrelatederror[26], places=5)
        self.assertAlmostEqual(1.52095, tmin.uncorrelatederror[166100], places=5)
        self.assertAlmostEqual(3.84890, tmax.uncorrelatederror[1], places=5)
        self.assertAlmostEqual(3.84652, tmax.uncorrelatederror[158878], places=5)
        self.assertAlmostEqual(1.51589, tmean.locallycorrelatederror[0][1522], places=5)
        self.assertAlmostEqual(1.51363, tmean.locallycorrelatederror[0][150244], places=5)
        self.assertAlmostEqual(1.81634, tmin.locallycorrelatederror[0][26], places=5)
        self.assertAlmostEqual(1.80916, tmin.locallycorrelatederror[0][166100], places=5)
        self.assertAlmostEqual(2.01594, tmax.locallycorrelatederror[0][1], places=5)
        self.assertAlmostEqual(2.00499, tmax.locallycorrelatederror[0][158878], places=5)

        # Check that values from northern hemisphere correspond to results of loading that file on its own
        # It has the last 40 degrees of latitudes
        # - corresponding to grid rows 560...719 at quarter degree resolution
        # - which should be indices 806400...1036799
        northern = ObservationSourceIceDMI(FILENAME_EXAMPLE_NH)
        northern_tmean = northern.observations(ObservationSource.TMEAN)
        numpy.testing.assert_almost_equal(result.observation_location_lookup()[:, 806400:], northern.observation_location_lookup())
        numpy.testing.assert_almost_equal(tmean.measurement[230400:], northern_tmean.measurement)
        numpy.testing.assert_almost_equal(tmean.uncorrelatederror[230400:], northern_tmean.uncorrelatederror)
        numpy.testing.assert_almost_equal(tmean.locallycorrelatederror[0][230400:], northern_tmean.locallycorrelatederror[0])
        numpy.testing.assert_almost_equal(result.observation_location_lookup()[:, tmean.location[230400:] ], northern.observation_location_lookup()[:, northern_tmean.location])
