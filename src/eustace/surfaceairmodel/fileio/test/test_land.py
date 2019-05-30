"""Tests for land surface model file I/O."""

import unittest
import numpy
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from ..land import ObservationSourceLandMetOfficeWP1

class TestObservationSourceLandMetOfficeWP1(unittest.TestCase):

    def test_load_example(self):

        # example filename to use
        # FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/MS8_version2/Land/data/eustace_satellite_3.000_20070207.nc'
        FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/MS9/Land/data/satellite_lsat/eustace_satellite_4.000_20070207.nc'

        # attempt load
        result = ObservationSourceLandMetOfficeWP1(FILENAME_EXAMPLE)

        # supports all default observations
        self.assertEqual(['Tmean', 'Tmin', 'Tmax'], result.observables())

        # get number of observations
        self.assertEqual(1036800, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 3 rows and one column per grid box
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual(numpy.float64, location_lookup.dtype)
        self.assertEqual((2,1036800), location_lookup.shape)

        # check ranges
        ranges = result.local_correlation_length_scale(ObservationSource.TMEAN)
        self.assertTrue(isinstance(ranges, list))

        # check shape and type of observations (Tmin)
        tmin = result.observations(ObservationSource.TMIN)
        self.assertTrue(isinstance(tmin, Observations))
        self.assertTrue(isinstance(tmin.mask, numpy.ndarray))
        self.assertTrue(tmin.mask.dtype == numpy.bool)
        self.assertEqual(numpy.int32, tmin.time.dtype)
        self.assertTrue(isinstance(tmin.location, numpy.ndarray))
        self.assertTrue(tmin.location.dtype == numpy.uint64)
        self.assertEqual((1036800,), tmin.location.shape)
        self.assertTrue(isinstance(tmin.measurement, numpy.ndarray))
        self.assertTrue(tmin.measurement.dtype == numpy.float32)
        self.assertEqual((1036800,), tmin.measurement.shape)
        self.assertTrue(isinstance(tmin.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmin.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((1036800,), tmin.uncorrelatederror.shape)
        self.assertEqual(1, len(tmin.locallycorrelatederror))
        self.assertTrue(isinstance(tmin.locallycorrelatederror[0], numpy.ndarray))
        self.assertTrue(tmin.locallycorrelatederror[0].dtype == numpy.float32)
        self.assertEqual((1036800,), tmin.locallycorrelatederror[0].shape)

        # check shape and type of observations (Tmax)
        tmax = result.observations(ObservationSource.TMAX)
        self.assertTrue(isinstance(tmax, Observations))
        self.assertTrue(isinstance(tmax.mask, numpy.ndarray))
        self.assertTrue(tmax.mask.dtype == numpy.bool)
        self.assertEqual(numpy.int32, tmax.time.dtype)
        self.assertTrue(isinstance(tmax.location, numpy.ndarray))
        self.assertTrue(tmax.location.dtype == numpy.uint64)
        self.assertEqual((1036800,), tmax.location.shape)
        self.assertTrue(isinstance(tmax.measurement, numpy.ndarray))
        self.assertTrue(tmax.measurement.dtype == numpy.float32)
        self.assertEqual((1036800,), tmax.measurement.shape)
        self.assertTrue(isinstance(tmax.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmax.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((1036800,), tmax.uncorrelatederror.shape)
        self.assertEqual(1, len(tmax.locallycorrelatederror))
        self.assertTrue(isinstance(tmax.locallycorrelatederror[0], numpy.ndarray))
        self.assertTrue(tmax.locallycorrelatederror[0].dtype == numpy.float32)
        self.assertEqual((1036800,), tmax.locallycorrelatederror[0].shape)

        # check shape and type of observations (Tmean)
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertTrue(isinstance(tmean, Observations))
        self.assertTrue(isinstance(tmean.mask, numpy.ndarray))
        self.assertTrue(tmean.mask.dtype == numpy.bool)
        self.assertEqual(numpy.int32, tmean.time.dtype)
        self.assertTrue(isinstance(tmean.location, numpy.ndarray))
        self.assertTrue(tmean.location.dtype == numpy.uint64)
        self.assertEqual((1036800,), tmean.location.shape)
        self.assertTrue(isinstance(tmean.measurement, numpy.ndarray))
        self.assertTrue(tmean.measurement.dtype == numpy.float32)
        self.assertEqual((1036800,), tmean.measurement.shape)
        self.assertTrue(isinstance(tmean.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmean.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((1036800,), tmean.uncorrelatederror.shape)
        self.assertEqual(1, len(tmean.locallycorrelatederror))
        self.assertTrue(isinstance(tmean.locallycorrelatederror[0], numpy.ndarray))
        self.assertTrue(tmean.locallycorrelatederror[0].dtype == numpy.float32)
        self.assertEqual((1036800,), tmean.locallycorrelatederror[0].shape)

        # Check values for tmin
        #
        # Used ncks to look at NetCDF and establish output values
        # ncks -v tasmin /gws/nopw/j04/eustace/data/internal/MS9/Land/data/satellite_lsat/eustace_satellite_4.000_20070207.nc  | grep -v _
        #
        # For tasmin:
        #
        # lat[40]=-79.875 lon[775]=13.875 time[0]=57380 tasmin[58375]=-7648 K
        # lat[40]=-79.875 lon[776]=14.125 time[0]=57380 tasmin[58376]=-7571 K
        # lat[40]=-79.875 lon[777]=14.375 time[0]=57380 tasmin[58377]=-7435 K
        # lat[40]=-79.875 lon[778]=14.625 time[0]=57380 tasmin[58378]=-7535 K
        # ...
        # lat[639]=69.875 lon[143]=-144.125 time[0]=57380 tasmin[920303]=-3641 K
        # lat[639]=69.875 lon[219]=-125.125 time[0]=57380 tasmin[920379]=-5632 K
        # lat[639]=69.875 lon[220]=-124.875 time[0]=57380 tasmin[920380]=-5672 K
        # lat[639]=69.875 lon[221]=-124.625 time[0]=57380 tasmin[920381]=-5727 K
        #
        # And total rows = 149044
        #
        # Note units are erroneous (must scale and offset before we get K)
        # Scale: 0.005
        # Offset: 273.15
        #

        # Check total valid obs
        self.assertEqual(149044, numpy.nonzero(tmin.mask == False)[0].shape[0])
        
        # Masks for invalid vs. valid ones
        self.assertEqual( True, tmin.mask[58374])
        self.assertEqual(False, tmin.mask[58375])
        self.assertEqual(False, tmin.mask[58376])
        self.assertEqual(False, tmin.mask[58377])
        self.assertEqual(False, tmin.mask[58378])
        self.assertEqual(False, tmin.mask[920303])
        self.assertEqual( True, tmin.mask[920304])
        self.assertEqual( True, tmin.mask[920378])
        self.assertEqual(False, tmin.mask[920379])
        self.assertEqual(False, tmin.mask[920380])
        self.assertEqual(False, tmin.mask[920381])

        # Data values
        self.assertEqual(-7648, numpy.round((tmin.measurement[ 58375] - 273.15) / 0.005))
        self.assertEqual(-7571, numpy.round((tmin.measurement[ 58376] - 273.15) / 0.005))
        self.assertEqual(-7435, numpy.round((tmin.measurement[ 58377] - 273.15) / 0.005))
        self.assertEqual(-7535, numpy.round((tmin.measurement[ 58378] - 273.15) / 0.005))
        self.assertEqual(-3641, numpy.round((tmin.measurement[920303] - 273.15) / 0.005))
        self.assertEqual(-5632, numpy.round((tmin.measurement[920379] - 273.15) / 0.005))
        self.assertEqual(-5672, numpy.round((tmin.measurement[920380] - 273.15) / 0.005))
        self.assertEqual(-5727, numpy.round((tmin.measurement[920381] - 273.15) / 0.005))


        # Check values for tmax
        #
        # Used ncks to look at NetCDF and establish output values
        # ncks -v tasmax /gws/nopw/j04/eustace/data/internal/MS9/Land/data/satellite_lsat/eustace_satellite_4.000_20070207.nc  | grep -v _
        #
        # For tasmax:
        #
        # lat[40]=-79.875 lon[775]=13.875 time[0]=57380 tasmax[58375]=-4605 K
        # lat[40]=-79.875 lon[776]=14.125 time[0]=57380 tasmax[58376]=-4536 K
        # lat[40]=-79.875 lon[777]=14.375 time[0]=57380 tasmax[58377]=-4413 K
        # lat[40]=-79.875 lon[778]=14.625 time[0]=57380 tasmax[58378]=-4501 K
        # ...
        # lat[639]=69.875 lon[143]=-144.125 time[0]=57380 tasmax[920303]=-896 K
        # lat[639]=69.875 lon[219]=-125.125 time[0]=57380 tasmax[920379]=-2714 K
        # lat[639]=69.875 lon[220]=-124.875 time[0]=57380 tasmax[920380]=-2829 K
        # lat[639]=69.875 lon[221]=-124.625 time[0]=57380 tasmax[920381]=-2924 K
        #
        # And total rows = 149044 (same as tasmin)
        #
        # Note units are erroneous (must scale and offset before we get K)
        # Scale: 0.005
        # Offset: 273.15
        #

        # Check total valid obs
        self.assertEqual(149044, numpy.nonzero(tmax.mask == False)[0].shape[0])
        
        # Masks for invalid vs. valid ones
        self.assertEqual( True, tmax.mask[58374])
        self.assertEqual(False, tmax.mask[58375])
        self.assertEqual(False, tmax.mask[58376])
        self.assertEqual(False, tmax.mask[58377])
        self.assertEqual(False, tmax.mask[58378])
        self.assertEqual(False, tmax.mask[920303])
        self.assertEqual( True, tmax.mask[920304])
        self.assertEqual( True, tmax.mask[920378])
        self.assertEqual(False, tmax.mask[920379])
        self.assertEqual(False, tmax.mask[920380])
        self.assertEqual(False, tmax.mask[920381])

        # Data values
        self.assertEqual(-4605, numpy.round((tmax.measurement[ 58375] - 273.15) / 0.005))
        self.assertEqual(-4536, numpy.round((tmax.measurement[ 58376] - 273.15) / 0.005))
        self.assertEqual(-4413, numpy.round((tmax.measurement[ 58377] - 273.15) / 0.005))
        self.assertEqual(-4501, numpy.round((tmax.measurement[ 58378] - 273.15) / 0.005))
        self.assertEqual( -896, numpy.round((tmax.measurement[920303] - 273.15) / 0.005))
        self.assertEqual(-2714, numpy.round((tmax.measurement[920379] - 273.15) / 0.005))
        self.assertEqual(-2829, numpy.round((tmax.measurement[920380] - 273.15) / 0.005))
        self.assertEqual(-2924, numpy.round((tmax.measurement[920381] - 273.15) / 0.005))


        # Uncorrelated uncertainty
        #
        # Offset: 32.767
        # Scale : 0.001
        #
        # unc_ran_tasmin
        #
        # lat[40]=-79.875 lon[775]=13.875 time[0]=57380 unc_rand_tasmin[58375]=-32668    (= 0.099)
        # ...
        # lat[639]=69.875 lon[221]=-124.625 time[0]=57380 unc_rand_tasmin[920381]=-32658 (= 0.109)
        #
        #
        # unc_ran_tasmax
        #
        # lat[40]=-79.875 lon[775]=13.875 time[0]=57380 unc_rand_tasmax[58375]=-32679    (= 0.088)
        # ...
        # lat[639]=69.875 lon[221]=-124.625 time[0]=57380 unc_rand_tasmax[920381]=-32670 (= 0.097)
        #
        self.assertAlmostEqual(0.099, tmin.uncorrelatederror[ 58375], places=4)
        self.assertAlmostEqual(0.109, tmin.uncorrelatederror[920381], places=4)
        self.assertAlmostEqual(0.088, tmax.uncorrelatederror[ 58375], places=4)
        self.assertAlmostEqual(0.097, tmax.uncorrelatederror[920381], places=4)

        # Locally correlated uncertainty
        #
        # Offset: 32.767
        # Scale : 0.001
        #
        # unc_corr_tasmin
        #
        # lat[40]=-79.875 lon[775]=13.875 time[0]=57380 unc_corr_tasmin[58375]=-29826    (=2.941)
        # ...
        # lat[639]=69.875 lon[221]=-124.625 time[0]=57380 unc_corr_tasmin[920381]=-29852 (=2.915)
        #
        #
        # unc_corr_tasmax
        #
        # lat[40]=-79.875 lon[775]=13.875 time[0]=57380 unc_corr_tasmax[58375]=-28763    (=4.004)
        # ...
        # lat[639]=69.875 lon[221]=-124.625 time[0]=57380 unc_corr_tasmax[920381]=-28778 (=3.989)
        #
        self.assertAlmostEqual(2.941, tmin.locallycorrelatederror[0][ 58375], places=4)
        self.assertAlmostEqual(2.915, tmin.locallycorrelatederror[0][920381], places=4)
        self.assertAlmostEqual(4.004, tmax.locallycorrelatederror[0][ 58375], places=4)
        self.assertAlmostEqual(3.989, tmax.locallycorrelatederror[0][920381], places=4)
