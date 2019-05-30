"""Tests for ocean surface model file I/O."""

import unittest
import numpy
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from ..ocean import ObservationSourceOceanMetOfficeWP1

class TestObservationSourceOceanMetOfficeWP1(unittest.TestCase):

    def test_load_example(self):

        # example filename to use
        # FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/MS8_version2/Ocean/data/2007/eustace_satellite_ocean_20070207.nc'
        FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/mat_from_sst/2007/200701071200-MAT_from_ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc'

        # attempt load
        result = ObservationSourceOceanMetOfficeWP1(FILENAME_EXAMPLE)

        # only supports tmean
        self.assertEqual([ ObservationSource.TMEAN ], result.observables())

        # get number of observations
        self.assertEqual(1036800, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # should be 3 rows and one column per grid box
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        self.assertEqual((2, 1036800), location_lookup.shape)

        # check ranges
        ranges = result.local_correlation_length_scale(ObservationSource.TMEAN)
        self.assertTrue(isinstance(ranges, list))

        # check shape and type of observations
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertTrue(isinstance(tmean, Observations))
        self.assertTrue(isinstance(tmean.mask, numpy.ndarray))
        self.assertTrue(tmean.mask.dtype == numpy.bool)
        self.assertEqual(numpy.float32, tmean.time.dtype)
        self.assertTrue(isinstance(tmean.location, numpy.ndarray))
        self.assertTrue(tmean.location.dtype == numpy.uint64)
        self.assertEqual((1036800,), tmean.location.shape)
        self.assertTrue(isinstance(tmean.measurement, numpy.ndarray))
        self.assertTrue(tmean.measurement.dtype == numpy.float64)
        self.assertEqual((1036800,), tmean.measurement.shape)
        self.assertTrue(isinstance(tmean.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmean.uncorrelatederror.dtype == numpy.float64)
        self.assertEqual((1036800,), tmean.uncorrelatederror.shape)
        self.assertEqual(1, len(tmean.locallycorrelatederror))
        self.assertTrue(isinstance(tmean.locallycorrelatederror[0], numpy.ndarray))
        self.assertTrue(tmean.locallycorrelatederror[0].dtype == numpy.float64)
        self.assertEqual((1036800,), tmean.locallycorrelatederror[0].shape)

        # Check values
        #
        # Used ncks to look at NetCDF and establish output values
        # ncks -v tas /gws/nopw/j04/eustace/data/internal/mat_from_sst/2007/200701071200-MAT_from_ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc | grep -v _
        # This gives:
        #
        # y[103] x[471] tas[148791]=17 K
        # y[103] x[472] tas[148792]=14 K
        # y[104] x[471] tas[150231]=22 K
        # y[104] x[472] tas[150232]=14 K
        # y[104] x[473] tas[150233]=12 K
        # y[105] x[473] tas[151673]=6 K
        # y[106] x[474] tas[153114]=10 K
        # y[106] x[475] tas[153115]=8 K
        #
        # ...
        #
        # y[663] x[805] tas[955525]=-52 K
        # y[663] x[828] tas[955548]=-53 K
        # y[663] x[832] tas[955552]=-37 K
        # y[663] x[833] tas[955553]=-46 K
        # y[663] x[834] tas[955554]=-53 K
        # y[664] x[827] tas[956987]=-59 K
        # y[664] x[833] tas[956993]=-75 K
        # y[664] x[834] tas[956994]=-78 K
        #
        # Note that the K is erroneous in the sense that these integers are unscaled
        # - must multiply by scale_factor 0.1 to get values
        # - and then the values are in deg C (not K)
        #
        # Total number of rows produced is: 57258 (looked at results with emacs to get this)
        #

        # Check total valid observations count
        self.assertEqual(57258, numpy.nonzero(tmean.mask == False)[0].shape[0])
        
        # Check masks of first and last three
        self.assertEqual(False, tmean.mask[148791])
        self.assertEqual(False, tmean.mask[148792])
        self.assertEqual(False, tmean.mask[150231])
        self.assertEqual(False, tmean.mask[956987])
        self.assertEqual(False, tmean.mask[956993])
        self.assertEqual(False, tmean.mask[956994])

        # Check data of first and last three
        self.assertAlmostEqual(274.85, tmean.measurement[148791], places=4)
        self.assertAlmostEqual(274.55, tmean.measurement[148792], places=4)
        self.assertAlmostEqual(275.35, tmean.measurement[150231], places=4)
        self.assertAlmostEqual(267.25, tmean.measurement[956987], places=4)
        self.assertAlmostEqual(265.65, tmean.measurement[956993], places=4)
        self.assertAlmostEqual(265.35, tmean.measurement[956994], places=4)

        # Check uncorrelated uncertainty
        # (Scale factor 0.001, offset 0.0)
        # y[103] x[471] unc_rand_tas[148791]=286
        # ...
        # y[664] x[834] unc_rand_tas[956994]=221
        self.assertAlmostEqual(0.286, tmean.uncorrelatederror[148791], places=4)
        self.assertAlmostEqual(0.221, tmean.uncorrelatederror[956994], places=4)

        # Check  locally correlated uncertainty
        # (Scale factor 0.001, offset 0.0)
        # y[103] x[471] unc_corr_tas[148791]=46
        # ...
        # y[664] x[834] unc_corr_tas[956994]=42
        self.assertAlmostEqual(0.046, tmean.locallycorrelatederror[0][148791], places=4)
        self.assertAlmostEqual(0.042, tmean.locallycorrelatederror[0][956994], places=4)
