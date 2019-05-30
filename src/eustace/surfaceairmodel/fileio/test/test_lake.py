"""Tests for lake surface model file I/O."""

import unittest
import numpy
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from ..lake import ObservationSourceLakeReading

class TestObservationSourceLakeReading(unittest.TestCase):

    def test_load_example(self):

        # example filename to use
        FILENAME_EXAMPLE = '/gws/nopw/j04/eustace/data/internal/MS8_version2/Lake/data/fTa.txt'

        # attempt load
        result = ObservationSourceLakeReading(FILENAME_EXAMPLE)

        # only supports tmean
        self.assertEqual([ ObservationSource.TMEAN ], result.observables())

        # get number of observations
        self.assertEqual(7300, result.number_of_observations())
        
        # get coordinates
        location_lookup = result.observation_location_lookup()

        # in this file there is just 1 unique location
        self.assertTrue(isinstance(location_lookup, numpy.ndarray))
        numpy.testing.assert_almost_equal([ [ 47.585 ], [-86.585 ] ], location_lookup, decimal=4)

        # check shape and type of observations
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertTrue(isinstance(tmean, Observations))
        self.assertTrue(isinstance(tmean.mask, numpy.ndarray))
        self.assertEqual(tmean.mask.dtype, numpy.bool)
        self.assertTrue(isinstance(tmean.time, numpy.ndarray))
        self.assertEqual(tmean.time.dtype, numpy.float32)
        self.assertEqual((7300,), tmean.time.shape)
        self.assertTrue(isinstance(tmean.location, numpy.ndarray))
        self.assertEqual(tmean.location.dtype, numpy.uint64)
        numpy.testing.assert_equal(numpy.zeros((7300,), numpy.uint64), tmean.location)
        self.assertTrue(isinstance(tmean.measurement, numpy.ndarray))
        self.assertEqual(tmean.measurement.dtype, numpy.float32)
        self.assertEqual((7300,), tmean.measurement.shape)
        self.assertTrue(isinstance(tmean.uncorrelatederror, numpy.ndarray))
        self.assertEqual(tmean.uncorrelatederror.dtype, numpy.float32)
        self.assertEqual((7300,), tmean.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmean.locallycorrelatederror, list))
        self.assertEqual(0, len(tmean.locallycorrelatederror))

        # check empty array of correlation length scales
        length_scale = result.local_correlation_length_scale(ObservationSource.TMEAN)
        self.assertTrue(isinstance(length_scale, numpy.ndarray))
        self.assertEqual(0, len(length_scale))

        # mask should be all in view
        self.assertTrue(numpy.all(numpy.logical_not(tmean.mask)))

        # check first two values
        self.assertAlmostEqual(267.9951, tmean.measurement[0], places=4)
        self.assertAlmostEqual(268.0673, tmean.measurement[1], places=4)

        # check last value
        self.assertAlmostEqual(271.9577, tmean.measurement[7299], places=4)

        # at present the uncorrelated error is set to 5 everywhere
        self.assertTrue(numpy.all(tmean.uncorrelatederror == numpy.float32(5.0)))
