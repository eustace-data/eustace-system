"""Test interface to insitu ocean observations."""

import unittest
import numpy
import os
import eustaceconfig
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from ..insitu_ocean import HadNMAT2Format
from ..insitu_ocean import ObservationSourceInsituOceanHadNMAT2

class TestHadNMAT2Format(unittest.TestCase):

    def test_init(self):

        FORMAT_EXAMPLE = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/Docs/format.txt')
        format =  HadNMAT2Format(FORMAT_EXAMPLE)
        self.assertEqual(8, len(format.usecols))
        self.assertTrue(isinstance( format.usecols[0], int))
        self.assertTrue(isinstance( format.usecols[1], int))
        self.assertTrue(isinstance( format.usecols[2], int))
        self.assertTrue(isinstance( format.usecols[3], int))
        self.assertTrue(isinstance( format.usecols[4], int))
        self.assertTrue(isinstance( format.usecols[5], int))
        self.assertTrue(isinstance( format.usecols[6], int))
        self.assertTrue(isinstance( format.usecols[7], int))

class TestObservationSourceInsituOceanHadNMAT2(unittest.TestCase):

    def test_load_example(self):

        # example filename to use
        FORMAT_EXAMPLE = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/Docs/format.txt')
        FILENAME_EXAMPLE = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/ascii/Corrected/193305.txt.gz')

        # attempt load format
        format = HadNMAT2Format(FORMAT_EXAMPLE)

        # attempt load
        result = ObservationSourceInsituOceanHadNMAT2(format, FILENAME_EXAMPLE)

        # only supports tmean
        self.assertEqual([ ObservationSource.TMEAN ], result.observables())

        # get number of observations
        self.assertEqual(18821, result.number_of_observations())
        
        # get coordinates
        coordinates = result.observation_location_lookup()

        # should be 3 rows with a column per observation
        self.assertTrue(isinstance(coordinates, numpy.ndarray))
        self.assertEqual((2, 18821), coordinates.shape)

        # first two values and last value
        numpy.testing.assert_almost_equal([ 51.4,  2.4 ], coordinates[:,0].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([ -6.5, 74.5 ], coordinates[:,1].transpose(), decimal=4)
        numpy.testing.assert_almost_equal([ 39.9,-67.8 ], coordinates[:,18820].transpose(), decimal=4)

        # check shape and type of observations
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertTrue(isinstance(tmean, Observations))
        self.assertTrue(isinstance(tmean.mask, numpy.ndarray))
        self.assertTrue(tmean.mask.dtype == numpy.bool)
        self.assertTrue(isinstance(tmean.measurement, numpy.ndarray))
        self.assertTrue(tmean.measurement.dtype == numpy.float32)
        self.assertEqual((18821,), tmean.measurement.shape)
        self.assertTrue(isinstance(tmean.uncorrelatederror, numpy.ndarray))
        self.assertTrue(tmean.uncorrelatederror.dtype == numpy.float32)
        self.assertEqual((18821,), tmean.uncorrelatederror.shape)
        self.assertTrue(isinstance(tmean.locallycorrelatederror, list))
        self.assertEqual(0, len(tmean.locallycorrelatederror))

        # mask should be all in view
        self.assertTrue(numpy.all(numpy.logical_not(tmean.mask)))

        # times (first two are same, last is later)
        self.assertAlmostEqual(30435.0, tmean.time[0], places=4)
        self.assertAlmostEqual(30435.0, tmean.time[1], places=4)
        self.assertAlmostEqual(30455.0, tmean.time[18820], places=4)

        # locations should be just contiguous
        numpy.testing.assert_equal(range(18821), tmean.location)

        # first two values and last value
        self.assertAlmostEqual(281.49, tmean.measurement[0], places=4)
        self.assertAlmostEqual(300.70, tmean.measurement[1], places=4)
        self.assertAlmostEqual(292.12, tmean.measurement[18820], places=4)

        # at present the uncorrelated error is set to 1.1 everywhere
        self.assertTrue(numpy.all(tmean.uncorrelatederror == numpy.float32(1.1)))

        # check empty array of correlation length scales
        length_scale = result.local_correlation_length_scale(ObservationSource.TMEAN)
        self.assertTrue(isinstance(length_scale, numpy.ndarray))
        self.assertEqual(0, len(length_scale))

    def test_load_example_with_hash(self):

        # Some files have a hashtag in the station ID
        # - must check this works and isn't interpreted as a comment character

        FORMAT_EXAMPLE = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/Docs/format.txt')
        FILENAME_EXAMPLE = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/ascii/Corrected/188709.txt.gz')
        result = ObservationSourceInsituOceanHadNMAT2(HadNMAT2Format(FORMAT_EXAMPLE), FILENAME_EXAMPLE)
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertAlmostEqual(277.02, tmean.measurement[7941], places=4)

    def test_load_example_with_blank_ids(self):

        # As a quick implementation the ID column was originally treated as an integer (and hence ignored because it has text),
        # but this causes some files to fail where the station ID text happens to look like NaN floating point
        # and a conversion to integer is attempted.  Checking the problematic files:

        FORMAT_EXAMPLE = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/Docs/format.txt')
        FILENAME_EXAMPLE = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/ascii/Corrected/198711.txt.gz')
        result = ObservationSourceInsituOceanHadNMAT2(HadNMAT2Format(FORMAT_EXAMPLE), FILENAME_EXAMPLE)
        tmean = result.observations(ObservationSource.TMEAN)
        self.assertAlmostEqual(287.25, tmean.measurement[0], places=4)
