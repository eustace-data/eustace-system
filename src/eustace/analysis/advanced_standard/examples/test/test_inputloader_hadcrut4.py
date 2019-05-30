"""Test reading of HadCRUT4 data."""

import unittest
import os.path
import numpy
from eustaceconfig import WORKSPACE_PATH
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.advanced_standard.examples.inputloader_hadcrut4 import AnalysisSystemInputLoaderHadCRUT4, ObservationSourceHadCRUT4

class TestInputLoaderHadCRUT4(unittest.TestCase):

    def test_load_source(self):

        # Input spec
        basepath = os.path.join(WORKSPACE_PATH, 'data/incoming/HadCRUT4.5.0.0')
        filenames = [
            'hadcrut4_median_netcdf.nc',
            'hadcrut4_uncorrelated_supplementary.nc',
            'hadcrut4_blended_uncorrelated.nc' ]

        monthnumber = 0

        # Configure load instance
        loader = AnalysisSystemInputLoaderHadCRUT4([ os.path.join(basepath, filename) for filename in filenames ])

        # Load it
        source = loader.load_observation_source(monthnumber)

        # Read obs
        obs = source.observations(ObservationSource.TMEAN)
        
        # Time just copied
        self.assertEqual(0, obs.time)

        # Check types
        self.assertIsInstance(obs.mask, numpy.ndarray)
        self.assertIsInstance(obs.location, numpy.ndarray)
        self.assertIsInstance(obs.measurement, numpy.ndarray)
        self.assertIsInstance(obs.uncorrelatederror, numpy.ndarray)
        self.assertIsNone(obs.locallycorrelatederror)

        # Check sizes
        self.assertEqual(2592, obs.number_of_observations())
        self.assertEqual((2592,), obs.mask.shape)
        self.assertEqual((2592,), obs.location.shape)
        self.assertEqual((2592,), obs.measurement.shape)
        self.assertEqual((2592,), obs.uncorrelatederror.shape)
        self.assertEqual(True, obs.mask[0])

        # Also read location lookup (should just be grid)
        lookup = source.observation_location_lookup()
        self.assertEqual((2, 2592), lookup.shape)
        numpy.testing.assert_almost_equal([ -87.5, -177.5 ], lookup[:, 0])
        numpy.testing.assert_almost_equal([ -87.5, -172.5 ], lookup[:, 1])
        numpy.testing.assert_almost_equal([ -82.5, -177.5 ], lookup[:, 72])
        numpy.testing.assert_almost_equal([  87.5,  177.5 ], lookup[:, 2591])
