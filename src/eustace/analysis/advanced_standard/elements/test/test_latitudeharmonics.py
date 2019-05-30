"""Tests for latitude harmonics element."""

import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeFunction
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsPrior
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElementDesign
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters

class TestLatitudeFunction(unittest.TestCase):
    """Test the generic function wrapper used by latitude element."""

    def test_init(self):

        f = LatitudeFunction(numpy.cos, 4.0)
        self.assertEqual(numpy.cos, f.numpymethod)

    def test_compute(self):

        numpy.testing.assert_almost_equal([ [ 0.0 ], [ -0.5 ] ], LatitudeFunction(numpy.sin, 2.0).compute( numpy.array([ -90.0, -15.0 ])))

class TestLatitudeHarmonicsElement(unittest.TestCase):

    class SimulatedObservationStructure(ObservationStructure):
        """Simulated data with 3 observations at different locations."""
                
        def location_polar_coordinates(self):
            """Array of [ [ latitude, longitude ] ] for each observation."""

            return numpy.array( [ 
                [ 15.0, -7.0 ],
                [  5.0, 100.0 ],
                [ -2.0, 103.2] ])
          
    def test_init(self):

        element = LatitudeHarmonicsElement()
        self.assertFalse(element.isnonlinear())

    def test_element_design(self):

        element = LatitudeHarmonicsElement()
        design = element.element_design(TestLatitudeHarmonicsElement.SimulatedObservationStructure())
        self.assertIsInstance(design, LatitudeHarmonicsElementDesign)
        self.assertFalse(design.isnonlinear())

    def test_element_prior(self):

        element = LatitudeHarmonicsElement()
        prior = element.element_prior(
            CombinationHyperparameters([ CovariateHyperparameters(c) for c in [ 1.0, 1.1, 1.2, 1.3 ]]))
        self.assertIsInstance(prior, LatitudeHarmonicsPrior)
        
        # As an example check precision - should be diagonals with exp(-2 x hyperparameter)
        precision = prior.prior_precision()
        self.assertEqual(SPARSEFORMAT, precision.getformat())
        self.assertEqual(4, precision.nnz)
        self.assertAlmostEqual(numpy.exp(-2.0), precision[0,0])
        self.assertAlmostEqual(numpy.exp(-2.2), precision[1,1])
        self.assertAlmostEqual(numpy.exp(-2.4), precision[2,2])
        self.assertAlmostEqual(numpy.exp(-2.6), precision[3,3])

class TestLatitudeHarmonicsElementDesign(unittest.TestCase):

    def test_init(self):

        design = LatitudeHarmonicsElementDesign(TestLatitudeHarmonicsElement.SimulatedObservationStructure())
        numpy.testing.assert_almost_equal([ 15.0, 5.0, -2.0 ], design.latitudes)

    def test_design_number_of_state_parameters(self):

        design = LatitudeHarmonicsElementDesign(TestLatitudeHarmonicsElement.SimulatedObservationStructure())
        self.assertEqual(4, design.design_number_of_state_parameters())
    
    def test_design_matrix(self):

        design = LatitudeHarmonicsElementDesign(TestLatitudeHarmonicsElement.SimulatedObservationStructure())
        m = design.design_matrix()
        self.assertEqual(SPARSEFORMAT, m.getformat())
        numpy.testing.assert_almost_equal(
            [ [ 0.86602540, 0.5,        0.5,        0.86602540 ],
              [ 0.98480775, 0.17364818, 0.93969262, 0.34202014 ],
              [ 0.99756405,-0.06975647, 0.99026807,-0.13917310 ] ],
            m.todense())
  