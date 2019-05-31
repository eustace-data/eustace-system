import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.latitudespline import PADDING_NODES
from eustace.analysis.advanced_standard.elements.latitudespline import LatitudeSplineElement
from eustace.analysis.advanced_standard.elements.latitudespline import LatitudeSplineElementDesign
from eustace.analysis.advanced_standard.elements.latitudespline import LatitudeSplineElementPrior
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT

class TestLatitudeSplineElement(unittest.TestCase):
    
    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 15.0, -7.0 ], [  5.0, 100.0 ] ])                 
    
    def test_init(self):
        
        alpha = 2
        n_nodes = 36
        overlap_factor = 2.5
        H = 1.0
        
        padded_nodes = n_nodes + 2 * PADDING_NODES 
        
        element = LatitudeSplineElement(alpha, n_nodes, overlap_factor, H)
        
        self.assertEqual(2, element.alpha)
        self.assertEqual(padded_nodes, element.spde.n_latent_variables())
        self.assertEqual(1.0, element.H)

    def test_element_design(self):

        alpha = 2
        n_nodes = 36
        overlap_factor = 2.5
        H = 1.0
        
        padded_nodes = n_nodes + 2 * PADDING_NODES 
        
        element = LatitudeSplineElement(alpha, n_nodes, overlap_factor, H)

        design = element.element_design(TestLatitudeSplineElement.SimulatedObservationStructure())
        self.assertTrue(isinstance(design, LatitudeSplineElementDesign))
        self.assertEqual(padded_nodes, design.design_number_of_state_parameters())

    def test_design_prior(self):

        alpha = 2
        n_nodes = 36
        overlap_factor = 2.5
        H = 1.0
        
        padded_nodes = n_nodes + 2 * PADDING_NODES 
        
        element = LatitudeSplineElement(alpha, n_nodes, overlap_factor, H)
        
        prior = element.element_prior(LocalHyperparameters(0.0, numpy.log(10.0 * numpy.pi/180)))
        self.assertTrue(isinstance(prior, LatitudeSplineElementPrior))
        numpy.testing.assert_almost_equal(prior.hyperparameters.get_array(), [ 0.0, numpy.log(10.0 * numpy.pi/180) ])
        self.assertEqual(padded_nodes, prior.prior_number_of_state_parameters())
        self.assertEqual(2, prior.alpha)

