import unittest
import numpy
import scipy.sparse

from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.local_view import LocalSubRegion
from eustace.analysis.advanced_standard.elements.local import LocalDesign
from eustace.analysis.advanced_standard.elements.local import LocalPrior
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
               
class TestLocalSubRegion(unittest.TestCase):
    
    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 15.0, -7.0 ], [  5.0, 100.0 ] ])                 
    
    def test_init(self):
        
        n_triangulation_divisions = 2
        neighbourhood_level = 1
        centre_index_at_level = 0
        localsubregion = LocalSubRegion(n_triangulation_divisions, neighbourhood_level, centre_index_at_level)
        self.assertEqual(2, localsubregion.spde.triangulation.level)

    def test_element_design(self):
        
        n_triangulation_divisions = 2
        neighbourhood_level = 1
        centre_index_at_level = 0
        localsubregion = LocalSubRegion(n_triangulation_divisions, neighbourhood_level, centre_index_at_level)
        design = localsubregion.element_design(TestLocalSubRegion.SimulatedObservationStructure())
        self.assertTrue(isinstance(design, LocalDesign))
        self.assertEqual(n_triangulation_divisions, design.spde.triangulation.level)

    def test_design_prior(self):
        
        n_triangulation_divisions = 2
        neighbourhood_level = 1
        centre_index_at_level = 0
        localsubregion = LocalSubRegion(n_triangulation_divisions, neighbourhood_level, centre_index_at_level)
        prior = localsubregion.element_prior(LocalHyperparameters(0.0, numpy.log(10.0 * numpy.pi/180)))
        self.assertTrue(isinstance(prior, LocalPrior))
        numpy.testing.assert_almost_equal(prior.hyperparameters.get_array(), [ 0.0, numpy.log(10.0 * numpy.pi/180) ])
        self.assertEqual(n_triangulation_divisions, prior.spde.triangulation.level)
        self.assertEqual(2, prior.alpha)