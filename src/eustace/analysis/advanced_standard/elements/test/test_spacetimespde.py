"""Tests for base classes of SPDE space-time approach."""

import unittest
import numpy
from eustace.analysis.advanced_standard.elements.spacetimespde import SpaceTimeSPDEElement
from eustace.analysis.advanced_standard.elements.spacetimespde import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.spacetimespde import SpaceTimeSPDEPrior
from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.advanced_standard.stats.spde.lattice import LatticeSPDE, WendlandC4Basis

class TestSpaceTimeSPDEHyperparameters(unittest.TestCase):

    def test_init(self):
        
        p = SpaceTimeSPDEHyperparameters(space_log_sigma=0.2, space_log_rho=0.8, time_log_rho=0.7)
        self.assertEqual(0.2, p.space_log_sigma)
        self.assertEqual(0.8, p.space_log_rho)
        self.assertEqual(0.7, p.time_log_rho)

    def test_get_array(self):

        numpy.testing.assert_equal([ 0.3, 0.4, 0.5 ], SpaceTimeSPDEHyperparameters(space_log_sigma=0.3, space_log_rho=0.4, time_log_rho=0.5).get_array())

    def test_set_array(self):

        p = SpaceTimeSPDEHyperparameters(space_log_sigma=0.3, space_log_rho=0.4, time_log_rho=0.5)
        p.set_array(numpy.array([ 23.2, 22.1, 20.5 ]))
        self.assertEqual(20.5, p.time_log_rho)
        self.assertEqual(23.2, p.space_log_sigma)
        self.assertEqual(22.1, p.space_log_rho)

class TestSpaceTimeSPDElement(unittest.TestCase):
 
    def test_init(self):
        spde_element = SpaceTimeSPDEElement(n_triangulation_divisions=3, alpha=2, starttime=0., endtime=10., n_nodes=11., overlap_factor=1, H=1)
        
        self.assertTrue(isinstance(spde_element.spatial_model, SphereMeshSPDE))
        self.assertEqual(spde_element.alpha, 2)
        self.assertTrue(isinstance(spde_element.temporal_model, LatticeSPDE))
        self.assertEqual(spde_element.H, 1)

        self.assertEqual(1, spde_element.temporal_model.lattice.n_dims)
        numpy.testing.assert_equal([[0.],[10.]], spde_element.temporal_model.lattice.node_bounds)
        self.assertEqual([11], spde_element.temporal_model.lattice.number_of_nodes_per_dimension)
        self.assertEqual([1], spde_element.temporal_model.lattice.node_spacing)

class TestSpaceTimeSPDEPrior(unittest.TestCase):

    def setUp(self):

        self.prior = SpaceTimeSPDEPrior(
            hyperparameters=SpaceTimeSPDEHyperparameters(numpy.log(1.0), numpy.log(1.1), numpy.log(1.2)),
            spatial_model=SphereMeshSPDE(level=1),
            alpha=2,
            temporal_model=LatticeSPDE.construct(
                dimension_specification = [(23, 27, 5)],
                basis_function=WendlandC4Basis(),
                overlap_factor=2.5),
            H=1.01)

    def test_init(self):

        numpy.testing.assert_equal(self.prior.hyperparameters.get_array(), [ numpy.log(1.0), numpy.log(1.1), numpy.log(1.2) ])
        self.assertEqual(1, self.prior.spatial_model.triangulation.level)
        self.assertEqual(2, self.prior.alpha)
        numpy.testing.assert_almost_equal(self.prior.temporal_model.lattice.axis_coordinates, [ [ 23.0 ], [ 24.0 ], [ 25.0 ], [ 26.0 ], [ 27.0 ] ])
        self.assertEqual(2.5, self.prior.temporal_model.lattice.basis_function.basis_span)
        self.assertEqual(1.01, self.prior.H)
