"""Tests for space-time factor analysis."""

import unittest
import numpy
import scipy.sparse
from datetime import datetime

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.mesh.geometry import cartesian_to_polar2d
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.spacetimespde import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeFactorElement
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeFactorDesign
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeFactorPrior

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.advanced_standard.stats.spde.lattice import LatticeSPDE, WendlandC4Basis

class TestSpaceTimeFactorElement(unittest.TestCase):

    class SimulatedObservationStructure(ObservationStructure):
                
        def location_polar_coordinates(self):

            return numpy.array( [ [ 15.0, -7.0 ], [  5.0, 100.0 ] ])      

        def time_datetime(self):

            return datetime(1850, 1, 25)

        def time_index(self):

            return 24

    def test_init(self):

        element = SpaceTimeFactorElement(n_triangulation_divisions=3, alpha=2, starttime=23, endtime=27, n_nodes=5, overlap_factor=2.5, H=1.2)
        self.assertEqual(3, element.spatial_model.triangulation.level)
        self.assertEqual(2, element.alpha)
        numpy.testing.assert_almost_equal(element.temporal_model.lattice.axis_coordinates, [ [ 23.0 ], [ 24.0 ], [ 25.0 ], [ 26.0 ], [ 27.0 ] ])
        self.assertEqual(2.5, element.temporal_model.lattice.basis_function.basis_span)
        self.assertEqual(1.2, element.H)

    def test_element_design(self):

        element = SpaceTimeFactorElement(n_triangulation_divisions=1, alpha=2, starttime=23, endtime=27, n_nodes=5, overlap_factor=2.5, H=1.2)
        design = element.element_design(TestSpaceTimeFactorElement.SimulatedObservationStructure())
        self.assertTrue(isinstance(design, SpaceTimeFactorDesign))
        self.assertEqual(2, design.alpha)
        self.assertEqual(1.2, design.H)
        self.assertEqual(1, element.spatial_model.triangulation.level)
        self.assertEqual(2.5, element.temporal_model.lattice.basis_function.basis_span)

    def test_element_prior(self):

        element = SpaceTimeFactorElement(n_triangulation_divisions=1, alpha=2, starttime=23, endtime=27, n_nodes=5, overlap_factor=2.5, H=1.2)
        prior = element.element_prior(SpaceTimeSPDEHyperparameters(0.0, 1.0, 1.0))
        self.assertTrue(isinstance(prior, SpaceTimeFactorPrior))
        

class TestSpaceTimeFactorDesign(unittest.TestCase):

    def test_init(self):

        design = SpaceTimeFactorDesign(
            observationstructure=TestSpaceTimeFactorElement.SimulatedObservationStructure(),
            spatial_model=SphereMeshSPDE(level=1),
            alpha=2,
            temporal_model=LatticeSPDE.construct(
                dimension_specification = [(23, 27, 5)],
                basis_function=WendlandC4Basis(),
                overlap_factor=2.5),
            H=1.01)

        self.assertEqual(2, design.alpha)
        self.assertEqual(1.01, design.H)
        self.assertEqual(1, design.spatial_model.triangulation.level)
        self.assertEqual(2.5, design.temporal_model.lattice.basis_function.basis_span)
        numpy.testing.assert_almost_equal(design.temporal_model.lattice.axis_coordinates, [ [ 23.0 ], [ 24.0 ], [ 25.0 ], [ 26.0 ], [ 27.0 ] ])

    def test_design_state_parameters(self):

        obs = TestSpaceTimeFactorElement.SimulatedObservationStructure()

        spde_space = SphereMeshSPDE(level=1)

        spde_time = LatticeSPDE.construct(
            dimension_specification = [(23, 27, 5)],
            basis_function=WendlandC4Basis(),
            overlap_factor=2.5)

        design = SpaceTimeFactorDesign(
            observationstructure=obs,
            spatial_model=spde_space,
            alpha=2,
            temporal_model=spde_time,
            H=1.01)
        
        self.assertEqual(47, design.design_number_of_state_parameters())       
        numpy.testing.assert_equal(numpy.zeros(42,), design.design_initial_state()[0:42])
        numpy.testing.assert_equal(numpy.ones(5,), design.design_initial_state()[42:])

    def test_design_function(self):

        obs = TestSpaceTimeFactorElement.SimulatedObservationStructure()

        spde_space = SphereMeshSPDE(level=1)

        spde_time = LatticeSPDE.construct(
            dimension_specification = [(23, 27, 5)],
            basis_function=WendlandC4Basis(),
            overlap_factor=2.5)

        design = SpaceTimeFactorDesign(
            observationstructure=obs,
            spatial_model=spde_space,
            alpha=2,
            temporal_model=spde_time,
            H=1.01)

        # Get and check indices of nonzero design elements
        A_space = spde_space.build_A(obs.location_polar_coordinates())
        observation_indices, vertex_indices = A_space.sorted_indices().nonzero()
        numpy.testing.assert_equal(observation_indices, [ 0, 0, 0, 1, 1, 1 ])
        self.assertEqual((6,), vertex_indices.shape)

        # Vertices of observations 0 and 1
        # (one row per vertex, one column per coordinate)
        vertices0 = spde_space.triangulation.points[vertex_indices[0:3],:]
        vertices1 = spde_space.triangulation.points[vertex_indices[3:6],:]

        # Multiply vertices by the weights from A and sum to get cartesian locations
        testpoint0 = A_space[0, vertex_indices[0:3] ] * vertices0
        testpoint1 = A_space[1, vertex_indices[3:6] ] * vertices1

        # Check results correspond to original polar coordinates
        numpy.testing.assert_almost_equal(cartesian_to_polar2d(testpoint0), [ [ 15.0,  -7.0 ] ])
        numpy.testing.assert_almost_equal(cartesian_to_polar2d(testpoint1), [ [  5.0, 100.0 ] ])

        # So the function with those indices set should give required values, provided it's evaluated at time index == 24
        state_vector = numpy.tile(777.7, (47,))
        state_vector[vertex_indices[0:3]] = 40.0
        state_vector[vertex_indices[3:6]] = 28.0
        state_vector[42:47] = [ 0.0, 0.75, 0.0, 0.0, 0.0 ]
        numpy.testing.assert_almost_equal(design.design_function(state_vector), [ 30.0 , 21.0 ])

        # But if we change to a time far away then we get nothing
        state_vector[42:47] = [ 0.0, 0.0, 0.0, 0.0, 111.1 ]
        numpy.testing.assert_almost_equal(design.design_function(state_vector), [ 0.0 , 0.0 ])
        
    def test_design_jacobian(self):

        obs = TestSpaceTimeFactorElement.SimulatedObservationStructure()

        spde_space = SphereMeshSPDE(level=1)

        spde_time = LatticeSPDE.construct(
            dimension_specification = [(23, 27, 5)],
            basis_function=WendlandC4Basis(),
            overlap_factor=2.5)

        design = SpaceTimeFactorDesign(
            observationstructure=obs,
            spatial_model=spde_space,
            alpha=2,
            temporal_model=spde_time,
            H=1.01)
    
        # Build a candidate state vector (as used above for function testing)
        A_space = spde_space.build_A(obs.location_polar_coordinates())
        observation_indices, vertex_indices = A_space.sorted_indices().nonzero()
        state_vector = numpy.tile(777.7, (47,))
        state_vector[vertex_indices[0:3]] = 40.0
        state_vector[vertex_indices[3:6]] = 28.0
        state_vector[42:47] = [ 0.0, 0.75, 0.0, 0.0, 0.0 ]

        # Numerical jacobian
        J_numerical = numpy.zeros((2, 47))
        epsilon = 0.01
        for parameter_index in range(47):

            x0 = numpy.copy(state_vector)
            x1 = numpy.copy(state_vector)
            x0[parameter_index] -= epsilon
            x1[parameter_index] += epsilon
            J_numerical[:, parameter_index] = (design.design_function(x1) - design.design_function(x0)) / (2.0 * epsilon)

        # Computed jacobian
        J_calculated = design.design_jacobian(state_vector)

        # Should be the same
        numpy.testing.assert_almost_equal(J_calculated.todense(), J_numerical)

        # And numbers of nonzeros also the same
        self.assertEqual(J_calculated.nnz, scipy.sparse.csc_matrix(J_numerical).nnz)

class TestSpaceTimeFactorPrior(unittest.TestCase):

    def setUp(self):

        self.prior = SpaceTimeFactorPrior(
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

    def test_prior_number_of_state_parameters(self):

        self.assertEqual(47, self.prior.prior_number_of_state_parameters())

    def test_prior_precision(self):

        Q = self.prior.prior_precision()

        # For now just check size and format of result
        # Should have dimension 47 (because spatial subdivision has dimension 42 and there are 5 time points)
        self.assertEqual(SPARSEFORMAT, Q.getformat())
        self.assertEqual((47,47), Q.shape)

    def test_prior_precision_derivative(self):

        dQ_0 = self.prior.prior_precision_derivative(0)
        dQ_1 = self.prior.prior_precision_derivative(1)
        dQ_2 = self.prior.prior_precision_derivative(2)

        # Check size and format of result
        # Should have dimension 47 (because spatial subdivision has dimension 42 and there are 5 time points)
        self.assertEqual(SPARSEFORMAT, dQ_0.getformat())
        self.assertEqual(SPARSEFORMAT, dQ_1.getformat())
        self.assertEqual(SPARSEFORMAT, dQ_2.getformat())
        self.assertEqual((47,47), dQ_0.shape)
        self.assertEqual((47,47), dQ_1.shape)
        self.assertEqual((47,47), dQ_2.shape)

        # Numerical derivative
        numerical = [ [ ], [ ], [ ] ]
        epsilon = 0.0001
        for parameter_index in range(3):

            for sign_index, sign in enumerate([ -epsilon, +epsilon ]):

                parameter_vector = numpy.array([numpy.log(1.0), numpy.log(1.1), numpy.log(1.2)])
                parameter_vector[parameter_index] += sign

                Q_numerical = SpaceTimeFactorPrior(
                    hyperparameters=SpaceTimeSPDEHyperparameters(parameter_vector[0], parameter_vector[1], parameter_vector[2]),
                    spatial_model=SphereMeshSPDE(level=1),
                    alpha=2,
                    temporal_model=LatticeSPDE.construct(
                        dimension_specification = [(23, 27, 5)],
                        basis_function=WendlandC4Basis(),
                        overlap_factor=2.5),
                    H=1.01).prior_precision()

                numerical[parameter_index].append(Q_numerical)

        numerical_dQ0 = ((numerical[0][1]) - (numerical[0][0])) / (2.0 * epsilon)
        numerical_dQ1 = ((numerical[1][1]) - (numerical[1][0])) / (2.0 * epsilon)
        numerical_dQ2 = ((numerical[2][1]) - (numerical[2][0])) / (2.0 * epsilon)

        # Check numerical derivative corresponds to computed one
        numpy.testing.assert_almost_equal(dQ_0.todense(), numerical_dQ0.todense())
        numpy.testing.assert_almost_equal(dQ_1.todense(), numerical_dQ1.todense())
        numpy.testing.assert_almost_equal(dQ_2.todense(), numerical_dQ2.todense(), decimal=7)
