"""Tests for lattice."""

import unittest
import numpy

from eustace.analysis.advanced_standard.stats.spde.lattice import LatticeMesh, LatticeSPDE, WendlandC4Basis

class TestLatticeMesh(unittest.TestCase):
    
    def test_locate_nearby_basis_indices(self):
        """ Test that specified indices are identified by locate_nearby_basis_indices
        
        Allows for additional redundent locations that are not included 
        in the expected test indices to also be identified .
        
        """
        
        obs_locations = numpy.array([[-0.5, 0.0, 0.5, 1.5, 8.5, 9.0, 9.5]]).T
        
        # Non wrapped dimension test
        lattice_1d = LatticeMesh( dimension_specification = [(0., 9., 10),],
                               basis_function=WendlandC4Basis(),
                               overlap_factor=2.51,
                               wrap_dimensions = None )
    
        obs_map, node_map = lattice_1d.locate_nearby_basis_indices(obs_locations)
    
        expected_output = [[0, 0],
                           [0, 1],
                           [0, 2],
                           [1, 0],
                           [1, 1],
                           [1, 2],
                           [2, 0],
                           [2, 1],
                           [2, 2],
                           [2, 3],
                           [3, 0],
                           [3, 1],
                           [3, 2],
                           [3, 3],
                           [3, 4],
                           [4, 6],
                           [4, 7],
                           [4, 8],
                           [4, 9],
                           [5, 7],
                           [5, 8],
                           [5, 9],
                           [6, 7],
                           [6, 8],
                           [6, 9]]
    
        output_mapping = numpy.vstack((obs_map, node_map)).T

        self.assertTrue( numpy.all([output in output_mapping for output in expected_output]) )

        
        # Wrapped dimension test
        lattice_1d_wrapped = LatticeMesh( dimension_specification = [(0., 10., 10),],
                                       basis_function=WendlandC4Basis(),
                                       overlap_factor=2.5,
                                       wrap_dimensions = [True,] )
                                       
        obs_map, node_map = lattice_1d_wrapped.locate_nearby_basis_indices(obs_locations)
        
        expected_output = [[0, 7],
                           [0, 8],
                           [0, 9],
                           [0, 0],
                           [0, 1],
                           [0, 2],
                           [1, 8],
                           [1, 9],
                           [1, 0],
                           [1, 1],
                           [1, 2],
                           [2, 8],
                           [2, 9],
                           [2, 0],
                           [2, 1],
                           [2, 2],
                           [2, 3],
                           [3, 9],
                           [3, 0],
                           [3, 1],
                           [3, 2],
                           [3, 3],
                           [3, 4],
                           [4, 6],
                           [4, 7],
                           [4, 8],
                           [4, 9],
                           [4, 0],
                           [4, 1],
                           [5, 7],
                           [5, 8],
                           [5, 9],
                           [5, 0],
                           [5, 1],
                           [6, 7],
                           [6, 8],
                           [6, 9],
                           [6, 0],
                           [6, 1],
                           [6, 2]]
                           
        self.assertTrue( numpy.all([output in output_mapping for output in expected_output]) )

    def test_build_A_1d(self):
        pass
    
    def test_build_A_2d(self):
        pass

    def test_build_A_1d_wrapping(self):
        pass

    def test_build_A_2d_wrapping(self):
        pass

class TestLatticeSPDE(unittest.TestCase):

    
    def test_derivative(self):

        spde = LatticeSPDE.construct(
                dimension_specification = [(5, 20, 15)],
                basis_function=WendlandC4Basis(),
                overlap_factor=2.5)

        A = spde.build_A(numpy.array( [ 10.2 ] ) )
        self.assertEqual((1,15), A.shape)

        Q = spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0), 2, 1.01, 'csc')
        self.assertEqual((15,15), Q.shape)

        epsilon = 0.000001

        expected_dQ0 = \
            ((spde.build_Q_stationary(numpy.log(1.0)+epsilon, numpy.log(1.0), 2, 1.01, 'csc') - 
              spde.build_Q_stationary(numpy.log(1.0)-epsilon, numpy.log(1.0), 2, 1.01, 'csc')) / (2.0 * epsilon)).todense()

        expected_dQ1 = \
            ((spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0)+epsilon, 2, 1.01, 'csc') - 
              spde.build_Q_stationary(numpy.log(1.0), numpy.log(1.0)-epsilon, 2, 1.01, 'csc')) / (2.0 * epsilon)).todense()

        dQ0 = spde.build_dQdp_stationary(numpy.log(1.0), numpy.log(1.0), 2, 1.01, 0, 'csc').todense()
        dQ1 = spde.build_dQdp_stationary(numpy.log(1.0), numpy.log(1.0), 2, 1.01, 1, 'csc').todense()
        
        numpy.testing.assert_almost_equal(dQ0, expected_dQ0, decimal=7)
        numpy.testing.assert_almost_equal(dQ1, expected_dQ1, decimal=7)
