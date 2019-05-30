"""Quick tests for extended cholmod wrapper features."""

from eustace.analysis.advanced_standard.linalg.pardisowrapper import PardisoWrapper

import unittest
import numpy
from scipy.sparse import csr_matrix

class TestPardisoWrapper(unittest.TestCase):

    def setUp(self):

        # -- Simulate problem

        # lower triangular 4x4 matrix
        L = csr_matrix(
            [ [ 2.5, 0.0, 0.0, 0.0 ],
              [ 0.3, 1.4, 0.0, 0.0 ],
              [ 0.0, 0.0, 3.8, 0.0 ],
              [ 0.2, 0.0, 0.0, 2.2 ] ] )

        # corresponding A = L . Lt
        A = L * L.T

        # upper triangular bit on its own
        self.A = A

        # simulate solution that had two right-hand-sides
        self.x = numpy.array([ [ 1.0,  900.0 ],
                               [ 2.0, -230.0 ],
                               [ 3.0,   22.2 ],
                               [ 4.0, -250.0 ] ])

        # simulate the right-hand sides
        self.b = A * self.x             # A.x  = b
        self.c = L * self.x             # L.x  = c
        self.d = L.T * self.x           # Lt.x = d

    def test_base(self):

        PardisoWrapper.cholesky(self.A)

    def test_solve_A(self):

        wrapper = PardisoWrapper.cholesky(self.A)
        numpy.testing.assert_almost_equal(self.x, wrapper.solve_A(self.b))
        
    def test_solve_backward_substitution(self):
      
        wrapper = PardisoWrapper.cholesky(self.A)
        numpy.testing.assert_almost_equal(self.x, wrapper.solve_backward_substitution(self.d))

    def test_solve_forward_substitution(self):
      
        wrapper = PardisoWrapper.cholesky(self.A)
        numpy.testing.assert_almost_equal(self.x, wrapper.solve_forward_substitution(self.c))

    def test_diag_inv_A(self):
	
	wrapper = PardisoWrapper.cholesky(self.A)
        expected_diagonal_inv = numpy.diagonal(numpy.linalg.inv(self.A.todense()))
	
	numpy.testing.assert_almost_equal(wrapper.diag_inv_A(max_blocksize = None), expected_diagonal_inv)
	numpy.testing.assert_almost_equal(wrapper.diag_inv_A(max_blocksize = 1), expected_diagonal_inv)
	numpy.testing.assert_almost_equal(wrapper.diag_inv_A(max_blocksize = 2), expected_diagonal_inv)
	numpy.testing.assert_almost_equal(wrapper.diag_inv_A(max_blocksize = 3), expected_diagonal_inv)
	numpy.testing.assert_almost_equal(wrapper.diag_inv_A(max_blocksize = 4), expected_diagonal_inv)
	numpy.testing.assert_almost_equal(wrapper.diag_inv_A(max_blocksize = 100), expected_diagonal_inv)

