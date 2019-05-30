import numpy as numpy
import scipy.sparse
import scipy.sparse.linalg
import sksparse.cholmod
from eustace.analysis.advanced_standard.linalg import qinv

import unittest

class TestPartialInverse(unittest.TestCase):
    
    def test_qinv_dense(self):
        """Check that qinv returns the correct full inverse for a dense Q matrix"""
        
        dim = 6
        nsamples = 30
        s = numpy.random.normal(0.0, 1.0, (dim,nsamples))
        S = numpy.dot(s, s.T)
        Q = scipy.sparse.csc_matrix( scipy.sparse.linalg.spsolve(S, numpy.eye(dim) ) )
        
        S_partial = qinv.Q_inv( Q )
        
        numpy.testing.assert_almost_equal(S_partial.todense(), S)

    def test_qinv_eye(self):
        """Check that qinv returns the correct diagonal for inverse of diagonal Q matrix"""
        
        dim = 5
        
        Q = scipy.sparse.eye( dim )
        S = scipy.sparse.linalg.spsolve(Q, numpy.eye(dim) )
        
        S_partial = qinv.Q_inv( Q )
        
        numpy.testing.assert_almost_equal(S_partial.diagonal(), S.diagonal())
    
    def test_qinv_sparse(self):
        """Check that non-zero elements of partial inverse match the full inverse at those elements"""
        
        Q = scipy.sparse.csc_matrix( [[3, -1, 0, 0, 0],
                                      [-1, 3, -1, 0, 0],
                                      [0, -1, 3, -1, 0],
                                      [0, 0, -1, 3, -1],
                                      [0, 0, 0, -1, 3],] )
        
        S = scipy.sparse.linalg.spsolve(Q, numpy.eye(Q.shape[0]) )
        
        S_partial = qinv.Q_inv( Q )        
        non_zero_indices = S_partial.nonzero()
        
        numpy.testing.assert_almost_equal(numpy.squeeze(numpy.asarray(S_partial[non_zero_indices])), \
                                          S[non_zero_indices])

    def test_qinv_sparsity(self):
        """Check that partial inverse has sparsity structure defined by cholesky factor.
        
        For the partial inverse, checks that:
        
            * the lower triangle of the (reordered) partial inverse has a sparsity pattern
              that matches that of L.
            * the upper triangle has a sparsity pattern matching that of L.T.
            * the partial inverse is symetric with values that match the full 
              inverse at indices where the partial inverse is non-zero (other
              indices should not be expected to match).
        
        """
        
        
        Q = scipy.sparse.csc_matrix( [[3, -1, 0, 0, 0],
                                      [-1, 3, -1, 0, 0],
                                      [0, -1, 3, -1, 0],
                                      [0, 0, -1, 3, -1],
                                      [0, 0, 0, -1, 3],] )
        
        fQ = sksparse.cholmod.cholesky( Q )

        # Cholesky factor and matrix reordering
        L = fQ.L()
        p = fQ.P()
        P = scipy.sparse.coo_matrix((numpy.ones(len(p)), (numpy.arange(len(p)),p))).tocsc()

        
        S_partial = qinv.Q_inv( Q )
        
        # Undo reordering before comparison of non-zero patter with that of  L
        S_partial_rotated = P * S_partial * P.T 

        # check lower triangle sparsity
        S_partial_lower_nonzero_inds = numpy.sort( numpy.array( scipy.sparse.tril(S_partial_rotated).nonzero() ) )
        L_nonzero_inds =  numpy.sort( numpy.array( L.nonzero() ) )        
        numpy.testing.assert_equal(S_partial_lower_nonzero_inds, L_nonzero_inds)
        
        # check upper triangle sparsity
        S_partial_upper_nonzero_inds = numpy.sort( numpy.array( scipy.sparse.triu(S_partial_rotated).nonzero() ) )
        Lt_nonzero_inds =  numpy.sort( numpy.array( L.T.nonzero() ) )        
        numpy.testing.assert_equal(S_partial_upper_nonzero_inds, Lt_nonzero_inds)
        
        # check symmetry
        numpy.testing.assert_almost_equal(S_partial.todense(), S_partial.T.todense())
    
        # check partial inverse correctness at non-zero locations
        S_full = fQ.solve_A( scipy.sparse.eye( Q.shape[0] ) )
        
        non_zero_indices = S_partial.nonzero()
        
        numpy.testing.assert_almost_equal(S_partial[non_zero_indices], \
                                          S_full[non_zero_indices])