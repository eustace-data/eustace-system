"""Sparse matrix operations"""

import numpy
import qinv

def trace_iQ_dQdp(Q, dQdp, fQ = None, iQp = None):
    """Compute trace( inverse(Q) * dQ/dp )
    
    Uses the property that a partial inverse of Q
    based on the Takahasi recursions has a non-
    zero pattern that encompases the non-zero
    pattern in dQdp.
    
    Args:
        * Q:
            Sparse matrix
        * dQdp:
            Sparse matrix that is the derivative of Q with respect to its
            parameters. Must have a non-zero pattern contained within that
            of Q.
    Kwargs:
        * fQ:
            Optional pre-factorised Q object returned from sksparse.cholesky(Q) 
        * iQp:
            Optional pre-computed partial inverse of Q.
    
    """
    
    if iQp is None:
        if fQ is None:
            fQ = cholmod.cholesky(Q)
        iQp = qinv.Q_inv(Q, fQ)
    
    trace = numpy.sum( (iQp * dQdp).diagonal() )
    
    return trace
    
def marginal_variance_Q_inv(Q, A = None):
    """Use partial inverse to compute the marginal variance"""
    
    Q_inv = qinv.Q_inv(Q)
    
    if A is not None:
        # The following operation is only correct if the non zero indices of A.T * A are a subset of
        # the non-zero elements of Q_inv so check that this condition is true before proceeding. 
        
        from sets import Set
        
        (P_i, P_j, P_v) = scipy.sparse.find(Q_inv)
        (AtA_i, AtA_j, AtA_v) = scipy.sparse.find(A.T * A)
        Q_inv_non_zero = Set( zip(P_i, P_j) )
        AtA_non_zero = Set( zip(AtA_i, AtA_j) )
        
        AtA_non_zero.issubset(Q_inv_non_zero)
        
        if AtA_non_zero.issubset(Q_inv_non_zero):
            return (A * Q_inv * A.T).diagonal()
        else:
            # Note that Q_inv could be extended to compute values at the non-zero elements of A too
            raise Exception('Marginal variance is only valid if the non-zero elements in A are a subset of those in Q_inv')
    
    else:
        return Q_inv.diagonal()