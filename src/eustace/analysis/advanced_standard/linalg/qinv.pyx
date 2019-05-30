import numpy as np
cimport numpy as np
import scipy.sparse
import sksparse.cholmod as cholmod

DTYPE =  np.float64
ctypedef np.int32_t cITYPE
ctypedef np.float64_t cDTYPE

def Q_inv(Q, fQ = None):
    """Compute partial inverse of Q
    
    Compute the values of the inverse of Q only at only the elements where the
    Q is non zero. Other elements are returned as zero.
    
    Marginal variances can be extracted from the diagonal of Q_inv(Q)
    where Q is a precision matrix. This can also useful for calculation of
    the trace of products such as trace( Q^-1 * A ) where the second 
    term has a non-zero pattern that does not require the full inv(Q). Typically
    this is useful in cases where A has the same sparsity pattern has a graph that
    is a sub-graph of that of Q, for example when computing trace( Q^-1 * dQdp )
    where dQdp is the derivative of Q with respect to some parameter of the elements
    of Q.
    
    Note that this function does *not* return the full inverse of Q.
    
    Basic recursion pseudo-code:
    
    for i = n-1, ..., 0:
        for j = n-1, ..., i:
            if (L[j,i] not known to be 0):
                S[j,i] = S[i,j] = (I(i==j)/L[i,i] - sum_{k=i+1}^{n-1} L[k,i] S[k,j] ) / L[i,i]
    
    """
    
    if fQ is None:
        fQ = cholmod.cholesky(Q)
    
    cdef cITYPE dim = Q.shape[0]
    
    Lp = fQ.L().tocsc() # get cholesky triangle as convert it to csc format
    
    if not Lp.has_sorted_indices:
        Lp.sort_indices()
    
    p  = fQ.P()
    P = scipy.sparse.coo_matrix((np.ones(len(p)), (np.arange(len(p)),p))).tocsr()
    
    # For csc matrix, item (i, j) can be accessed as data[indptr[j]+k] where k is
    # the position of i in the indices variable at indices[indptr[j]:indptr[j+1]]
    cdef cDTYPE [:] L_data    = Lp.data
    cdef cITYPE [:] L_indices = Lp.indices
    cdef cITYPE [:] L_indptr  = Lp.indptr

    # initialise output non-zero pattern based on that of the cholesky factor
    Sigma_tilde = Lp + Lp.T
    if not Sigma_tilde.has_sorted_indices:
        Sigma_tilde.sort_indices()
    
    cdef cDTYPE [:] S_data = np.zeros(Sigma_tilde.data.shape)
    cdef cITYPE [:] S_indices = Sigma_tilde.indices
    cdef cITYPE [:] S_indptr = Sigma_tilde.indptr
    
    cdef cITYPE [:] S_col_I, S_col_J, L_col_I
    
    cdef cITYPE i, j
    cdef cITYPE L_col_I_index, S_col_I_index, S_col_J_index
    
    # direct indices to data array for indicated columns
    cdef cITYPE L_data_I_index, S_data_I_index, S_data_J_index 
    
    cdef cDTYPE value
    
    for i in range(dim-1,-1,-1):
        
        # step backwards through non-missing indices in column i of S
        for S_col_I_index in range(S_indptr[i+1]-S_indptr[i]-1, -1, -1):
            # computing S element (i,j) and also filling (j,i)
            S_data_I_index = S_indptr[i] + S_col_I_index
            j = S_indices[S_data_I_index]
            
            if j < i:
                continue
            
            value = 0.0
            
            # start loop from last element in columns i and j and decrement finding matching i and j, terminating when i == 0
            L_data_I_index = L_indptr[i+1] - 1
            S_data_J_index = S_indptr[j+1] - 1
            
            while L_indices[L_data_I_index] > i:

                # decrement S_col_J_index until S_col_J[S_col_J_index] = L_col_I[L_col_I_index]
                while (S_data_J_index > 0 and (L_indices[L_data_I_index] < S_indices[S_data_J_index])):
                    S_data_J_index-= 1
            
                if L_indices[L_data_I_index] == S_indices[S_data_J_index]:
                    
                    value -= L_data[L_data_I_index] * S_data[S_data_J_index] 
                    S_data_J_index -= 1
                    
                L_data_I_index -= 1
                
            if i ==j:

                value += 1.0 / L_data[L_data_I_index]
                value /= L_data[L_data_I_index]

                S_data[S_data_I_index] = value
            else:

                value /= L_data[L_data_I_index]

                S_data[S_data_I_index] = value
            
            # now fill in S[j,i] in upper triangle
            while S_indices[S_data_J_index] > i:
                
                S_data_J_index -= 1

            S_data[S_data_J_index] = value

    Sigma_tilde = scipy.sparse.csc_matrix( (S_data, S_indices, S_indptr) )
    
    return P.T * Sigma_tilde.tocsc() * P
