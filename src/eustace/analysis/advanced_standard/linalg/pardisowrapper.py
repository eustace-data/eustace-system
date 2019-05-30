"""
Work directly with Intel MKL PARDISO library.

This Python file is inspired by code from Adrian Haas and ETH Zurich,
though it has been modified a lot from the original.
A copy of the license associated with that code is below
as required by the license agreement.

The original was dependent on an Anaconda library which allowed
programmatic setting of number of threads to use. We use
instead the use of the MKL_NUM_THREADS environment variable
to choose number of threads before calling.

--
Copyright (c) 2016, Adrian Haas and ETH Zurich
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this 
list of conditions and the following disclaimer. Redistributions in binary 
form must reproduce the above copyright notice, this list of conditions and the 
following disclaimer in the documentation and/or other materials provided 
with the distribution.
Neither the name of ETH Zurich nor the names of its contributors may be used 
to endorse or promote products derived from this software without specific 
prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
--
"""

import ctypes
import hashlib
import numpy
import warnings
import scipy.sparse
import sys
from scipy.sparse import SparseEfficiencyWarning
from time import time
from wrapper import SolverLibraryWrapper

class PardisoError(Exception):
    """Errors during PARDISO solves."""
    
    def __init__(self, value):
        """Initialise with error code."""

        self.value = value
    

    def __str__(self):
        """Return a string saying that it's a Pardiso error, followed by error code."""

        return ('The Pardiso solver failed with error code {}. '
                'See Pardiso documentation for details.'.format(self.value))
      

class PardisoWrapper(SolverLibraryWrapper):
    """
    Python interface to the Intel MKL PARDISO library for solving large sparse linear systems of equations Ax=b.
    
    Pardiso documentation: https://software.intel.com/en-us/node/470282
    """

    """Shared library names for different operating systems."""
    PLATFORM_MKL_LIBRARY = \
    {
        'linux2' : 'libmkl_rt.so',
        'win32' : 'mkl_rt.dll'
    }

    """Matrix type constant: real positive-definite """
    MATRIX_TYPE_REAL_POSITIVE_DEFINITE=2

    """Calculation: do analysis and factorise"""
    PHASE_ANALYSE_AND_FACTORISE = 12

    """Calculation: do solve"""
    PHASE_SOLVE = 33

    """Calculation: backward substitution to solve L.T x = b where A = L.LT"""
    PHASE_FORWARD_SUBSTITUTION = 331

    """Calculation: backward substitution to solve L.T x = b where A = L.LT"""
    PHASE_BACKWARD_SUBSTITUTION = 333

    """Phase to free all internal memory"""
    PHASE_FREE_ALL_MEMORY = -1


    """Set this iparm index to nonzero to make use of any other iparm."""
    IPARM_IN_USE = 0

    """Set this iparm index to one of ORDERING_MINIMUM_DEGREE, ORDERING_METIS, ORDERING_PARMETIS."""
    IPARM_ORDERING = 1

    """Set this iparm index to one of PERMUTATION_IGNORE, PERMUTATION_SET, PERMUTATION_GET."""
    IPARM_PERMUTATION = 4

    """Use zero-based ordering."""
    IPARM_ZEROBASED = 34
    
    ORDERING_MINIMUM_DEGREE = 0
    ORDERING_METIS = 2
    ORDERING_PARMETIS = 3

    """Means don't use permutation info (Pardiso will work internally with permutations)."""
    PERMUTATION_IGNORE = 0

    """Set permutation to one we already know."""
    PERMUTATION_SET = 1

    """Get the permutation produced by an analyse operation."""
    PERMUTATION_GET = 2

    """In core or out of core mode"""
    IPARM_CORE_MODE = 59
    
    IPARM_IN_CORE = 0     # In core mode (default)
    IPARM_AUTO_CORE = 1   # Pardiso decides whether to use in core or out of core
    IPARM_OUT_OF_CORE = 2 # Out of core

    @staticmethod
    def cholesky(A, mtype=MATRIX_TYPE_REAL_POSITIVE_DEFINITE, printstats=False):

        return PardisoWrapper(A, mtype, printstats)
    
    def __init__(self, A, mtype, printstats, permutation=None, out_of_core=False):

        # Check it's square
        if A.shape[0] != A.shape[1]:
            raise ValueError('Matrix A needs to be square, but has shape: {}'.format(A.shape))
        
        # Check of expected sparse format
        if not scipy.sparse.isspmatrix_csr(A):
            raise ValueError('Only CSR is supported')

        # scipy allows csr matrices with empty rows. a square matrix with an empty row is singular. calling 
        # pardiso with a matrix A that contains empty rows leads to a segfault, same applies for csc with 
        # empty columns
        if not numpy.diff(A.indptr).all():
            raise ValueError('Matrix A is singular, because it contains empty row(s)'.format(row_col))

        # check floating point is 64-bit
        if A.dtype != numpy.float64:
            raise TypeError('Pardiso wrapper currently only supports float64, but matrix A has dtype: {}'.format(A.dtype))

        # Always use just upper-triangular part of A
        # NOTE: this adds unnecessary expense but at time of making this
        #       the inputs make fully-populated symmetrix matrices which then crash Pardiso
        self.A = scipy.sparse.triu(A, format='csr')

        # store it
        self.mtype = mtype
        self.printstats = printstats

        # get pardiso library
        self.libmkl = ctypes.CDLL( PardisoWrapper.PLATFORM_MKL_LIBRARY[sys.platform] )
        self._mkl_pardiso = self.libmkl.pardiso
        
        # determine 32bit or 64bit architecture
        if ctypes.sizeof(ctypes.c_void_p) == 8:
            self._pt_type = (ctypes.c_int64, numpy.int64)
        else:
            self._pt_type = (ctypes.c_int32, numpy.int32)

        # declare arguments for native call to pardiso
        self._mkl_pardiso.argtypes = [ctypes.POINTER(self._pt_type[0]),    # pt
                                      ctypes.POINTER(ctypes.c_int32),      # maxfct
                                      ctypes.POINTER(ctypes.c_int32),      # mnum
                                      ctypes.POINTER(ctypes.c_int32),      # mtype
                                      ctypes.POINTER(ctypes.c_int32),      # phase
                                      ctypes.POINTER(ctypes.c_int32),      # n
                                      ctypes.POINTER(None),                # a
                                      ctypes.POINTER(ctypes.c_int32),      # ia
                                      ctypes.POINTER(ctypes.c_int32),      # ja
                                      ctypes.POINTER(ctypes.c_int32),      # perm
                                      ctypes.POINTER(ctypes.c_int32),      # nrhs
                                      ctypes.POINTER(ctypes.c_int32),      # iparm
                                      ctypes.POINTER(ctypes.c_int32),      # msglvl
                                      ctypes.POINTER(None),                # b
                                      ctypes.POINTER(None),                # x
                                      ctypes.POINTER(ctypes.c_int32)]      # error

        # handle used to refer to this instance of data on native side
        self.pt = numpy.zeros(64, dtype=self._pt_type[1])

        # permutation information (or None if none yet exists)
        self.permutation = permutation

        # run in out of core mode
        self.out_of_core = out_of_core

        # ensure indices sorted (NumPy allows this not to be the case)
        if not self.A.has_sorted_indices:
            self.A.sort_indices()

        # run factorisation
        b = numpy.zeros((self.A.shape[0],1))
        self._call_pardiso(b, phase=PardisoWrapper.PHASE_ANALYSE_AND_FACTORISE)
    
    def solve_A(self, b):
        """
        solve Ax=b for x
        
        --- Parameters ---
        b: numpy ndarray
           right-hand side(s), b.shape[0] needs to be the same as A.shape[0]
           
        --- Returns ---
        x: numpy ndarray
           solution of the system of linear equations, same shape as input b
        """


        # Check shape
        if b.shape[0] != self.A.shape[0]:
            raise ValueError("Dimension mismatch: Matrix A {} and array b {}".format(self.A.shape, b.shape))

        # Check it's a dense vector or matrix
        if scipy.sparse.isspmatrix(b):
            warnings.warn('Pardiso requires the right-hand side b to be a dense array for maximum efficiency', 
                          SparseEfficiencyWarning)
            b = b.todense()
        
        # Pardiso expects fortran (column-major) order if b is a matrix
        if b.ndim == 2:
            b = numpy.asfortranarray(b)
                    
        # Type conversion
        if b.dtype != numpy.float64:
            raise TypeError('Dtype {} for array b is not supported'.format(str(b.dtype)))

        # Solve it
        x = self._call_pardiso(b, phase=PardisoWrapper.PHASE_SOLVE)

        return SolverLibraryWrapper.return_like(b, x)

    def diag_inv_A(self, max_blocksize = 1000):
        """
        compute diagonal elements of A^(-1), where
        solve AA^(-1)=I
           
        --- Returns ---
        diag(A^(-1)): numpy ndarray
           the diagonal elements of the inverse of A
        """    
        print(self.A.shape[0])
	identy_matrix = scipy.sparse.eye(self.A.shape[0], format = 'csc' )
    
	if max_blocksize is None:
	    # Compute full diagonal of covariance matrix at once. 
	    L_inverse = self.solve_forward_substitution(identy_matrix)
	    inverse_diagonal = numpy.multiply(L_inverse,L_inverse).sum(0)
	    
	else:
	    # Compute diagonal elements for at most max_blocksize elements at a time. 
	    number_of_state_parameters = identy_matrix.shape[0]
           
	    inverse_diagonal = numpy.zeros(number_of_state_parameters)   
	    counter = 0
	    t0 =time()
	    while counter < number_of_state_parameters:
		index_1 = counter
		if counter+max_blocksize < number_of_state_parameters:
		    index_2 = counter+max_blocksize
		else:
		    index_2 = number_of_state_parameters
            
		sub_identy_matrix = identy_matrix[index_1:index_2,:].T
    
		sub_L_inverse = self.solve_forward_substitution(sub_identy_matrix)
    
		#if index_2-index_1 == 1:
		 #   x = scipy.sparse.csc_matrix(x[:,np.newaxis])
                
		#sub_P = np.asarray(x.multiply(x).sum(0)).ravel()

		inverse_diagonal[index_1:index_2] = numpy.multiply(sub_L_inverse,sub_L_inverse).sum(0).ravel()
		del sub_L_inverse
		counter+=max_blocksize
	    print('Diagonal elements computation done in '+str(time()-t0)+' seconds')
	return inverse_diagonal

    def solve_backward_substitution(self, b):
        """
        solve LTx=b for x where A = L.LT
        """

        # Check shape
        if b.shape[0] != self.A.shape[0]:
            raise ValueError("Dimension mismatch: Matrix A {} and array b {}".format(self.A.shape, b.shape))

        # Check it's a dense vector or matrix
        if scipy.sparse.isspmatrix(b):
            warnings.warn('Pardiso requires the right-hand side b to be a dense array for maximum efficiency', 
                          SparseEfficiencyWarning)
            b = b.todense()
        
        # Pardiso expects fortran (column-major) order if b is a matrix
        if b.ndim == 2:
            b = numpy.asfortranarray(b)
                    
        # Type conversion
        if b.dtype != numpy.float64:
            raise TypeError('Dtype {} for array b is not supported'.format(str(b.dtype)))

        # Solve it
        x = self._call_pardiso(b, phase=PardisoWrapper.PHASE_BACKWARD_SUBSTITUTION)

        return SolverLibraryWrapper.return_like(b, x)

    def solve_forward_substitution(self, b):
        """
        solve Lx=b for x where A = L.LT
        """

        # Check shape
        if b.shape[0] != self.A.shape[0]:
            raise ValueError("Dimension mismatch: Matrix A {} and array b {}".format(self.A.shape, b.shape))

        # Check it's a dense vector or matrix
        if scipy.sparse.isspmatrix(b):
            warnings.warn('Pardiso requires the right-hand side b to be a dense array for maximum efficiency', 
                          SparseEfficiencyWarning)
            b = b.todense()
        
        # Pardiso expects fortran (column-major) order if b is a matrix
        if b.ndim == 2:
            b = numpy.asfortranarray(b)
                    
        # Type conversion
        if b.dtype != numpy.float64:
            raise TypeError('Dtype {} for array b is not supported'.format(str(b.dtype)))

        # Solve it
        x = self._call_pardiso(b, phase=PardisoWrapper.PHASE_FORWARD_SUBSTITUTION)

        return SolverLibraryWrapper.return_like(b, x)


    def _call_pardiso(self, b, phase):
        """Internal call to Pardiso shared library."""
        
        x = numpy.zeros_like(b)
        pardiso_error = ctypes.c_int32(0)
        c_int32_p = ctypes.POINTER(ctypes.c_int32)
        c_float64_p = ctypes.POINTER(ctypes.c_double)

        # Buffers
        ia = self.A.indptr
        ja = self.A.indices

        # Parameters - these have meanings as defined here:
        # https://software.intel.com/en-us/mkl-developer-reference-fortran-pardiso-iparm-parameter
        iparm = numpy.zeros(64, dtype=numpy.int32)
        iparm[PardisoWrapper.IPARM_IN_USE] = 1
        iparm[PardisoWrapper.IPARM_ZEROBASED] = 1

        # use par-metis
        iparm[PardisoWrapper.IPARM_ORDERING] = PardisoWrapper.ORDERING_PARMETIS

        # run out of core
        if self.out_of_core:
            iparm[PardisoWrapper.IPARM_CORE_MODE] = PardisoWrapper.IPARM_OUT_OF_CORE
        
        # populate perm
        if self.permutation is None:

            iparm[PardisoWrapper.IPARM_PERMUTATION] = PardisoWrapper.PERMUTATION_GET
            self.permutation = numpy.zeros((self.A.shape[0],), dtype=numpy.int32)

        else:

            iparm[PardisoWrapper.IPARM_PERMUTATION] = PardisoWrapper.PERMUTATION_SET
        
        self._mkl_pardiso(self.pt.ctypes.data_as(ctypes.POINTER(self._pt_type[0])), # pt
                          ctypes.byref(ctypes.c_int32(1)), # maxfct
                          ctypes.byref(ctypes.c_int32(1)), # mnum
                          ctypes.byref(ctypes.c_int32(self.mtype)), # mtype -> 2 for real positive-definite
                          ctypes.byref(ctypes.c_int32(phase)), # phase -> 13 
                          ctypes.byref(ctypes.c_int32(self.A.shape[0])), #N -> number of equations/size of matrix
                          self.A.data.ctypes.data_as(c_float64_p), # A -> non-zero entries in matrix
                          ia.ctypes.data_as(c_int32_p), # ia -> csr-indptr
                          ja.ctypes.data_as(c_int32_p), # ja -> csr-indices
                          self.permutation.ctypes.data_as(c_int32_p), # get or set permutation as required
                          ctypes.byref(ctypes.c_int32(1 if b.ndim == 1 else b.shape[1])), # nrhs
                          iparm.ctypes.data_as(c_int32_p), # iparm-array
                          ctypes.byref(ctypes.c_int32(self.printstats)), # msg-level -> 1: statistical info is printed
                          b.ctypes.data_as(c_float64_p), # b -> right-hand side vector/matrix
                          x.ctypes.data_as(c_float64_p), # x -> output
                          ctypes.byref(pardiso_error)) # pardiso error
        
        if pardiso_error.value != 0:
            raise PardisoError(pardiso_error.value)
        else:
            return numpy.ascontiguousarray(x) # change memory-layout back from fortran to c order
               
    def free_memory(self):
        """Request pardiso to free all internal memory."""

        self.A = scipy.sparse.csr_matrix((0,0))
        b = numpy.zeros(0)
        phase = PardisoWrapper.PHASE_FREE_ALL_MEMORY
        self._call_pardiso(b, phase)
