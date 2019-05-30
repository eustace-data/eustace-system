# cython: profile=True
"""Binning operation of sum of squared value."""

import numpy
cimport numpy
import cython
import sys
from operation import GridBinOperation

ctypedef numpy.int32_t indextype
ctypedef numpy.float32_t valuetype

class GridBinOperationSumSq(GridBinOperation):
    """Compute sum of squares on each bin."""

    OUTPUT = 'sumsq'
    INPUT = 'inputvalue'

    def __init__(self):
        super(GridBinOperationSumSq, self).__init__(
            outputids=[GridBinOperationSumSq.OUTPUT],
            inputids=[GridBinOperationSumSq.INPUT])

    @cython.boundscheck(False)
    def operate(self, numpy.ndarray[indextype, ndim=1] mapfrom, numpy.ndarray[indextype, ndim=1] mapto, numpy.ndarray[valuetype, ndim=1] sumsq, numpy.ndarray[valuetype, ndim=1] inputvalue):
        """Perform binning operation to compute sum of squares."""

        cdef int n = mapto.shape[0]
        cdef indextype indexto
        cdef indextype indexfrom
        cdef valuetype value
        for i in range(n):
            indexto = mapto[i]
            indexfrom = mapfrom[i]
            value = inputvalue[indexfrom]
            sumsq[indexto] += value * value

class GridBinOperationSumConstSq(GridBinOperation):
    """Compute sum of squares on each bin where input is a constant value."""

    OUTPUT = 'sumsq'
    INPUT = 'singlevalue'

    def __init__(self):
        super(GridBinOperationSumConstSq, self).__init__(
            outputids=[GridBinOperationSumConstSq.OUTPUT],
            inputids=[GridBinOperationSumConstSq.INPUT])

    @cython.boundscheck(False)
    def operate(self, numpy.ndarray[indextype, ndim=1] mapfrom, numpy.ndarray[indextype, ndim=1] mapto, numpy.ndarray[valuetype, ndim=1] sumsq, numpy.ndarray[valuetype, ndim=1] singlevalue):
        """Perform binning operation to compute sum of squares."""

        cdef valuetype sq = singlevalue[0] * singlevalue[0]
        # sys.stderr.write('gridbin.operation.sumsq.GridBinOperationSumConstSq: {0}\n'.format(sq))
        cdef int n = mapto.shape[0]
        cdef indextype indexto
        for i in range(n):
            indexto = mapto[i]
            sumsq[indexto] += sq
