# cython: profile=True
"""Occupancy count operation."""

import numpy
cimport numpy
from operation import GridBinOperation

ctypedef numpy.int32_t indextype
ctypedef numpy.int32_t counttype

class GridBinOperationCount(GridBinOperation):
    """Sum inputs onto bins and keep a count of total per bin."""

    OUTPUT = 'count'
    INPUT = 'inputvalue'

    def __init__(self):
        super(GridBinOperationCount, self).__init__(
            outputids=[GridBinOperationCount.OUTPUT],
            inputids=[GridBinOperationCount.INPUT])

    def operate(self, numpy.ndarray[indextype, ndim=1] mapfrom, numpy.ndarray[indextype, ndim=1] mapto, numpy.ndarray[counttype, ndim=1] count, numpy.ndarray inputvalue):
        """Perform binning operation to compute total and count of input values."""

        cdef int n = mapto.shape[0]
        cdef indextype indexto
        cdef indextype numberone = 1
        for i in range(n):
            indexto = mapto[i]
            count[indexto] += numberone
