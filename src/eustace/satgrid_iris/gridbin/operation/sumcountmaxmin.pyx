# cython: profile=True
"""Define simultaneous sum, count, max, min."""

import numpy
cimport numpy
from operation import GridBinOperation

ctypedef numpy.int32_t indextype
ctypedef numpy.int32_t counttype
ctypedef numpy.float32_t valuetype

class GridBinOperationSumCountMaxMin(GridBinOperation):
    """Sum inputs onto bins and keep a count of total per bin."""

    OUTPUTTOTAL = 'total'
    OUTPUTCOUNT = 'count'
    OUTPUTMAX = 'maximum'
    OUTPUTMIN = 'minimum'
    INPUT = 'inputvalue'

    def __init__(self):
        super(GridBinOperationSumCountMaxMin, self).__init__(
            outputids=[GridBinOperationSumCountMaxMin.OUTPUTTOTAL, 
                       GridBinOperationSumCountMaxMin.OUTPUTCOUNT,
                       GridBinOperationSumCountMaxMin.OUTPUTMAX,
                       GridBinOperationSumCountMaxMin.OUTPUTMIN],
            inputids=[GridBinOperationSumCountMaxMin.INPUT])

    def operate(self, 
        numpy.ndarray[indextype, ndim=1] mapfrom,
        numpy.ndarray[indextype, ndim=1] mapto,
        numpy.ndarray[valuetype, ndim=1] total,
        numpy.ndarray[counttype, ndim=1] count,
        numpy.ndarray[valuetype, ndim=1] maximum,
        numpy.ndarray[valuetype, ndim=1] minimum,
        numpy.ndarray[valuetype, ndim=1] inputvalue):
        """Perform binning operation to compute total and count of input values."""

        cdef int n = mapto.shape[0]
        cdef indextype indexto
        cdef indextype indexfrom
        for i in range(n):
            indexto = mapto[i]
            indexfrom = mapfrom[i]
            value = inputvalue[indexfrom]
            total[indexto] += value
            count[indexto] += 1
            maximum[indexto] = max(maximum[indexto], value)
            minimum[indexto] = min(minimum[indexto], value)
