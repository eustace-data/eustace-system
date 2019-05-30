# cython: profile=True
"""Binning operation of sum of squared deviation from some mean value which already exists on the grid."""

import numpy
cimport numpy
from operation import GridBinOperation

ctypedef numpy.int32_t indextype
ctypedef numpy.float32_t valuetype

class GridBinOperationSumSqDev(GridBinOperation):
    """Sum inputs onto bins and keep a count of total per bin."""

    OUTPUT = 'sumsqdev'
    EXISTINGMEAN = 'existingmean'
    INPUTVALUE = 'inputvalue'

    def __init__(self):
        super(GridBinOperationSumSqDev, self).__init__(
            outputids=[GridBinOperationSumSqDev.OUTPUT, GridBinOperationSumSqDev.EXISTINGMEAN],
            inputids=[GridBinOperationSumSqDev.INPUTVALUE])

    def operate(self, numpy.ndarray[indextype, ndim=1] mapfrom, numpy.ndarray[indextype, ndim=1] mapto, numpy.ndarray[valuetype, ndim=1] sumsqdev, numpy.ndarray[valuetype, ndim=1] existingmean, numpy.ndarray[valuetype, ndim=1] inputvalue):
        """Perform binning operation to compute total and count of input values."""

        cdef int n = mapto.shape[0]
        cdef indextype indexto
        cdef indextype indexfrom
        cdef valuetype dev
        for i in range(n):
            indexto = mapto[i]
            indexfrom = mapfrom[i]
            dev = inputvalue[indexfrom] - existingmean[indexto]
            sumsqdev[indexto] += dev * dev
