# cython
"""Aggregation operation for diagnostic grids."""

import numpy
cimport numpy

ctypedef numpy.int32_t indextype
ctypedef numpy.int64_t counttype
ctypedef numpy.float64_t valuetype

def diagnosticgrid_aggregate(
    numpy.ndarray[indextype, ndim=1] input_indices, 
    numpy.ndarray[valuetype, ndim=1] input_observation,
    numpy.ndarray[valuetype, ndim=1] input_uncertainty,
    numpy.ndarray[counttype, ndim=1] count,
    numpy.ndarray[valuetype, ndim=1] observation_min,
    numpy.ndarray[valuetype, ndim=1] observation_max,
    numpy.ndarray[valuetype, ndim=1] uncertainty_min,
    numpy.ndarray[valuetype, ndim=1] uncertainty_max,
    numpy.ndarray[valuetype, ndim=1] observation_weightedsum,
    numpy.ndarray[valuetype, ndim=1] precision_sum):
    """Aggregation operation."""

    cdef int n = input_indices.shape[0]
    cdef indextype gridindex
    cdef indextype numberone = 1
    for i in range(n):
        gridindex = input_indices[i]
        obs = input_observation[i]
        unc = input_uncertainty[i]
        count[gridindex] += numberone
        if obs < observation_min[gridindex]:
            observation_min[gridindex] = obs
        if obs > observation_max[gridindex]:
            observation_max[gridindex] = obs
        if obs < uncertainty_min[gridindex]:
            uncertainty_min[gridindex] = unc
        if unc > uncertainty_max[gridindex]:
            uncertainty_max[gridindex] = unc
        observation_weightedsum[gridindex] += obs / unc
        precision_sum[gridindex] += 1.0 / unc
