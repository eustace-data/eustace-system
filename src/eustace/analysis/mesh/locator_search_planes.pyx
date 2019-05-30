# cython
"""Check of 3D points against array of normals."""

import numpy
cimport numpy

ctypedef numpy.int64_t indextype
ctypedef numpy.float64_t valuetype

def locator_search_planes_exhaustive(
    numpy.ndarray[indextype, ndim=1] indices,
    numpy.ndarray[valuetype, ndim=2] normal0,
    numpy.ndarray[valuetype, ndim=2] normal1,
    numpy.ndarray[valuetype, ndim=2] normal2,
    numpy.ndarray[valuetype, ndim=2] points):
    """Perform simple exhaustive search and return result in indices array."""

    cdef int m = normal0.shape[0] # m normals
    cdef int n = points.shape[0]  # n points to evaluate (should equal number of indices to fill)

    # Loop over points
    for j in range(n):

    	# Get 3D vector for one point only
        p = points[j,:]

	# Loop over all normals and check occupancy
        for i in range(m):
            n0 = normal0[i,:]
            n1 = normal1[i,:]
            n2 = normal2[i,:]
            result0 = p[0]*n0[0] + p[1]*n0[1] + p[2]*n0[2]
            result1 = p[0]*n1[0] + p[1]*n1[1] + p[2]*n1[2]
            result2 = p[0]*n2[0] + p[1]*n2[1] + p[2]*n2[2]
            if (result0 <= 0.0) and (result1 <= 0.0) and (result2 <= 0.0):
                indices[j] = i
                break
