"""Location within a mesh."""

import numpy
from eustace.analysis.mesh.locator_search_planes import locator_search_planes_exhaustive
from geometry import triangle_area_spherical
from geometry import triangle_area_euclidean
from geometry import dividerows
from geometry import barycentric_coordinates_planar
from geometry import barycentric_coordinates_spherical
from mesh import MeshIcosahedron

class TriangleLocator(object):
    """Find triangles on a mesh."""

    def __init__(self, mesh):
        """Construct with reference to given mesh."""

        # store mesh reference
        self.mesh = mesh

        # outward-facing normals to planes bounding the triangle
        self.normal0 = numpy.cross(mesh.points[mesh.triangles[:,1],:], mesh.points[mesh.triangles[:,2],:], axis=1)
        self.normal1 = numpy.cross(mesh.points[mesh.triangles[:,2],:], mesh.points[mesh.triangles[:,3],:], axis=1)
        self.normal2 = numpy.cross(mesh.points[mesh.triangles[:,3],:], mesh.points[mesh.triangles[:,1],:], axis=1)

        # normalise
        self.normal0 = dividerows(self.normal0, numpy.linalg.norm(self.normal0, axis=1))
        self.normal1 = dividerows(self.normal1, numpy.linalg.norm(self.normal1, axis=1))
        self.normal2 = dividerows(self.normal2, numpy.linalg.norm(self.normal2, axis=1))

    def locate_vector_indices_exhaustive(self, v):
        """Find indices of triangles that the given points are on.
           Note this must be shifted if the triangle number is required.
           v should have 3 columns and one row per point to evaluate."""
        
        # check v as expected
        if len(v.shape) != 2:
            raise ValueError('Expected matrix with 1 row per point')
        if v.shape[1] != 3:
            raise ValueError('Expected 3D coordinate per point')

        # allocate array for indices (init to invalid value of -1)
        indices = -numpy.ones((v.shape[0],), numpy.int64)

        # call Cython method to compute values of indices
        locator_search_planes_exhaustive(indices, self.normal0, self.normal1, self.normal2, v)

        # result
        return indices

    def barycentric_coordinates_at_indices(self, indices, v, evaluation_method=barycentric_coordinates_planar):
        """Compute barycentric coordinates for given vector and corresponding indices.
           By default this uses computes barycentres in the plane defined by the triangle points.
           To use spherical triangular geometry set evaluation_method=barycentric_coordinates_spherical."""

        # check v and indices as expected
        if len(v.shape) != 2:
            raise ValueError('Expected matrix with 1 row per point')
        if v.shape[1] != 3:
            raise ValueError('Expected 3D coordinate per point')
        if len(indices.shape) != 1:
            raise ValueError('Vector of indices expected to be one-dimensional')
        if v.shape[0] != indices.shape[0]:
            raise ValueError('Expected one index value per point')

        # get point values
        point0 = self.mesh.points[self.mesh.triangles[indices,1],:]
        point1 = self.mesh.points[self.mesh.triangles[indices,2],:]
        point2 = self.mesh.points[self.mesh.triangles[indices,3],:]

        # retrieve result
        return evaluation_method(point0, point1, point2, v)


class TriangleLocatorHierarchical(TriangleLocator):
    """Hierarchical method for triangle location assuming subdivided icosahedron structure."""

    def __init__(self, toplevelmesh):
        """Initialise using the top level icosahedron."""
        
        super(TriangleLocatorHierarchical, self).__init__(toplevelmesh)

        
    def locate_vector_hierarchical(self, v, levels=29, evaluation_method=barycentric_coordinates_planar):
        """Get triangle number of given points after subdivision, and corresponding barycentres.
           By default this uses computes barycentres in the plane defined by the triangle points.
           To use spherical triangular geometry set evaluation_method=barycentric_coordinates_spherical."""

        # first use top level mesh to get the super-triangle numbers
        top_level_indices = self.locate_vector_indices_exhaustive(v)

        # obtain super-triangle points surrounding each input vector
        point0 = self.mesh.points[self.mesh.triangles[top_level_indices,1],:]
        point1 = self.mesh.points[self.mesh.triangles[top_level_indices,2],:]
        point2 = self.mesh.points[self.mesh.triangles[top_level_indices,3],:]

        # initialise triangle numbers
        trianglenumbers = numpy.uint64(top_level_indices) << MeshIcosahedron.BASE_HASH_SHIFT

        # evaluate barycentres (planar geometry)
        barycentres = evaluation_method(point0, point1, point2, v)

        # and now perform repeated subdivision to get finest possible triangle number
        for i in range(1, int(levels) + 1):

            # find the subdivision number for each barycentre
            subnumbers = TriangleLocatorHierarchical.triangle_barycentre_subdivision_number(barycentres)

            # update points accordingly
            TriangleLocatorHierarchical.triangle_points_to_subtriangle(point0, point1, point2, subnumbers)

            # re-evaluate barycentres with new points
            barycentres = evaluation_method(point0, point1, point2, v)

            # update triangle numbers
            trianglenumbers = trianglenumbers | (subnumbers << (MeshIcosahedron.BASE_HASH_SHIFT -  numpy.uint64(2*i)))
    
        return trianglenumbers, barycentres

    @staticmethod
    def triangle_barycentre_subdivision_number(barycentres):
        """Subdivide the triangle into 4. Return indices 0...3 of the triangles
           in which given barrycentres are found."""

        # Subdivision looks like this:
        #
        #                   p1
        #                   /\
        #                  /  \
        #                 /    \
        #                /  t1  \
        #               /        \
        #         p2'  x----l1----x  p0'
        #             / \         /\
        #            /   \  t3   /  \
        #           /    l0    l2    \
        #          /  t0   \   /  t2  \
        #         /         \ /        \
        #    p0  x-----------x----------x  p2
        #                   p1'
        #
        # The lines l0, l1, l2 delineate the corner triangles t0, t1, t2.
        # Barrycentric coordinates (x0, x1, x2) on these lines satisfy:
        #
        #   l0:   x0 - x1 - x2 = 0
        #   l1:   x1 - x2 - x0 = 0
        #   l2:   x2 - x0 - x1 = 0
        #
        # And testing >= 0 tells us if we are in one of those corner triangles.
        #

        # evaluate inclusion for each sub-triangle, checking that points on 
        # boundaries are not assigned to multiple sub-triangles.

        inclusiontest = numpy.zeros((barycentres.shape[0],4), dtype=numpy.bool)
        
        inclusiontest[:,0] = (barycentres[:,0] - barycentres[:,1] - barycentres[:,2]) >= 0.0
        notfoundindicator = numpy.logical_not(numpy.any(inclusiontest, axis = 1))
        
        inclusiontest[:,1] = numpy.logical_and(notfoundindicator,(barycentres[:,1] - barycentres[:,2] - barycentres[:,0]) >= 0.0)
        notfoundindicator = numpy.logical_not(numpy.any(inclusiontest, axis = 1))
        
        inclusiontest[:,2] = numpy.logical_and(notfoundindicator,(barycentres[:,2] - barycentres[:,0] - barycentres[:,1]) >= 0.0)
        notfoundindicator = numpy.logical_not(numpy.any(inclusiontest, axis = 1))
        
        inclusiontest[:,3] = notfoundindicator
        
        # one way to get the triangle numbers
        return numpy.uint64( numpy.mod(numpy.nonzero(inclusiontest.ravel())[0], 4) )

    @staticmethod
    def triangle_points_to_subtriangle(point0, point1, point2, subnumbers):
        """Using subdivision number to convert points to corresponding sub-triangle."""

        # select indices for each possible triangle number
        indices_t0 = numpy.nonzero(subnumbers == 0)
        indices_t1 = numpy.nonzero(subnumbers == 1)
        indices_t2 = numpy.nonzero(subnumbers == 2)
        indices_t3 = numpy.nonzero(subnumbers == 3)

        # t0
        point1[indices_t0,:] = 0.5*(point0[indices_t0,:] + point1[indices_t0,:])
        point2[indices_t0,:] = 0.5*(point2[indices_t0,:] + point0[indices_t0,:])

        # t1
        point0[indices_t1,:] = 0.5*(point0[indices_t1,:] + point1[indices_t1,:])
        point2[indices_t1,:] = 0.5*(point1[indices_t1,:] + point2[indices_t1,:])

        # t2
        point0[indices_t2,:] = 0.5*(point2[indices_t2,:] + point0[indices_t2,:])
        point1[indices_t2,:] = 0.5*(point1[indices_t2,:] + point2[indices_t2,:])

        # t3
        # (note require intermediate variables here as point arrays used multiple times)
        update0 = 0.5*(point1[indices_t3,:] + point2[indices_t3,:])
        update1 = 0.5*(point2[indices_t3,:] + point0[indices_t3,:])
        update2 = 0.5*(point0[indices_t3,:] + point1[indices_t3,:])
        point0[indices_t3,:] = update0
        point1[indices_t3,:] = update1
        point2[indices_t3,:] = update2

        # normalise
        point0[:,:] = dividerows(point0, numpy.linalg.norm(point0, axis=1))
        point1[:,:] = dividerows(point1, numpy.linalg.norm(point1, axis=1))
        point2[:,:] = dividerows(point2, numpy.linalg.norm(point2, axis=1))

    @staticmethod
    def triangle_barycentre_subdivide(barycentres, subnumbers):

        # To transform barrycentric coordinates from the whole triangle
        # to coordinates to subtriangle t0, consider a coordinate system with p0 as
        # origin, and note p2' = 0.5 p1, p1' = 0.5 p2.  So coefficients
        # of p0, p2', p1' are (x0', x1', x2') where:
        #
        # [ x0' ] = [ x0 -  x1 - x2  ]
        # [ x1' ]   [      2.x1      ]
        # [ x2' ]   [      2.x2      ]
        #
        # And similarly for t1, t2.
        #
        # Similar process for t3 gives:
        #
        # [ x0' ] = [ -x0 + x1 + x2 ]
        # [ x1' ] = [  x0 - x1 + x2 ]
        # [ x2' ] = [  x0 + x1 - x2 ]

        # select indices for each possible triangle number
        indices_t0 = numpy.nonzero(subnumbers == 0)
        indices_t1 = numpy.nonzero(subnumbers == 1)
        indices_t2 = numpy.nonzero(subnumbers == 2)
        indices_t3 = numpy.nonzero(subnumbers == 3)

        # build new coordinates
        newcoords = numpy.zeros(shape=barycentres.shape)

        # Triangle t0
        newcoords[indices_t0,0] = barycentres[indices_t0,0] - barycentres[indices_t0,1] - barycentres[indices_t0,2]
        newcoords[indices_t0,1] = barycentres[indices_t0,1]*2.0
        newcoords[indices_t0,2] = barycentres[indices_t0,2]*2.0

        # Triangle t1
        newcoords[indices_t1,1] = barycentres[indices_t1,1] - barycentres[indices_t1,2] - barycentres[indices_t1,0]
        newcoords[indices_t1,2] = barycentres[indices_t1,2]*2.0
        newcoords[indices_t1,0] = barycentres[indices_t1,0]*2.0

        # Triangle t2
        newcoords[indices_t2,2] = barycentres[indices_t2,2] - barycentres[indices_t2,0] - barycentres[indices_t2,1]
        newcoords[indices_t2,0] = barycentres[indices_t2,0]*2.0
        newcoords[indices_t2,1] = barycentres[indices_t2,1]*2.0

        # Triangle t3
        newcoords[indices_t3,0] = -barycentres[indices_t3,0] + barycentres[indices_t3,1] + barycentres[indices_t3,2]
        newcoords[indices_t3,1] =  barycentres[indices_t3,0] - barycentres[indices_t3,1] + barycentres[indices_t3,2]
        newcoords[indices_t3,2] =  barycentres[indices_t3,0] + barycentres[indices_t3,1] - barycentres[indices_t3,2]

        # return result
        return newcoords
