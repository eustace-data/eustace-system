import unittest
from ..mesh import MeshIcosahedron
from ..mesh import MeshIcosahedronSubdivision
import numpy

class TestMeshIcosahedron(unittest.TestCase):

    def test_build(self):

        m = MeshIcosahedron.build()

        # Check data structure sizes
        self.assertEqual((12,3), m.points.shape)
        self.assertEqual((20,4), m.triangles.shape)
        self.assertEqual((30,6), m.edges.shape)

        # Types
        self.assertEqual(numpy.float64, m.points.dtype)
        self.assertEqual(numpy.uint64, m.triangles.dtype)
        self.assertEqual(numpy.uint64, m.edges.dtype)
        self.assertEqual(numpy.uint64, m.level.dtype)             

        # Must be level 0
        self.assertEqual(numpy.uint64(0), m.level)

        # all points at unit radius
        numpy.testing.assert_almost_equal(numpy.linalg.norm(m.points, axis=1), 1.0)

        # Check on clockwiseness: normal to clockwise vectors around triangle
        # should have direction the same as vector from origin through triangle centre
        v1 = m.points[m.triangles[:,2],:] - m.points[m.triangles[:,1],:]
        v2 = m.points[m.triangles[:,3],:] - m.points[m.triangles[:,2],:]
        v3 = m.points[m.triangles[:,1],:] - m.points[m.triangles[:,3],:]
        centre =  m.points[m.triangles[:,1],:] + m.points[m.triangles[:,2],:] + m.points[m.triangles[:,3],:]
        normal = numpy.cross(v2, v1)
        centre /= numpy.transpose(numpy.tile(numpy.linalg.norm(centre, axis=1), (3,1)))
        normal /= numpy.transpose(numpy.tile(numpy.linalg.norm(normal, axis=1), (3,1)))
        numpy.testing.assert_almost_equal(centre, normal)

        # also all distances should be the same
        r = (1.0 + numpy.sqrt(5.0)) * 0.5
        sidelength = 2.0 / numpy.sqrt(1.0 + r*r)
        self.assertAlmostEqual(1.05146222424, sidelength)
        numpy.testing.assert_equal(numpy.linalg.norm(v1, axis=1), sidelength)
        numpy.testing.assert_equal(numpy.linalg.norm(v2, axis=1), sidelength)
        numpy.testing.assert_equal(numpy.linalg.norm(v3, axis=1), sidelength)

        # and each point should appear 5 times in triangle list
        allpoints = numpy.sort(m.triangles[:,1:].ravel())
        expectedpoints = numpy.transpose(numpy.tile(range(12), (5, 1))).ravel()
        numpy.testing.assert_equal(allpoints, expectedpoints)

        # check first edge
        self.assertEqual(3, m.edges[0,0]) # start point in clockwise order for first triangle
        self.assertEqual(0, m.edges[0,1]) # end point in clockwise order for first triangle
        self.assertEqual(5, int(m.edges[0,2]) >> 58) # first triangle number
        self.assertEqual(2, m.edges[0,3]) # index in clockwise ordering for first triangle
        self.assertEqual(4, int(m.edges[0,4]) >> 58) # second triangle number
        self.assertEqual(0, m.edges[0,5]) # index in clockwise ordering for second triangle

        # check last edge
        self.assertEqual(10, m.edges[29,0]) # start point in clockwise order for first triangle
        self.assertEqual( 9, m.edges[29,1]) # end point in clockwise order for first triangle
        self.assertEqual(10, int(m.edges[29,2]) >> 58) # first triangle number
        self.assertEqual( 2, m.edges[29,3]) # index in clockwise ordering for first triangle
        self.assertEqual(11, int(m.edges[29,4]) >> 58) # second triangle number
        self.assertEqual( 0, m.edges[29,5]) # index in clockwise ordering for second triangle

        # get the three edges for every triangle (will have duplicates because triangles share common edges)
        edges_with_duplicates = numpy.vstack((m.edges[:,0:4], numpy.hstack((m.edges[:,0:2], m.edges[:,4:6]))))

        # sort according to triangle number
        sort_by_triangle = numpy.argsort(edges_with_duplicates[:,2] | edges_with_duplicates[:,3])
        sorted_edges_with_duplicates = edges_with_duplicates[sort_by_triangle,:]

        # should get three per triangle
        numpy.testing.assert_equal(sorted_edges_with_duplicates[0::3,3],0)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[1::3,3],1)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[2::3,3],2)        


class TestMeshIcosahedronSubdivision(unittest.TestCase):

    def test_subdivide_once(self):

        parent = MeshIcosahedron.build()
        m = MeshIcosahedronSubdivision.subdivide(parent, project_onto_sphere=False)

        # Types
        self.assertEqual(numpy.float64, m.points.dtype)
        self.assertEqual(numpy.uint64, m.triangles.dtype)
        self.assertEqual(numpy.uint64, m.edges.dtype)
        self.assertEqual(numpy.uint64, m.level.dtype)             

        # Check data structure sizes
        self.assertEqual(( 42,3), m.points.shape)
        self.assertEqual(( 80,4), m.triangles.shape)
        self.assertEqual((120,6), m.edges.shape)
             
        # Must be level 1
        self.assertEqual(numpy.uint64(1), m.level)

        # First split of first edge
        self.assertEqual( 3, m.edges[0,0]) # start point in clockwise order for first triangle
        self.assertEqual(12, m.edges[0,1]) # new mid point
        self.assertEqual((4*5+2) << 56, m.edges[0,2]) # first triangle number
        self.assertEqual( 2, m.edges[0,3]) # index in clockwise ordering for first triangle

        # get the three edges for every triangle (will have duplicates because triangles share common edges)
        edges_with_duplicates = numpy.vstack((m.edges[:,0:4], numpy.hstack((m.edges[:,0:2], m.edges[:,4:6]))))

        # sort according to triangle number
        sort_by_triangle = numpy.argsort(edges_with_duplicates[:,2] | edges_with_duplicates[:,3])
        sorted_edges_with_duplicates = edges_with_duplicates[sort_by_triangle,:]
 
        # should get three per triangle
        numpy.testing.assert_equal(sorted_edges_with_duplicates[0::3,3],0)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[1::3,3],1)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[2::3,3],2)

        # Check on clockwiseness: normal to clockwise vectors around triangle
        # should have positive dot-product with line from origin through triangle centre
        v1 = m.points[m.triangles[:,2],:] - m.points[m.triangles[:,1],:]
        v2 = m.points[m.triangles[:,3],:] - m.points[m.triangles[:,2],:]
        v3 = m.points[m.triangles[:,1],:] - m.points[m.triangles[:,3],:]
        centre =  m.points[m.triangles[:,1],:] + m.points[m.triangles[:,2],:] + m.points[m.triangles[:,3],:]
        normal = numpy.cross(v2, v1)
        centre /= numpy.transpose(numpy.tile(numpy.linalg.norm(centre, axis=1), (3,1)))
        normal /= numpy.transpose(numpy.tile(numpy.linalg.norm(normal, axis=1), (3,1)))
        dotprod = numpy.sum(centre*normal,1)
        # numpy.testing.assert_almost_equal(dotprod, 1.0, decimal=3)

        # first 20 points should appear 5 times in triangle list
        # and last 22 points 6 times
        allpoints = numpy.sort(m.triangles[:,1:].ravel())
        expectedbasepoints = numpy.transpose(numpy.tile(range(0,12), (5, 1))).ravel()
        expectedsubpoints = numpy.transpose(numpy.tile(range(12,42), (6, 1))).ravel()
        numpy.testing.assert_equal(allpoints[0:60], expectedbasepoints)
        numpy.testing.assert_equal(allpoints[60:], expectedsubpoints)

        # side lengths should be identical
        r = (1.0 + numpy.sqrt(5.0)) * 0.5
        sidelength = (2.0 / numpy.sqrt(1.0 + r*r)) / 2.0
        numpy.testing.assert_equal(numpy.linalg.norm(v1, axis=1), sidelength)
        numpy.testing.assert_equal(numpy.linalg.norm(v2, axis=1), sidelength)
        numpy.testing.assert_equal(numpy.linalg.norm(v3, axis=1), sidelength)

        # projecting onto sphere should give same but normalised points
        msphere = MeshIcosahedronSubdivision.subdivide(parent, project_onto_sphere=True)
        numpy.testing.assert_equal(m.triangles, msphere.triangles)
        numpy.testing.assert_equal(m.edges, msphere.edges)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(msphere.points, axis=1), 1.0)
        numpy.testing.assert_almost_equal(numpy.cross(m.points, msphere.points, axis=1), numpy.zeros((42,3)))


    def test_subdivide_twice(self):

        # Run subdivision twice
        m2 = MeshIcosahedronSubdivision.build(2, project_onto_sphere=False)

        # Types
        self.assertEqual(numpy.float64, m2.points.dtype)
        self.assertEqual(numpy.uint64, m2.triangles.dtype)
        self.assertEqual(numpy.uint64, m2.edges.dtype)

        # Check data structure sizes
        self.assertEqual((162,3), m2.points.shape)
        self.assertEqual((320,4), m2.triangles.shape)
        self.assertEqual((480,6), m2.edges.shape)

        # get the three edges for every triangle (will have duplicates because triangles share common edges)
        edges_with_duplicates = numpy.vstack((m2.edges[:,0:4], numpy.hstack((m2.edges[:,0:2], m2.edges[:,4:6]))))

        # sort according to triangle number
        sort_by_triangle = numpy.argsort(edges_with_duplicates[:,2] | edges_with_duplicates[:,3])
        sorted_edges_with_duplicates = edges_with_duplicates[sort_by_triangle,:]
 
        # should get three per triangle
        numpy.testing.assert_equal(sorted_edges_with_duplicates[0::3,3],0)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[1::3,3],1)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[2::3,3],2)

        # first 20 points should appear 5 times in triangle list
        # and last 142 points 6 times
        allpoints = numpy.sort(m2.triangles[:,1:].ravel())
        expectedbasepoints = numpy.transpose(numpy.tile(range(0,12), (5, 1))).ravel()
        expectedsubpoints = numpy.transpose(numpy.tile(range(12,162), (6, 1))).ravel()
        numpy.testing.assert_equal(allpoints[0:60], expectedbasepoints)
        numpy.testing.assert_equal(allpoints[60:], expectedsubpoints)

        # side lengths should be identical
        v1 = m2.points[m2.triangles[:,2],:] - m2.points[m2.triangles[:,1],:]
        v2 = m2.points[m2.triangles[:,3],:] - m2.points[m2.triangles[:,2],:]
        v3 = m2.points[m2.triangles[:,1],:] - m2.points[m2.triangles[:,3],:]
        r = (1.0 + numpy.sqrt(5.0)) * 0.5
        sidelength = (2.0 / numpy.sqrt(1.0 + r*r)) / numpy.power(2.0, 2.0)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(v1, axis=1), sidelength)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(v2, axis=1), sidelength)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(v3, axis=1), sidelength)


    def test_subdivide_x7(self):

        # do subdivision seven times
        m = MeshIcosahedronSubdivision.build(number_of_subdivisions=7, project_onto_sphere=False)

        # Compute expected numbers for outcome of subdivision
        p = 12
        t = 20
        e = 30
        for subdivision in range(7):
            p = p + e
            e = 3*t + 2*e
            t = 4*t

        # Should satisfy euler characteristic
        self.assertEqual(2, p - e + t)

        # Total number of triangles = 20*4^7 = 327680
        self.assertEqual(327680, t)
        self.assertEqual(163842, p)
        self.assertEqual(491520, e)

        # Types
        self.assertEqual(numpy.float64, m.points.dtype)
        self.assertEqual(numpy.uint64, m.triangles.dtype)
        self.assertEqual(numpy.uint64, m.edges.dtype)
        self.assertEqual(numpy.uint64, m.level.dtype)

        # Ranges
        self.assertEqual((p,3), m.points.shape)
        self.assertEqual((t,4), m.triangles.shape)
        self.assertEqual((e,6), m.edges.shape)

        # get the three edges for every triangle (will have duplicates because triangles share common edges)
        edges_with_duplicates = numpy.vstack((m.edges[:,0:4], numpy.hstack((m.edges[:,0:2], m.edges[:,4:6]))))

        # sort according to triangle number
        sort_by_triangle = numpy.argsort(edges_with_duplicates[:,2] | edges_with_duplicates[:,3])
        sorted_edges_with_duplicates = edges_with_duplicates[sort_by_triangle,:]
 
        # should get three per triangle
        numpy.testing.assert_equal(sorted_edges_with_duplicates[0::3,3],0)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[1::3,3],1)
        numpy.testing.assert_equal(sorted_edges_with_duplicates[2::3,3],2)

        # Also if right-shift the triangle numbers should get contiguous integers
        # (the shift amount is: BASE_HASH_SHIFT - 2*subdivisions = 44 in this case)
        uniqueids = sorted_edges_with_duplicates[:,2] >> numpy.uint64(44)
        numpy.testing.assert_array_equal(uniqueids[0::3], range(327680))
        numpy.testing.assert_array_equal(uniqueids[1::3], range(327680))
        numpy.testing.assert_array_equal(uniqueids[2::3], range(327680))

        # first 20 points should appear 5 times in triangle list
        # the rest should appear 6 times
        allpoints = numpy.sort(m.triangles[:,1:].ravel())
        expectedbasepoints = numpy.transpose(numpy.tile(range(0,12), (5, 1))).ravel()
        expectedsubpoints = numpy.transpose(numpy.tile(range(12, p), (6, 1))).ravel()
        numpy.testing.assert_equal(allpoints[0:60], expectedbasepoints)
        numpy.testing.assert_equal(allpoints[60:], expectedsubpoints)

        # side lengths should be identical
        v1 = m.points[m.triangles[:,2],:] - m.points[m.triangles[:,1],:]
        v2 = m.points[m.triangles[:,3],:] - m.points[m.triangles[:,2],:]
        v3 = m.points[m.triangles[:,1],:] - m.points[m.triangles[:,3],:]
        r = (1.0 + numpy.sqrt(5.0)) * 0.5
        sidelength = (2.0 / numpy.sqrt(1.0 + r*r)) / numpy.power(2.0, 7.0)
        self.assertAlmostEqual(0.0082145486, sidelength)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(v1, axis=1), sidelength)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(v2, axis=1), sidelength)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(v3, axis=1), sidelength)


    def test_subdivide_x3_radius(self):
       
        m = MeshIcosahedronSubdivision.build(number_of_subdivisions=3, project_onto_sphere=True)
        numpy.testing.assert_almost_equal(numpy.linalg.norm(m.points, axis=1), numpy.ones(642))
