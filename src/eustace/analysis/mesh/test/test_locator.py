"""Test mesh locator tools."""

from ..locator import TriangleLocator
from ..locator import TriangleLocatorHierarchical
from ..mesh import MeshIcosahedron
from ..mesh import MeshIcosahedronSubdivision
from ..geometry import barycentric_coordinates_planar
from ..geometry import barycentric_coordinates_spherical
import unittest
import numpy

class TestTriangleLocator(unittest.TestCase):

    def test_locate_vector_indices_exhaustive(self):
        
        mesh = MeshIcosahedronSubdivision.build(0)

        locator = TriangleLocator(mesh)

        surfacepoints = numpy.array( [ [ 1.0, 0.0, 0.0 ],
                                       [ 0.0, 1.0, 0.0 ],
                                       [ 0.0, 0.0, 1.0 ] ] )

        result = locator.locate_vector_indices_exhaustive(surfacepoints)

        numpy.testing.assert_equal([ 4, 8, 0 ], result)


    def test_barycentric_coordinates_at_indices(self):

        mesh0 = MeshIcosahedron.build()
        mesh1 = MeshIcosahedronSubdivision.subdivide(mesh0)
        mesh2 = MeshIcosahedronSubdivision.subdivide(mesh1)
        mesh3 = MeshIcosahedronSubdivision.subdivide(mesh2)

        locator = TriangleLocator(mesh0)

        indices1 = locator.locate_vector_indices_exhaustive(mesh1.points)
        points1 = locator.barycentric_coordinates_at_indices(indices1, mesh1.points, barycentric_coordinates_planar)

        indices2 = locator.locate_vector_indices_exhaustive(mesh2.points)
        points2 = locator.barycentric_coordinates_at_indices(indices2, mesh2.points)


class TestTriangleLocatorHierarchical(unittest.TestCase):

    def test_locator_hierarchical(self):

        # invent some surface points
        surfacepoints = numpy.array( [ [ numpy.cos(numpy.radians(5))*numpy.cos(numpy.radians(10)),
                                         numpy.sin(numpy.radians(5))*numpy.cos(numpy.radians(10)),
                                         numpy.sin(numpy.radians(10)) ],
                                       
                                       [ numpy.cos(numpy.radians(15))*numpy.cos(numpy.radians(10)),
                                         numpy.sin(numpy.radians(15))*numpy.cos(numpy.radians(10)),
                                         numpy.sin(numpy.radians(10)) ] ] )


        # reference against top-level icosahedral mesh
        toplevelmesh = MeshIcosahedron.build()        
        locator = TriangleLocatorHierarchical(toplevelmesh)

        # use hierarchical method to compute triangle numbers
        result_trianglenumbers, result_barycentres = locator.locate_vector_hierarchical(surfacepoints, levels=5)

        # compare with a level 5 mesh direct exhaustive computation
        mesh5 = MeshIcosahedronSubdivision.build(5)
        indices5 = TriangleLocator(mesh5).locate_vector_indices_exhaustive(surfacepoints)
        numbers5 = mesh5.triangles[indices5,0]

        # check triangle numbers
        numpy.testing.assert_equal(numbers5, result_trianglenumbers)
        numpy.testing.assert_equal(indices5, result_trianglenumbers >> numpy.uint64(MeshIcosahedron.BASE_HASH_SHIFT - 2*5))

        # check barycentres correspond
        expected_barycentres = barycentric_coordinates_planar(
                mesh5.points[mesh5.triangles[indices5,1],:],
                mesh5.points[mesh5.triangles[indices5,2],:],
                mesh5.points[mesh5.triangles[indices5,3],:],
                surfacepoints)
        numpy.testing.assert_almost_equal(expected_barycentres, result_barycentres)
        

    def test_triangle_barycentre_subdivision_number(self):

        testdata = numpy.array( [ [ 0.7 , 0.2 , 0.1  ],
                                  [ 0.32, 0.32, 0.36 ],
                                  [ 0.1 , 0.8 , 0.1  ],
                                  [ 0.05, 0.05, 0.9  ] ] )
        
        indices = TriangleLocatorHierarchical.triangle_barycentre_subdivision_number(testdata)

        numpy.testing.assert_equal([ 0, 3, 1, 2 ], indices)

    def test_triangle_barycentre_subdivision_number_at_vertex(self):

        testdata = numpy.array( [ [ 0.5 , 0.5 , 0.0  ],
                                  [ 0.0 , 0.5 , 0.5  ],
                                  [ 0.5 , 0.0 , 0.5  ] ])
        indices = TriangleLocatorHierarchical.triangle_barycentre_subdivision_number(testdata)
        numpy.testing.assert_equal([ 0, 1, 0 ], indices)

    def test_triangle_points_to_subtriangle(self):

        point0 = numpy.array( [ [ 1.0, 1.0, 1.0 ],
                                [ 2.0, 2.0, 2.0 ] ] )
                              
        point1 = numpy.array( [ [ 4.0, 7.0, 1.0 ],
                                [ 8.0, 15.0, 2.0 ] ] )

        point2 = numpy.array( [ [ 7.0, 1.0, 1.0 ],
                                [16.0, 2.0, 2.0 ] ] )

        TriangleLocatorHierarchical.triangle_points_to_subtriangle(point0, point1, point2, numpy.array( [ 1, 3 ] ))

        test00 = numpy.array([ 2.5, 4.0, 1.0 ]) / numpy.sqrt(2.5*2.5 + 4.0*4.0 + 1.0*1.0)
        test01 = numpy.array([12.0, 8.5, 2.0 ]) / numpy.sqrt(12.0*12.0 + 8.5*8.5 + 2.0*2.0)

        test10 = numpy.array([ 4.0, 7.0, 1.0 ]) / numpy.sqrt(4.0*4.0 + 7.0*7.0 + 1.0*1.0)
        test11 = numpy.array([ 9.0, 2.0, 2.0 ]) / numpy.sqrt(9.0*9.0 + 2.0*2.0 + 2.0*2.0)

        test20 = numpy.array([ 5.5, 4.0, 1.0 ]) / numpy.sqrt(5.5*5.5 + 4.0*4.0 + 1.0*1.0)
        test21 = numpy.array([ 5.0, 8.5, 2.0 ]) / numpy.sqrt(5.0*5.0 + 8.5*8.5 + 2.0*2.0)

        numpy.testing.assert_almost_equal(numpy.vstack((test00, test01)), point0)
        numpy.testing.assert_almost_equal(numpy.vstack((test10, test11)), point1)
        numpy.testing.assert_almost_equal(numpy.vstack((test20, test21)), point2)



    def test_triangle_barycentre_subdivide(self):

        testdata = numpy.array( [ [ 0.7 , 0.2 , 0.1  ],
                                  [ 0.32, 0.32, 0.36 ],
                                  [ 0.1 , 0.8 , 0.1  ],
                                  [ 0.05, 0.05, 0.9  ] ] )
        
        indices = TriangleLocatorHierarchical.triangle_barycentre_subdivision_number(testdata)

        barycentres = TriangleLocatorHierarchical.triangle_barycentre_subdivide(testdata, indices)
        
        # Invent a basis for checking purposes
        P = numpy.array( [ [ -1.0, 0.0        ],
                           [  0.0, 1.73205081 ],
                           [  1.0, 0.0        ] ] )

        # Coordinates with respect to this basis
        bigtrianglecoords = numpy.matmul(testdata, P)

        # Now make some sub-triangle bases

        Psub0 = numpy.array([ [ -1.0, 0.0        ],
                              [ -0.5, 0.8660254  ],
                              [  0.0, 0.0        ] ])

        Psub1 = numpy.array([ [ -0.5, 0.8660254  ],
                              [  0.0, 1.73205081 ],
                              [  0.5, 0.8660254  ] ])

        Psub2 = numpy.array([ [  0.0, 0.0        ],
                              [  0.5, 0.8660254  ],
                              [  1.0, 0.0        ] ])

        Psub3 = numpy.array([ [  0.5, 0.8660254  ],
                              [  0.0, 0.0        ],
                              [ -0.5, 0.8660254  ] ])

        # And check expression w.r.t. sub-bases using sub-coordinates is correct
        numpy.testing.assert_almost_equal(bigtrianglecoords[0,:], numpy.matmul(barycentres[0,:], Psub0))
        numpy.testing.assert_almost_equal(bigtrianglecoords[1,:], numpy.matmul(barycentres[1,:], Psub3))
        numpy.testing.assert_almost_equal(bigtrianglecoords[2,:], numpy.matmul(barycentres[2,:], Psub1))
        numpy.testing.assert_almost_equal(bigtrianglecoords[3,:], numpy.matmul(barycentres[3,:], Psub2))


    def test_face_centre_at_poles(self):

        poles =  numpy.array( [ [ 0.0, 0.0, 1.0 ],
                                [ 0.0, 0.0,-1.0 ] ] )        
        mesh = MeshIcosahedronSubdivision.build(0, face_centre_at_poles=True)
        result_trianglenumbers, result_barycentres = TriangleLocatorHierarchical(mesh).locate_vector_hierarchical(poles, levels=0)
        third = 1.0 / 3.0
        numpy.testing.assert_almost_equal(
            result_barycentres,
            [ [ third, third, third ],
              [ third, third, third ] ])
        
