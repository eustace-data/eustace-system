"""Test geometric methods"""

from ..geometry import triangle_circumcentre
from ..geometry import norm2
from ..geometry import scalerows
from ..geometry import dividerows
from ..geometry import triangle_area_spherical
from ..geometry import triangle_area_euclidean
from ..geometry import normalised_cross_product
from ..geometry import barycentric_coordinates_planar
from ..geometry import cartesian_to_polar2d
from ..geometry import polar2d_to_cartesian
from ..geometry import vector_angle

import unittest
import numpy

class TestGeometry(unittest.TestCase):


    def test_triangle_circumcentre(self):

        # a triangle as 3 lots of 3 coordinates
        triangles = numpy.array([ [ 0.0, 0.0, 1.0,
                                    4.0, 0.0, 1.0,
                                    0.0, 3.0, 1.0 ],

                                  [ 0.0, 0.0,-1.0,
                                   -4.0, 0.0,-1.0,
                                    0.0, 3.0,-1.0 ] ])

        # call the method
        centres = triangle_circumcentre(triangles)
        
        # compute radius vectors
        radius_vectors = numpy.reshape(triangles - numpy.tile(centres, (1,3)), (6, 3))

        # get radii
        radii = numpy.linalg.norm(radius_vectors, axis=1)

        # check all the same as expected (2.5)
        self.assertEqual((6,), radii.shape)
        numpy.testing.assert_almost_equal(radii, 2.5)

        
    def test_norm2(self):

        a = numpy.array ( [ [ 3, 4, 5 ],
                            [ 0.5, 1.2, 2.2 ] ])
        result = norm2(a)
        numpy.testing.assert_almost_equal(result, [ 50, 6.53 ])

    def test_scalerows(self):

        a = numpy.array( [ [ -2.0, 6.0, 8.0 ],
                           [ 3.33333333, -1.0, 1000.0 ] ] )
        s = numpy.array( [ 7.0 , 3.0 ] )
        result = scalerows(s, a)
        numpy.testing.assert_almost_equal(result, [ [ -14.0, 42.0, 56.0 ], [ 10.0, -3.0, 3000.0 ] ])

    def test_dividerows(self):

        a = numpy.array( [ [ -56.0, 42.0, 3.5 ],
                           [  22.0,  1.0,-888.8 ] ] )
        s = numpy.array( [ 7.0, 2.0 ] )
        result = dividerows(a, s)
        numpy.testing.assert_almost_equal(result, [ [ -8.0, 6.0, 0.5 ], [ 11.0, 0.5, -444.4 ] ])

    def test_normalised_cross_product(self):

        a = numpy.array( [ [ -9.0, 3.0, 2.0 ],
                           [  4.0, 0.0, 5.0 ] ] )

        b = numpy.array( [ [  8.0, 6.0, 1.5 ],
                           [  1.2, 2.2, 3.2 ] ] )

        result = normalised_cross_product(a, b)

        numpy.testing.assert_almost_equal(result,
            [ [-0.08957499,  0.35232829, -0.93157989],
              [-0.70322362, -0.43472006,  0.5625789 ] ])
              

    def test_triangle_area_spherical(self):

        # Create test data
        p0 = numpy.array( [ [ numpy.cos(numpy.radians(5))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(5))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(10)) ],

                            [ numpy.cos(numpy.radians(15))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(15))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(10)) ] ] )


        p1 = numpy.array( [ [ numpy.cos(numpy.radians(35))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(35))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(10)) ],

                            [ numpy.cos(numpy.radians(45))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(45))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(10)) ] ] )

        p2 = numpy.array( [ [ numpy.cos(numpy.radians(20))*numpy.cos(numpy.radians(-35)),
                              numpy.sin(numpy.radians(20))*numpy.cos(numpy.radians(-35)),
                              numpy.sin(numpy.radians(-35)) ],

                            [ numpy.cos(numpy.radians(30))*numpy.cos(numpy.radians(-35)),
                              numpy.sin(numpy.radians(30))*numpy.cos(numpy.radians(-35)),
                              numpy.sin(numpy.radians(-35)) ] ] )

        # Check these are set up correctly with unit radius
        numpy.testing.assert_almost_equal(numpy.linalg.norm(p0, axis=1), [ 1.0, 1.0 ])
        numpy.testing.assert_almost_equal(numpy.linalg.norm(p1, axis=1), [ 1.0, 1.0 ])
        numpy.testing.assert_almost_equal(numpy.linalg.norm(p2, axis=1), [ 1.0, 1.0 ])

        # Outward-facing normals as unit vectors
        n0 = numpy.cross(p0, p1, axis=1)
        n1 = numpy.cross(p1, p2, axis=1)
        n2 = numpy.cross(p2, p0, axis=1)
        n0 = dividerows(n0, numpy.linalg.norm(n0, axis=1))
        n1 = dividerows(n1, numpy.linalg.norm(n1, axis=1))
        n2 = dividerows(n2, numpy.linalg.norm(n2, axis=1))

        # Check these also unit vectors
        numpy.testing.assert_almost_equal(numpy.linalg.norm(n0, axis=1), [ 1.0, 1.0 ])
        numpy.testing.assert_almost_equal(numpy.linalg.norm(n1, axis=1), [ 1.0, 1.0 ])
        numpy.testing.assert_almost_equal(numpy.linalg.norm(n2, axis=1), [ 1.0, 1.0 ])

        # Look at angles in degrees
        self.assertAlmostEqual(75.8649314014, numpy.degrees(numpy.arccos(-numpy.sum(n0[0,:] * n1[0,:]))))
        self.assertAlmostEqual(40.6644678599, numpy.degrees(numpy.arccos(-numpy.sum(n1[0,:] * n2[0,:]))))
        self.assertAlmostEqual(75.8649314014, numpy.degrees(numpy.arccos(-numpy.sum(n2[0,:] * n0[0,:]))))

        # Second set should match first because they are just rotated around z-axis by 10 degrees
        self.assertAlmostEqual(75.8649314014, numpy.degrees(numpy.arccos(-numpy.sum(n0[1,:] * n1[1,:]))))
        self.assertAlmostEqual(40.6644678599, numpy.degrees(numpy.arccos(-numpy.sum(n1[1,:] * n2[1,:]))))
        self.assertAlmostEqual(75.8649314014, numpy.degrees(numpy.arccos(-numpy.sum(n2[1,:] * n0[1,:]))))

        # Excess in degrees =  75.8649314014 + 40.6644678599 + 75.8649314014 - 180 = 12.39433066
        # Excess in radians = 0.21632188
        result = triangle_area_spherical(n0, n1, n2)
        numpy.testing.assert_almost_equal(result, [ 0.21632188, 0.21632188 ])
        
    def test_triangle_area_euclidean(self):

        # Create test data
        p0 = numpy.array( [ [ numpy.cos(numpy.radians(5))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(5))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(10)) ],

                            [ 0, 0, 0 ] ] )


        p1 = numpy.array( [ [ numpy.cos(numpy.radians(35))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(35))*numpy.cos(numpy.radians(10)),
                              numpy.sin(numpy.radians(10)) ],

                            [ 4, 0, 0 ] ] )

        p2 = numpy.array( [ [ numpy.cos(numpy.radians(20))*numpy.cos(numpy.radians(-35)),
                              numpy.sin(numpy.radians(20))*numpy.cos(numpy.radians(-35)),
                              numpy.sin(numpy.radians(-35)) ],

                            [ 4, 3, 0 ] ] )

        # vectors for first row so we can check values
        a = p1[0,:] - p0[0,:]
        b = p2[0,:] - p0[0,:]
        numpy.testing.assert_almost_equal(a, [-0.17435298, 0.47903087, 0 ])
        numpy.testing.assert_almost_equal(b, [-0.21130913, 0.19433485, -0.74722461])

        # do the computation
        result = triangle_area_euclidean(p0, p1, p2)

        # area for first triangle should be  (1/2) | a x b |
        #                   
        #   = 0.5 * sqrt( (0.47903087*-0.74722461)**2 + (-0.17435298*-0.74722461)**2 + (-0.17435298*0.19433485 - 0.47903087*-0.21130913)**2 )
        #   = 0.19341118
        #
        # and 2nd triangle we should know the area is 6
        #
        # so:
        numpy.testing.assert_almost_equal(result, [ 0.19341118, 6.00000000 ])


    def test_barycentric_coordinates_planar(self):

        point0 = numpy.array([ [ 1.0, 1.0, 1.0 ],
                               [ 2.0, 2.0, 2.0 ] ])
        point1 = numpy.array([ [ 3.0, 5.0, 1.0 ],
                               [ 6.0, 2.0, 7.0 ] ])
        point2 = numpy.array([ [ 6.0, 1.0, 1.0 ],
                               [10.0, 2.0, 2.0 ] ])

        testvalues = numpy.vstack( ( 0.2 * point0[0, :] +
                                     0.7 * point1[0, :] +
                                     0.1 * point2[0, :],

                                     0.4 * point0[1, :] +
                                     0.3 * point1[1, :] +
                                     0.3 * point2[1, :]
                                     ) )

        result = barycentric_coordinates_planar(point0, point1, point2, testvalues)

        numpy.testing.assert_almost_equal( [ [ 0.2, 0.7, 0.1 ], [ 0.4, 0.3, 0.3 ] ], result)    
        
    def test_polar2d_convention(self):

        # Check convention
        # x-axis is 0E 0N, y-axis is 90E 0N, and z is 0E 90N
        basis = cartesian_to_polar2d( numpy.array( [ [ 1.0, 0.0, 0.0 ],
                                                     [ 0.0, 1.0, 0.0 ],
                                                     [ 0.0, 0.0, 1.0 ] ] ) )
        numpy.testing.assert_almost_equal(basis, [ [  0.0,  0.0 ],
                                                   [  0.0, 90.0 ],
                                                   [ 90.0,  0.0 ] ])
        
    def test_cartesian_to_polar2d_roundtrip(self):

        # some unit vectors made using
        # v = numpy.random.randn(5,3) ; v / numpy.tile(numpy.linalg.norm(v, axis=1).T, (3,1)).T
        testdata = numpy.array([[-0.31416367, -0.16220734,  0.93540899],
                                [ 0.22179823,  0.11004287, -0.96886331],
                                [ 0.18196861,  0.54028381, -0.82157217],
                                [-0.48406711, -0.87112947, -0.08253771],
                                [ 0.02503293, -0.08542615,  0.99602998]])

        # do round trip and check result
        polar2d = cartesian_to_polar2d(testdata)
        cartesian = polar2d_to_cartesian(polar2d)
        numpy.testing.assert_almost_equal(cartesian, testdata)
        
        # also check it works if we wrongly scale the unit vectors
        alternative_polar2d = cartesian_to_polar2d(123.0*testdata)
        numpy.testing.assert_almost_equal(alternative_polar2d, polar2d)

    def test_vector_angle(self):

        # Something obviously at 30 degrees
        numpy.testing.assert_almost_equal(vector_angle(numpy.array([ 1.0, 0.0, 0.0 ]), numpy.array([ numpy.cos(numpy.radians(30)), numpy.sin(numpy.radians(30)), 0.0 ])), 30.0)

        numpy.testing.assert_almost_equal(
            vector_angle(
                numpy.array([ [ 1.0, 0.0, 0.0 ], [ 3.0, 0.0, 0.0 ] ]),
                numpy.array([ [ numpy.cos(numpy.radians(30)), numpy.sin(numpy.radians(30)), 0.0 ], [ 20.0*numpy.cos(numpy.radians(-173.0)), 0.0, 20.0*numpy.sin(numpy.radians(-173.0)) ] ]) ),
            [ 30.0, 173.0 ])
        
