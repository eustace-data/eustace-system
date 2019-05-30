"""Geometric methods on a mesh."""

import numpy

def triangle_circumcentre(points):
    """Compute circumcentres for array of triangles (assumed rows are triangles and columns are 9 coordinates).
       This is given by:

                  ( |a|^2 b  - |b|^2  a ) x ( a x b )
       p = c  +   -----------------------------------
                            2 | a x b |^2

       Where c is vertex 0
             a is vertex 1 - vertex 0
             b is vertex 2 - vertex 0
    """

    c = points[:, 0:3]
    a = points[:, 3:6] - c
    b = points[:, 6:9] - c    
    a_cross_b = numpy.cross(a, b)
    numerator = numpy.cross(scalerows(norm2(a), b) - scalerows(norm2(b), a), a_cross_b)
    denominator = 2 * norm2(a_cross_b)
    return c + dividerows(numerator, denominator)

def norm2(v):
    """Square of 2-norm of points."""

    return numpy.sum(numpy.multiply(v, v), axis=1)

def scalerows(s, v):
    """Multiply row-wise."""

    return numpy.multiply(numpy.transpose(numpy.tile(s, (3,1))), v)

def dividerows(v, s):
    """Divide row-wise."""

    return numpy.divide(v, numpy.transpose(numpy.tile(s, (3,1))))

def normalised_cross_product(a, b):
    """Compute (a x b) / | a x b | for one or more rows of 3D vectors."""

    c = numpy.cross(a, b, axis=1)
    c = dividerows(c, numpy.linalg.norm(c, axis=1))
    return c

def triangle_area_euclidean(point0, point1, point2):
    """Area of triangles using usual euclidean geometry.
       Points can have multiple rows for calculation of multiple triangles."""

    # Vector for first side
    side0 = point1 - point0

    # Vector for second side
    side1 = point2 - point0

    # Area using cross-product formula
    area = 0.5 * numpy.linalg.norm(numpy.cross(side0, side1, axis=1), axis=1)
    return area

def triangle_area_spherical(n0, n1, n2):
    """Area of triangles on unit sphere with sides defined by the great circles
       whose normals are n0, n1, n2.  Normals can have multiple rows
       for calculation of multiple triangles."""

    # compute dot product for inside angles (requires negative of outside dot product)
    dot_prod = numpy.empty(shape=n0.shape, dtype=n0.dtype)
    dot_prod[:,0] = -numpy.sum(numpy.multiply(n0, n1), axis=1)
    dot_prod[:,1] = -numpy.sum(numpy.multiply(n1, n2), axis=1)
    dot_prod[:,2] = -numpy.sum(numpy.multiply(n2, n0), axis=1)
    
    # compute angles
    angles = numpy.arccos(dot_prod)

    # area on unit sphere is equal to the "angular excess" (sum of angles minus pi)
    area = numpy.sum(angles, axis=1) - numpy.pi

    # result
    return area
        
def barycentric_coordinates_planar(point0, point1, point2, v):
    """Barycentric coordinates in plane defined by points."""

    # rescale along length so it lies on plane
    # n = (p1 - p0) x (p2 - p0)
    # want k s.t. (k v - p0 ).n = 0
    # => k  = (p0.n) / (v.n)

    n = numpy.cross(point1 - point0, point2 - point0, axis=1)
    k = numpy.sum(point0 * n, axis=1) / numpy.sum(v * n, axis=1)
    vp = scalerows(k, v)

    # Compute areas
    areas = numpy.transpose(numpy.vstack( (triangle_area_euclidean(point1, point2, vp), 
                                           triangle_area_euclidean(point2, point0, vp),
                                           triangle_area_euclidean(point0, point1, vp) ) ) )

    # Normalise to get result
    barycentres = dividerows(areas, numpy.sum(areas, axis=1))

    return barycentres

def barycentric_coordinates_spherical(point0, point1, point2, v):
    """Barycentric coordinates using areas on unit sphere through points."""

    # Consider the triangles:
    # Triangle A: point1 --- point2 --- v
    # Triangle B: point2 --- point0 --- v
    # Triangle C: point0 --- point1 --- v

    normalA0 = normalised_cross_product(point1, point2)
    normalA1 = normalised_cross_product(point2, v)
    normalA2 = normalised_cross_product(v, point1)

    normalB0 = normalised_cross_product(point2, point0)
    normalB1 = normalised_cross_product(point0, v)
    normalB2 = normalised_cross_product(v, point2)

    normalC0 = normalised_cross_product(point0, point1)
    normalC1 = normalised_cross_product(point1, v)
    normalC2 = normalised_cross_product(v, point0)

    # Compute areas on unit sphere
    areas = numpy.transpose(numpy.vstack( (triangle_area_spherical(normalA0, normalA1, normalA2), 
                                           triangle_area_spherical(normalB0, normalB1, normalB2),
                                           triangle_area_spherical(normalC0, normalC1, normalC2) ) ) )

    barycentres = dividerows(areas, numpy.sum(areas, axis=1))

    return barycentres

def cartesian_to_polar2d(v):
    """Convert unit vector v to polar coordinates using convention that
       x-axis is 0E 0N, y-axis is 90E 0N, and z is 0E 90N ."""
       
    radius = numpy.linalg.norm(v, axis=1)
    latitude = numpy.arcsin(v[:,2] / radius)
    longitude = numpy.arctan2(v[:,1], v[:,0])
    return numpy.degrees(numpy.vstack( (latitude, longitude) ).T)

def polar2d_to_cartesian(polar2d):
    """Convert polar coordinates polar2d to cartesian unit vector."""

    z = numpy.sin(numpy.radians(polar2d[:,0]))
    minor_r = numpy.cos(numpy.radians(polar2d[:,0]))
    x = numpy.cos(numpy.radians(polar2d[:,1])) * minor_r
    y = numpy.sin(numpy.radians(polar2d[:,1])) * minor_r
    return numpy.vstack( (x, y, z) ).T

def vector_angle(a, b):
    """Angle in degrees between two arrays of vectors (rows are vectors, columns are XYZ)."""

    last_axis = a.ndim - 1
    dotprod = (a * b).sum(axis=last_axis)
    mod_a2 = (a * a).sum(axis=last_axis)
    mod_b2 = (b * b).sum(axis=last_axis)
    return numpy.degrees(numpy.arccos(dotprod / numpy.sqrt(mod_a2 * mod_b2)))

