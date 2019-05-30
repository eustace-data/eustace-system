"""Tests for grid module."""
# pylint: disable=missing-docstring,bad-whitespace

__version__ = "$Revision: 395 $"
__author__ = "Joel R. Mitchelson"

import unittest
import numpy
from ..grid import GridAxis
from ..grid import GridLatLon

class TestGridAxis(unittest.TestCase):

    def test_compute_indices(self):

        # test data
        testvalues = numpy.array([0.25, 0.65, 0.8, 1.6, 1.65], numpy.float32)

        # regrid onto grid boxes starting at -0.5, 0.0, 0.5, 1.0, 1.5
        axis = GridAxis(start=-0.5, increment=0.5, length=5)
        result = axis.compute_indices(testvalues)

        # check
        self.assertEqual([1, 2, 2, 4, 4], result.tolist())

    def test_compute_indices_error_check(self):

        # check for underflow (this would regrid onto negative indices)
        self.assertRaises(ValueError, GridAxis.compute_indices,
                          GridAxis(start=-1.5, increment=0.5, length=1000),
                          numpy.array([-3.0, -1.0, 2.1, 5.9], numpy.float32))

        # check for overflow (this would regrid outside the axis length)
        self.assertRaises(ValueError, GridAxis.compute_indices,
                          GridAxis(start=1.0, increment=1.0, length=3),
                          numpy.array([2.1, 5.9], numpy.float32))

    def test_compute_indices_wrap(self):

        testvalues = numpy.array([0.0, -3.2, 4.0, 2.2, 3.99], numpy.float32)
        # without wrap this should have an error due to 4.0 being on the
        # boundary
        self.assertRaises(ValueError, GridAxis.compute_indices,
                          GridAxis(start=-4.0, increment=1.0, length=8),
                          testvalues)
        # but with wrap we should get a result with 4.0 mapped to grid zero
        result = GridAxis(start=-4.0, increment=1.0, length=8,
                          wrap=True).compute_indices(testvalues)
        numpy.testing.assert_equal(result, [4, 0, 0, 6, 7])


class TestGridLatLon(unittest.TestCase):

    def test_init(self):
        result = GridLatLon(
            axis_lat=GridAxis(start=-90.0, increment=2.0, length=88),
            axis_lon=GridAxis(start=-180.0, increment=1.0, length=360))
        self.assertEqual(-90.0, result.axis_lat.start)
        self.assertEqual(2.0, result.axis_lat.increment)
        self.assertEqual(88, result.axis_lat.length)
        self.assertEqual(-180.0, result.axis_lon.start)
        self.assertEqual(1.0, result.axis_lon.increment)
        self.assertEqual(360, result.axis_lon.length)

    def test_compute_indices(self):

        # make a tiny test grid
        grid = GridLatLon(
            axis_lat=GridAxis(start=-2.0, increment=2.0, length=2),
            axis_lon=GridAxis(start=-2.5, increment=1.0, length=5))

        # make test data, including edge cases
        # grid box (row,col):     (0,4) (0,3) (1,4) (1,4) (1,0)
        # grid flat index   :       4     3     9     9     5
        coord_lat = numpy.array([-0.2, -1.2,  1.8,  1.7,  0.0], numpy.float32)  # noqa (allow multiple spaces)
        coord_lon = numpy.array([ 1.6,  1.4,  2.3,  1.5, -2.5], numpy.float32)  # noqa (allow multiple spaces)

        # compute
        result = grid.compute_indices(coord_lat, coord_lon)

        # check
        self.assertEqual((5,), result.shape)
        self.assertEqual(numpy.int32, result.dtype)
        self.assertEqual([4, 3, 9, 9, 5], result.tolist())
