"""Tests for aggregate_grid module."""
# pylint: disable=missing-docstring, invalid-name, bad-whitespace

import unittest
import tempfile
import numpy
from ..grid import GridAxis
from ..grid import GridLatLon
from ..aggregate_grid import AggregateGridCount
from ..aggregate_grid import AggregateGridSum
from ..aggregate_grid import AggregateGridSumDevSq
from ..aggregate_grid import AggregateGridArray
from ..aggregate_grid import AggregateGridCountSumMinMax

__version__ = "$Revision: 472 $"
__author__ = "Joel R. Mitchelson"


class TestAggregateGridCount(unittest.TestCase):

    def test_aggregate_at_indices(self):

        # practice grid of 3 rows x 5 columns (starts and increments not
        # important)
        grid = GridLatLon(
            axis_lat=GridAxis(start=-20000.0, increment=100.0, length=3),
            axis_lon=GridAxis(start=250000.0, increment=1000.0, length=5))

        # aggregation object
        a = AggregateGridCount(grid.dimensions)

        # the indices to accumulate
        indices = numpy.array([0, 0, 6, 14, 14, 14], numpy.int32)

        # process
        a.aggregate_at_indices(indices)

        # check result
        numpy.testing.assert_equal(
            actual=a.count,
            desired=numpy.array(
                [[2, 0, 0, 0, 0],
                 [0, 1, 0, 0, 0],
                 [0, 0, 0, 0, 3]],
                numpy.int32))

    def test_aggregate_at_coords(self):

        # test grid
        grid = GridLatLon(
            axis_lat=GridAxis(start=-2.0, increment=2.0, length=2),
            axis_lon=GridAxis(start=-2.5, increment=1.0, length=5))

        # test data, including edge cases
        # grid box (row,col):     (0,4) (0,3) (1,4) (1,4) (1,0)
        # grid flat index   :       4     3     9     9     5
        coord_lat = numpy.array([-0.2, -1.2,  1.8,  1.7,  0.0], numpy.float32)
        coord_lon = numpy.array([1.6,  1.4,  2.3,  1.5, -2.5], numpy.float32)

        # aggregation object
        a = AggregateGridCount(grid.dimensions)

        # process
        a.aggregate_at_indices(grid.compute_indices(coord_lat, coord_lon))

        # check count
        numpy.testing.assert_equal(
            actual=a.count,
            desired=numpy.array(
                [[0, 0, 0, 1, 1],
                 [1, 0, 0, 0, 2]]))

        # add some data
        # grid box (row,col):         (0,4) (0,2) (1,4)
        # grid flat index   :           4     0     9
        add_coord_lat = numpy.array([-0.2, -2.0,  1.8], numpy.float32)
        add_coord_lon = numpy.array([1.6, -2.5,  2.3], numpy.float32)

        # process increment
        a.aggregate_at_indices(grid.compute_indices(
            add_coord_lat, add_coord_lon))

        # check count
        numpy.testing.assert_equal(
            actual=a.count,
            desired=numpy.array(
                [[1, 0, 0, 1, 2],
                 [1, 0, 0, 0, 3]]))


class TestAggregateGridSum(unittest.TestCase):

    def test_init(self):
        result = AggregateGridSum((88, 360))
        self.assertEqual((88, 360), result.count.shape)
        self.assertEqual(numpy.int32, result.count.dtype)
        self.assertEqual((88, 360), result.sum.shape)
        self.assertEqual(numpy.float32, result.sum.dtype)

    def test_aggregate_at_indices(self):

        # practice grid of 3 rows x 5 columns (starts and increments not
        # important)
        grid = GridLatLon(
            axis_lat=GridAxis(start=-20000.0, increment=100.0, length=3),
            axis_lon=GridAxis(start=250000.0, increment=1000.0, length=5))

        # aggregation object
        a = AggregateGridSum(grid.dimensions)

        # the indices to accumulate
        indices = numpy.array([0, 0, 6, 14, 14, 14], numpy.int32)

        # values to accumulate
        values = numpy.array(
            [6.0, 8.0, -329.875, -0.9, 0.3, 0.3], numpy.float32)

        # process
        a.aggregate_at_indices(indices, values)

        # check result
        numpy.testing.assert_almost_equal(
            actual=a.sum,
            desired=numpy.array(
                [[14.0,     0.0,  0.0, 0.0, 0.0],
                 [0.0, -329.875, 0.0, 0.0, 0.0],
                 [0.0,    0.0,   0.0, 0.0, -0.3]],
                numpy.float32),
            decimal=7)

    def test_aggregate_at_coords(self):

        # test grid
        grid = GridLatLon(
            axis_lat=GridAxis(start=-2.0, increment=2.0, length=2),
            axis_lon=GridAxis(start=-2.5, increment=1.0, length=5))

        # aggregation object
        a = AggregateGridSum(grid.dimensions)

        # test data, including edge cases
        # grid box (row,col):     (0,4) (0,3) (1,4) (1,4) (1,0)
        # grid flat index   :       4     3     9     9     5
        coord_lat = numpy.array([-0.2, -1.2,  1.8,  1.7,  0.0], numpy.float32)
        coord_lon = numpy.array([1.6,  1.4,  2.3,  1.5, -2.5], numpy.float32)
        values = numpy.array([-10.4, 29.88, -3.0, -4.0,  1.1], numpy.float32)

        # process
        a.aggregate_at_indices(grid.compute_indices(
            coord_lat, coord_lon), values)

        # check values

        numpy.testing.assert_almost_equal(
            actual=a.sum,
            desired=numpy.array(
                [[0.0, 0.0, 0.0, 29.88, -10.4],
                 [1.1, 0.0, 0.0,  0.0,   -7.0]],
                numpy.float32),
            decimal=6)

        numpy.testing.assert_equal(
            actual=a.count,
            desired=numpy.array(
                [[0, 0, 0, 1, 1],
                 [1, 0, 0, 0, 2]]))

        # add some data
        # grid box (row,col):         (0,4) (0,2) (1,4)
        # grid flat index   :           4     0     9
        add_coord_lat = numpy.array([-0.2, -2.0,  1.8], numpy.float32)
        add_coord_lon = numpy.array([1.6, -2.5,  2.3], numpy.float32)
        add_values = numpy.array([0.1, 88.8,  0.5], numpy.float32)

        # process increment
        a.aggregate_at_indices(grid.compute_indices(
            add_coord_lat, add_coord_lon), add_values)

        # check added values

        numpy.testing.assert_almost_equal(
            actual=a.sum,
            desired=numpy.array(
                [[88.8, 0.0, 0.0, 29.88, -10.3],
                 [1.1, 0.0, 0.0,  0.0,   -6.5]],
                numpy.float32),
            decimal=6)

        numpy.testing.assert_equal(
            actual=a.count,
            desired=numpy.array(
                [[1, 0, 0, 1, 2],
                 [1, 0, 0, 0, 3]]))

    def test_get_mean(self):

        # temp file (usually in /tmp folder): will be destroyed when
        # this variable goes out of scope (even if tests fail)
        outputfile = tempfile.NamedTemporaryFile(
            prefix='eustace.', suffix='.nc')

        # make a tiny test grid as before and aggregate
        grid = GridLatLon(
            axis_lat=GridAxis(start=-2.0, increment=2.0, length=2),
            axis_lon=GridAxis(start=-2.5, increment=1.0, length=5))

        # aggregation object
        a = AggregateGridSum(grid.dimensions)

        # grid box (row,col):     (0,4) (0,3) (1,4) (1,4) (1,0)
        # grid flat index   :       4     3     9     9     5
        test_lat = numpy.array([-0.2, -1.2,  1.8,  1.7,  0.0], numpy.float32)
        test_lon = numpy.array([1.6,  1.4,  2.3,  1.5, -2.5], numpy.float32)
        test_values = numpy.array([1.0,  3.0,  7.0, 13.0,  9.0], numpy.float32)
        a.aggregate_at_indices(observations=test_values, indices=grid.compute_indices(
            coord_lat=test_lat, coord_lon=test_lon))
        # compute mean
        values = a.get_mean(flag_invalid=numpy.float32(-1000.0))

        # save
        # grid.save(outputfile.name, 'testresult')
        # load in again
        # r = Dataset(outputfile.name)
        # lat = r.variables['latitude'][:]
        # lon = r.variables['longitude'][:]
        # values = r.variables['testresult'][:]

        # check mean is correct
        # the -1000 figures are for grid boxes without data (fill value for invalid items)
        # numpy.testing.assert_almost_equal(lat, [ -1.0, 1.0 ])
        # numpy.testing.assert_almost_equal(lon, [ -2.0, -1.0, 0.0, 1.0, 2.0 ])
        numpy.testing.assert_almost_equal(
            values,
            numpy.array([[-1000.0, -1000.0, -1000.0,     3.0,   1.0],
                         [9.0, -1000.0, -1000.0, -1000.0,  10.0]], numpy.float32))

        # file is automatically deleted anyway
        # but this avoids a style error due to unused variable
        del outputfile


class TestAggregateGridSumDevSq(unittest.TestCase):

    def test_init(self):
        result = AggregateGridSumDevSq(dimensions=(34, 22))
        self.assertEqual((34, 22), result.dimensions)
        self.assertEqual((34, 22), result.sum_devsq.shape)
        self.assertEqual(numpy.float32, result.sum_devsq.dtype)

    def test_aggregate_at_indices(self):

        a = AggregateGridSumDevSq((2, 3))

        mean = numpy.array(
            [[5.0, 6.0, 9.0],
             [1.0, 8.0, 4.0]], numpy.float32)

        # indices to place
        indices = numpy.array([1, 3, 1, 3, 3], numpy.int32)

        # values to accumulate into devsq
        values = numpy.array([3.0, -2.0, 8.0, 7.0, -1.0], numpy.float32)

        # do aggregation
        a.aggregate_at_indices(mean, indices, values)

        # check
        # element 1: (3-6)^2 + (8-6)^2 = 3^2 + 2^2 = 13
        # element 3: (-2-1)^2 + (7-1)^2 + (-1-1)^2 = 3^2 + 6^2 + 2^2 = 49
        numpy.testing.assert_almost_equal(
            numpy.array([[0.0, 13.0, 0.0],
                         [49.0,  0.0, 0.0]]),
            a.sum_devsq)

        # repeat aggregation
        a.aggregate_at_indices(mean, indices, values)

        # should double
        numpy.testing.assert_almost_equal(
            numpy.array([[0.0, 26.0, 0.0],
                         [98.0,  0.0, 0.0]]),
            a.sum_devsq)


class TestAggregateGridArray(unittest.TestCase):

    def test_init(self):
        result = AggregateGridArray(dimensions=(88, 360), max_bin_obs=43)
        self.assertEqual((88, 360), result.count.shape)
        self.assertEqual(numpy.int32, result.count.dtype)
        self.assertEqual((88, 360, 43), result.values.shape)
        self.assertEqual(numpy.float32, result.values.dtype)
        self.assertEqual((88, 360), result.uncertainty_random_sum_sq.shape)
        self.assertEqual((88, 360), result.uncertainty_local_atm_sum_sq.shape)
        self.assertEqual((88, 360), result.uncertainty_local_sfc_sum_sq.shape)
        self.assertEqual((88, 360), result.uncertainty_systematic_sum_sq.shape)

    def test_compute_fields(self):

        # temp file (usually in /tmp folder): will be destroyed when
        # this variable goes out of scope (even if tests fail)
        outputfile = tempfile.NamedTemporaryFile(
            prefix='eustace.', suffix='.nc')

        # make a tiny test grid as before and aggregate
        grid = GridLatLon(
            axis_lat=GridAxis(start=-2.0, increment=2.0, length=2),
            axis_lon=GridAxis(start=-2.5, increment=1.0, length=5))

        # aggregation object under test
        a = AggregateGridArray(dimensions=grid.dimensions, max_bin_obs=5)

        # grid box (row,col):     (0,4) (0,3) (1,4) (1,4) (1,0)
        # grid flat index   :       4     3     9     9     5
        test_lat = numpy.array([-0.2, -1.2,  1.8,  1.7,  0.0], numpy.float32)
        test_lon = numpy.array([1.6,  1.4,  2.3,  1.5, -2.5], numpy.float32)
        test_values = numpy.array([1.0,  3.0,  7.0, 13.0,  9.0], numpy.float32)
        test_uncertainty_systematic = numpy.array([1.0], numpy.float32)
        test_uncertainty_random = numpy.array([3.5, 1.0, 2.5, 9.3, 1.0], numpy.float32)
        test_uncertainty_local_atm = numpy.array([4.0, 6.5, 1.2, 1.1, 3.9], numpy.float32)
        test_uncertainty_local_sfc = numpy.array([8.0,13.0, 2.4, 2.2, 7.8], numpy.float32)
        a.aggregate_at_indices(
            observations=test_values, indices=grid.compute_indices(test_lat, test_lon),
            uncertainty_systematic=test_uncertainty_systematic,
            uncertainty_random=test_uncertainty_random,
            uncertainty_local_atm=test_uncertainty_local_atm,
            uncertainty_local_sfc=test_uncertainty_local_sfc)

        # compute output fields (dictionary of results)
        fields = a.compute_fields(flag_invalid=numpy.float32(-1000.0))

        # check mean is correct
        # the -1000 figures are for grid boxes without data (fill value for
        # invalid items)
        numpy.testing.assert_almost_equal(
            fields['tsmean'],
            numpy.array([[-1000.0, -1000.0, -1000.0,     3.0,   1.0],
                         [9.0, -1000.0, -1000.0, -1000.0,  10.0]], numpy.float32))

        # min
        numpy.testing.assert_almost_equal(
            fields['tsmin'],
            numpy.array([[-1000.0, -1000.0, -1000.0,     3.0,   1.0],
                         [9.0, -1000.0, -1000.0, -1000.0,   7.0]], numpy.float32))

        # max
        numpy.testing.assert_almost_equal(
            fields['tsmax'],
            numpy.array([[-1000.0, -1000.0, -1000.0,     3.0,   1.0],
                         [9.0, -1000.0, -1000.0, -1000.0,  13.0]], numpy.float32))

        # number of observations
        numpy.testing.assert_equal(
            fields['ts_number_of_observations'],
            numpy.array([[0, 0, 0, 1, 1],
                         [1, 0, 0, 0, 2]], numpy.int32))

        # file is automatically deleted anyway
        # but this avoids a style error due to unused variable
        del outputfile


class TestAggregateGridCountSumMinMax(unittest.TestCase):

    def test_init(self):
        result = AggregateGridCountSumMinMax((88, 360))
        self.assertEqual((88, 360), result.count.shape)
        self.assertEqual(numpy.int32, result.count.dtype)
        self.assertEqual((88, 360), result.values_sum.shape)
        self.assertEqual(numpy.float32, result.values_sum.dtype)
        self.assertEqual((88, 360), result.count.shape)
        self.assertEqual(numpy.int32, result.count.dtype)
        self.assertEqual((88, 360), result.values_min.shape)
        self.assertEqual(numpy.float32, result.values_min.dtype)
        self.assertEqual((88, 360), result.values_min.shape)
        self.assertEqual(numpy.float32, result.values_min.dtype)

    def test_compute_fields(self):

        # temp file (usually in /tmp folder): will be destroyed when
        # this variable goes out of scope (even if tests fail)
        outputfile = tempfile.NamedTemporaryFile(
            prefix='eustace.', suffix='.nc')

        # make a tiny test grid as before and aggregate
        grid = GridLatLon(
            axis_lat=GridAxis(start=-2.0, increment=2.0, length=2),
            axis_lon=GridAxis(start=-2.5, increment=1.0, length=5))

        # aggregation object under test
        a = AggregateGridCountSumMinMax(grid.dimensions)

        # grid box (row,col):     (0,4) (0,3) (1,4) (1,4) (1,0)
        # grid flat index   :       4     3     9     9     5
        test_lat = numpy.array([-0.2, -1.2,  1.8,  1.7,  0.0], numpy.float32)
        test_lon = numpy.array([1.6,  1.4,  2.3,  1.5, -2.5], numpy.float32)
        test_values = numpy.array([1.0,  3.0,  7.0, 13.0,  9.0], numpy.float32)
        test_uncertainty_systematic = numpy.array([1.0], numpy.float32)
        test_uncertainty_random = numpy.array(
            [3.5, 1.0, 2.5, 9.3, 1.0], numpy.float32)
        test_uncertainty_local_atm = numpy.array([4.0, 6.5, 1.2, 1.1, 3.9], numpy.float32)
        test_uncertainty_local_sfc = numpy.array([8.0,13.0, 2.4, 2.2, 7.8], numpy.float32)
        a.aggregate_at_indices(
            observations=test_values, indices=grid.compute_indices(
                coord_lat=test_lat, coord_lon=test_lon),
            uncertainty_systematic=test_uncertainty_systematic,
            uncertainty_random=test_uncertainty_random,
            uncertainty_local_atm=test_uncertainty_local_atm,
            uncertainty_local_sfc=test_uncertainty_local_sfc)

        # compute output fields (dictionary of results)
        fields = a.compute_fields(flag_invalid=numpy.float32(-1000.0))

        # check mean is correct
        # the -1000 figures are for grid boxes without data (fill value for
        # invalid items)
        numpy.testing.assert_almost_equal(
            fields['tsmean'],
            numpy.array([[-1000.0, -1000.0, -1000.0,     3.0,   1.0],
                         [9.0, -1000.0, -1000.0, -1000.0,  10.0]], numpy.float32))

        # min
        numpy.testing.assert_almost_equal(
            fields['tsmin'],
            numpy.array([[-1000.0, -1000.0, -1000.0,     3.0,   1.0],
                         [9.0, -1000.0, -1000.0, -1000.0,   7.0]], numpy.float32))

        # max
        numpy.testing.assert_almost_equal(
            fields['tsmax'],
            numpy.array([[-1000.0, -1000.0, -1000.0,     3.0,   1.0],
                         [9.0, -1000.0, -1000.0, -1000.0,  13.0]], numpy.float32))

        # number of observations
        numpy.testing.assert_equal(
            fields['ts_number_of_observations'],
            numpy.array([[0, 0, 0, 1, 1],
                         [1, 0, 0, 0, 2]], numpy.int32))

        # file is automatically deleted anyway
        # but this avoids a style error due to unused variable
        del outputfile
