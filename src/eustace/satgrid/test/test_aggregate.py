"""Unit tests for native aggregation code."""
# pylint: disable=missing-docstring, invalid-name, bad-whitespace

import unittest
import numpy
from eustace.satgrid import aggregate


class TestAggregate(unittest.TestCase):

    def test_aggregate_sum(self):

        result = numpy.zeros((3, 2), numpy.float32)

        a = numpy.array([[1.0, 2.0], [3.0, 2.0]], numpy.float32)
        b = numpy.array([[5,  1], [1,  3]], numpy.int32)

        aggregate.aggregate_sum(result, a, b)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[0.0, 5.0],
                 [0.0, 2.0],
                 [0.0, 1.0]], numpy.float32),
            result)

        aggregate.aggregate_sum(result, a, b)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[0.0, 10.0],
                 [0.0,  4.0],
                 [0.0,  2.0]], numpy.float32),
            result)

    def test_aggregate_constant(self):

        result = numpy.zeros((3, 2), numpy.float32)

        aggregate.aggregate_constant(result, numpy.float32(
            2.5), numpy.array([1, 1, 1, 0, 4], numpy.int32))

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[2.5,  7.5],
                 [0.0,  0.0],
                 [2.5,  0.0]], numpy.float32),
            result)

    def test_aggregate_count_sum_min_max(self):

        result_count = numpy.zeros((3, 2), numpy.int32)
        result_sum = numpy.zeros((3, 2), numpy.float32)
        result_min = numpy.ones((3, 2), numpy.float32) * numpy.float32(1.0E12)
        result_max = numpy.ones((3, 2), numpy.float32) * numpy.float32(-1.0E12)

        a = numpy.array([[1.0, 2.0], [3.0, 2.0]], numpy.float32)
        b = numpy.array([[5,  1], [1,  3]], numpy.int32)
        c = numpy.array([23.0, -2.0, 4.0, 0.33], numpy.float32)
        d = numpy.array([0,   1,  2,  4], numpy.int32)

        # first aggregation
        aggregate.aggregate_count_sum_min_max(
            result_count, result_sum, result_min, result_max, a, b)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[0.0, 5.0],
                 [0.0, 2.0],
                 [0.0, 1.0]], numpy.float32),
            result_sum)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[0, 2],
                 [0, 1],
                 [0, 1]], numpy.int32),
            result_count)

        # should only be affected at non-empty bins
        numpy.testing.assert_almost_equal(
            numpy.array(
                [[1.0E12, 2.0],
                 [1.0E12, 2.0],
                 [1.0E12, 1.0]], numpy.float32),
            result_min)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1.0E12, 3.0],
                 [-1.0E12, 2.0],
                 [-1.0E12, 1.0]], numpy.float32),
            result_max)

        # second aggregation
        aggregate.aggregate_count_sum_min_max(
            result_count, result_sum, result_min, result_max, c, d)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[23.0,  3.0],
                 [4.0,  2.0],
                 [0.33, 1.0]], numpy.float32),
            result_sum)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[1, 3],
                 [1, 1],
                 [1, 1]], numpy.int32),
            result_count)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[23.0, -2.0],
                 [4.0, 2.0],
                 [0.33, 1.0]], numpy.float32),
            result_min)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[23.0, 3.0],
                 [4.0, 2.0],
                 [0.33, 1.0]], numpy.float32),
            result_max)

    def test_aggregate_sum_dev_sq(self):

        result = numpy.zeros((2, 3), numpy.float32)
        mean = numpy.array([[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]], numpy.float32)
        values0 = numpy.array([[2.0, 3.0], [4.0, 5.0]], numpy.float32)
        indices0 = numpy.array([3, 1, 3, 3], numpy.int32)
        values1 = numpy.array([[7.0, 8.0]], numpy.float32)
        indices1 = numpy.array([2, 1], numpy.int32)

        # first aggregation
        aggregate.aggregate_sum_dev_sq(result, mean, values0, indices0)

        # element 1: (3-2)^2 = 1
        # element 3: (2-1)^2 + (4-1)^2 + (5-1)^2 = 1 + 9 + 16 = 26
        numpy.testing.assert_almost_equal(
            numpy.array([[0.0, 1.0, 0.0], [26.0, 0.0, 0.0]], numpy.float32),
            result)

        # second aggregation
        aggregate.aggregate_sum_dev_sq(result, mean, values1, indices1)

        # element 1: 1 + (8-2)^2 = 37
        # element 2: 0 + (7-3)^2 = 16
        # element 3: 26 (as before)
        numpy.testing.assert_almost_equal(
            numpy.array([[0.0, 37.0, 16.0], [26.0, 0.0, 0.0]], numpy.float32),
            result)

    def test_aggregate_count(self):

        result = numpy.zeros((2, 4), numpy.int32)
        aggregate.aggregate_count(result, numpy.array(
            [3, 2, 3, 3, 4, 7], numpy.int32))
        numpy.testing.assert_equal(result, numpy.array(
            [[0, 0, 1, 3], [1, 0, 0, 1]], numpy.int32))
        aggregate.aggregate_count(result, numpy.array([3], numpy.int32))
        numpy.testing.assert_equal(result, numpy.array(
            [[0, 0, 1, 4], [1, 0, 0, 1]], numpy.int32))
