"""Tests for definition of epoch."""
# pylint: disable=missing-docstring

import unittest
import numpy
from eustace.timeutils.epoch import EPOCH
from eustace.timeutils.epoch import days_since_epoch
from eustace.timeutils.epoch import epoch_plus_days
from eustace.timeutils.epoch import seconds_since_epoch
from eustace.timeutils.epoch import epoch_plus_seconds
from datetime import datetime


class TestEpoch(unittest.TestCase):

    def test_epoch(self):
        self.assertEqual(datetime(1850, 1, 1), EPOCH)

    def test_days_since_epoch_none_fails(self):
        self.assertRaises(TypeError, days_since_epoch, None)

    def test_days_since_epoch_zero(self):
        result = days_since_epoch(datetime(1850, 1, 1))
        self.assertTrue(isinstance(result, float))
        self.assertAlmostEqual(0.0, result)

    def test_days_since_epoch_negative(self):
        self.assertAlmostEqual(-1.0, days_since_epoch(datetime(1849, 12, 31)))

    def test_days_since_epoch_20151105(self):
        result = days_since_epoch(datetime(2015, 11, 05))
        self.assertTrue(isinstance(result, float))
        self.assertAlmostEqual(60573, result)

    def test_days_since_epoch_20151105_153247(self):
        result = days_since_epoch(datetime(2015, 11, 05, 15, 32, 47))
        self.assertTrue(isinstance(result, float))
        self.assertAlmostEqual(60573.647766204, result)

    def test_epoch_plus_days(self):
        result = epoch_plus_days(60573.647766204)
        self.assertTrue(isinstance(result, datetime))
        self.assertEqual(2015, result.year)
        self.assertEqual(11, result.month)
        self.assertEqual(5, result.day)
        self.assertEqual(15, result.hour)
        self.assertEqual(32, result.minute)
        self.assertEqual(47, result.second)
        # Note microsecond is not tested
        # - this is not precise to microsecond resolution

    def test_seconds_since_epoch_none_fails(self):
        self.assertRaises(TypeError, seconds_since_epoch, None)

    def test_seconds_since_epoch_zero(self):
        result = seconds_since_epoch(datetime(1850, 1, 1))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertEqual(0, result)

    def test_seconds_since_epoch_negative(self):
        result = seconds_since_epoch(datetime(1849, 12, 31))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertEqual(-86400, result)

    def test_seconds_since_epoch_20151105(self):
        result = seconds_since_epoch(datetime(2015, 11, 05))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertEqual(5233507200, result)

    def test_seconds_since_epoch_20151105_153247(self):
        result = seconds_since_epoch(datetime(2015, 11, 05, 15, 32, 47))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertEqual(5233563167, result)

    def test_epoch_plus_seconds(self):
        result = epoch_plus_seconds(numpy.uint64(5233563167))
        self.assertTrue(isinstance(result, datetime))
        self.assertEqual(datetime(2015, 11, 05, 15, 32, 47), result)
