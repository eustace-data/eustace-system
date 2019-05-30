"""Tests for timebase offsets."""

# pylint: disable=missing-docstring,invalid-name

import unittest
import numpy
from datetime import datetime
from ..timebase import TimeBaseDays
from ..timebase import TimeBaseSeconds


class TestTimeBaseDays(unittest.TestCase):

    def test_init(self):
        t = TimeBaseDays(datetime(1952, 3, 22))
        self.assertEqual(datetime(1952, 3, 22), t.epoch)

    def test_datatype(self):
        self.assertEqual(TimeBaseDays(
            datetime(1850, 1, 1)).datatype(), numpy.float32)

    def test_units(self):
        self.assertEqual('days since 1827-03-22 23:29:07 UTC',
                         TimeBaseDays(datetime(1827, 3, 22, 23, 29, 7)).units())

    def test_days_since_epoch_none_fails(self):
        self.assertRaises(TypeError, TimeBaseDays.datetime_to_number,
                          TimeBaseDays(datetime(1970, 1, 1)), None)

    def test_days_since_epoch_zero(self):
        result = TimeBaseDays(datetime(1842, 4, 21)).datetime_to_number(
            datetime(1842, 4, 21))
        self.assertTrue(isinstance(result, float))
        self.assertAlmostEqual(0.0, result)

    def test_days_since_epoch_negative(self):
        self.assertAlmostEqual(-1.0, TimeBaseDays(datetime(1850, 1, 1)).datetime_to_number(datetime(1849, 12, 31)))

    def test_days_since_epoch_20151105(self):
        result = TimeBaseDays(datetime(1850, 1, 1)).datetime_to_number(datetime(2015, 11, 05))
        self.assertTrue(isinstance(result, float))
        self.assertAlmostEqual(60573, result)

    def test_days_since_epoch_20151105_153247(self):
        result = TimeBaseDays(datetime(1850, 1, 1)).datetime_to_number(datetime(2015, 11, 05, 15, 32, 47))
        self.assertTrue(isinstance(result, float))
        self.assertAlmostEqual(60573.647766204, result)

    def test_epoch_plus_days(self):
        result = TimeBaseDays(datetime(1850, 1, 1)).number_to_datetime(60573.647766204)
        self.assertTrue(isinstance(result, datetime))
        self.assertEqual(2015, result.year)
        self.assertEqual(11, result.month)
        self.assertEqual(5, result.day)
        self.assertEqual(15, result.hour)
        self.assertEqual(32, result.minute)
        self.assertEqual(47, result.second)
        # Note microsecond is not tested
        # - this is not precise to microsecond resolution


class TestTimeBaseSeconds(unittest.TestCase):

    def test_init(self):
        t = TimeBaseSeconds(datetime(1952, 3, 22, 11, 18, 52))
        self.assertEqual(datetime(1952, 3, 22, 11, 18, 52), t.epoch)

    def test_datatype(self):
        self.assertEqual(TimeBaseSeconds(
            datetime(1850, 1, 1)).datatype(), numpy.int64)

    def test_units(self):
        self.assertEqual('seconds since 1827-03-22 23:29:07 UTC',
                         TimeBaseSeconds(datetime(1827, 3, 22, 23, 29, 7)).units())

    def test_days_since_epoch_none_fails(self):
        self.assertRaises(TypeError, TimeBaseDays.datetime_to_number,
                          TimeBaseSeconds(datetime(1970, 1, 1)), None)

    def test_days_since_epoch_zero(self):
        result = TimeBaseSeconds(datetime(1842, 4, 21, 12, 13, 14)).datetime_to_number(datetime(1842, 4, 21, 12, 13, 14))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertAlmostEqual(0, result)

    def test_seconds_since_epoch_negative(self):
        result = TimeBaseSeconds(datetime(1850, 1, 1)).datetime_to_number(datetime(1849, 12, 31))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertEqual(-86400, result)

    def test_seconds_since_epoch_20151105(self):
        result = TimeBaseSeconds(datetime(1850, 1, 1)).datetime_to_number(datetime(2015, 11, 05))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertEqual(5233507200, result)

    def test_seconds_since_epoch_20151105_153247(self):
        result = TimeBaseSeconds(datetime(1850, 1, 1)).datetime_to_number(datetime(2015, 11, 05, 15, 32, 47))
        self.assertTrue(isinstance(result, numpy.int64))
        self.assertEqual(5233563167, result)

    def test_epoch_plus_seconds(self):
        result = TimeBaseSeconds(datetime(1850, 1, 1)).number_to_datetime(numpy.uint64(5233563167))
        self.assertTrue(isinstance(result, datetime))
        self.assertEqual(datetime(2015, 11, 05, 15, 32, 47), result)
