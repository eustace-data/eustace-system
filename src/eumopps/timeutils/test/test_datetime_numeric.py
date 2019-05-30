"""Test expression of datetime as string of numerals."""

# pylint: disable=missing-docstring

import unittest
from .. import datetime_numeric
from datetime import datetime


class TestDateTimeNumeric(unittest.TestCase):

    def test_parse(self):
        self.assertEqual(datetime(1873, 03, 19, 15, 45), datetime_numeric.parse('18730319154500'))
        self.assertRaises(TypeError, datetime_numeric.parse, None)
        self.assertRaises(ValueError, datetime_numeric.parse, 'Bob')

    def test_build(self):
        self.assertEqual('18540908203201', datetime_numeric.build(datetime(1854, 9, 8, 20, 32, 1)))

    def test_build_from_pattern(self):

        self.assertEqual('some_02_time_03_1854_ago',
                         datetime_numeric.build_from_pattern('some_%d_time_%m_%Y_ago', datetime(1854, 3, 2)))

        self.assertEqual('Date 1982-03-01 Time 14:32:19 Go!',
                         datetime_numeric.build_from_pattern('Date %Y-%m-%d Time %H:%M:%S Go!', datetime(1982, 3, 1, 14, 32, 19)))
                         
