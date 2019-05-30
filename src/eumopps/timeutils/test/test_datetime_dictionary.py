"""Test representation as dictionary of strings."""

# pylint: disable=missing-docstring

import unittest
from .. import datetime_dictionary
from datetime import datetime


class TestDateTimeDictionary(unittest.TestCase):

    def test_keys(self):
        self.assertEqual('year', datetime_dictionary.YEAR)
        self.assertEqual('month', datetime_dictionary.MONTH)
        self.assertEqual('day', datetime_dictionary.DAY)
        self.assertEqual('hour', datetime_dictionary.HOUR)
        self.assertEqual('minute', datetime_dictionary.MINUTE)
        self.assertEqual('second', datetime_dictionary.SECOND)

    def test_key_list(self):
        self.assertEqual(6, len(datetime_dictionary.FIELDS))
        self.assertTrue('year' in datetime_dictionary.FIELDS)
        self.assertTrue('month' in datetime_dictionary.FIELDS)
        self.assertTrue('day' in datetime_dictionary.FIELDS)
        self.assertTrue('hour' in datetime_dictionary.FIELDS)
        self.assertTrue('minute' in datetime_dictionary.FIELDS)
        self.assertTrue('second' in datetime_dictionary.FIELDS)

    def test_parse(self):
        testvalues = {'year': '1852', 'month': '08', 'day': '22',
                      'hour': '18', 'minute': '53', 'second': '07'}
        result = datetime_dictionary.parse(testvalues)
        self.assertEqual(datetime(1852, 8, 22, 18, 53, 7), result)
