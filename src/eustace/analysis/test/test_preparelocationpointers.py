"""Tests for (some parts of) location pointer preparation."""

import unittest
from ..preparelocationpointers import extractdaynumber

class TestExtractDayNumber(unittest.TestCase):

    def test_extract_daynumber(self):

        self.assertEqual(44810, extractdaynumber('bob_19720908.bin'))
        self.assertEqual(3, extractdaynumber('bob_18500104.bin'))

    def test_bad_date_fails(self):

        with self.assertRaises(ValueError):

            extractdaynumber('bob_ABCDEFGH.bin')

    def test_bad_extension_fails(self):

        with self.assertRaises(ValueError):

            extractdaynumber('bob_19720908.who')
