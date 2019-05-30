"""Tests for Checksum computation."""

# pylint: disable=missing-docstring,invalid-name

import unittest
import os
import numpy
from ..checksum import Checksum


class TestChecksum(unittest.TestCase):

    def test_checksum(self):
        pathname = os.path.join(os.path.dirname(
            __file__), 'checksumexample.txt')
        result = Checksum(pathname)
        self.assertEqual('3469475163', result.checksum)
        self.assertEqual(1431, result.size)
        self.assertTrue(isinstance(result.size, numpy.int64))
