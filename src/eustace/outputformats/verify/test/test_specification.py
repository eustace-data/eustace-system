"""Tests of specification structures used as template for format verification."""

# pylint: disable=missing-docstring,invalid-name

__version__ = "$Revision: 550 $"
__author__ = "Joel R. Mitchelson"

import unittest
import numpy
from ..specification import VariableSpecification


class TestVariableSpecification(unittest.TestCase):

    def test_init(self):
        v = VariableSpecification('MyName', numpy.float64)
        self.assertEqual('MyName', v.name)
        self.assertEqual(numpy.float64, v.dtype)
        self.assertEqual({}, v.metadata)

    def test_add_metadata(self):
        v = VariableSpecification('FruitTree', numpy.float32)
        self.assertEqual('FruitTree', v.name)
        self.assertEqual(numpy.float32, v.dtype)
        self.assertEqual({}, v.metadata)
        returnvalue1 = v.add_metadata('Variety', 'Apple')
        returnvalue2 = v.add_metadata('Size', 'Massive')
        self.assertEqual('FruitTree', v.name)
        self.assertEqual(numpy.float32, v.dtype)
        self.assertEqual(2, len(v.metadata))
        self.assertEqual('Apple', v.metadata['Variety'])
        self.assertEqual('Massive', v.metadata['Size'])
        self.assertTrue(v == returnvalue1)
        self.assertTrue(v == returnvalue2)

    def test_str(self):
        v = VariableSpecification('Traffic', numpy.int32)
        v.add_metadata('Variety', 'Cessna')
        v.add_metadata('Size', 'OnlyLittle')
        self.assertEqual(
            "--Traffic <type 'numpy.int32'>\n    Size: OnlyLittle\n    Variety: Cessna", str(v))
