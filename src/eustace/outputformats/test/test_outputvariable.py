"""Checks that OutputVariable and OutputVariableTemplate construction behaves as expected."""

import unittest
import numpy

from ..outputvariable import OutputVariableTemplate
from ..outputvariable import OutputVariable

class TestOutputVariableTemplate(unittest.TestCase):

    def test_init(self):
        
        v = OutputVariableTemplate(numpy.int16, -3200, units='mm')
        self.assertEqual(numpy.int16, v.dtype)
        self.assertEqual(-3200, v.fill_value)
        self.assertEqual('mm', v.units)

    def test_extend(self):

        v = OutputVariableTemplate(numpy.int16, -3200, units='mm')
        w = OutputVariableTemplate.extend(v, standard_name='Bob', scale_factor=3.25, length_scale=4,length_scale_units='km')
        self.assertEqual(numpy.int16, w.dtype)
        self.assertEqual(-3200, w.fill_value)
        self.assertEqual('mm', w.units)
        self.assertEqual('Bob', w.standard_name)
        self.assertEqual(3.25, w.scale_factor)
        self.assertEqual(4, w.length_scale)
        self.assertEqual('km', w.length_scale_units)



class TestOutputVariable(unittest.TestCase):

    def test_from_template(self):

        t = OutputVariableTemplate(numpy.int16, -3200, units='mm')
        v = OutputVariable.from_template(t, 'bob', 'Bob\'s House', ancillary_variables='Bob\'s Shed')
        self.assertEqual(numpy.int16, v.dtype)
        self.assertEqual(-3200, v.fill_value)
        self.assertEqual('mm', v.units)
        self.assertEqual('bob', v.name)
        self.assertEqual('Bob\'s House', v.long_name)
        self.assertEqual('Bob\'s Shed', v.ancillary_variables)

    def test_from_template_with_long_name(self):

        t = OutputVariableTemplate(numpy.int16, -3200, '{person}\'s {dwelling} near {surrounds}.', units='mm')
        v = OutputVariable.from_template(t, 'bob', person='bob', dwelling='Cottage', surrounds='Exeter', ancillary_variables='Bob\'s Shed')
        self.assertEqual(numpy.int16, v.dtype)
        self.assertEqual(-3200, v.fill_value)
        self.assertEqual('mm', v.units)
        self.assertEqual('bob', v.name)
        self.assertEqual('Bob\'s Cottage near Exeter.', v.long_name)
        self.assertEqual('Bob\'s Shed', v.ancillary_variables)
