"""Tests for underlying filebuilder class."""

import unittest
import tempfile
import numpy
from ..filebuilder import FileBuilder
from ..outputvariable import OutputVariable
import netCDF4

class TestFileBuilder(unittest.TestCase):

    def test_all_methods(self):

        # Test data
        testfile = tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.test_filebuilder.', suffix='.nc')
        testdata_values = numpy.array(
            [ [ 5.0, 2.1, 9.0, 22.0, 1.3 ],
              [ 1.1,99.9, 3.0,  1.0, 2.2 ],
              [23.1, 1.8, 8.8,  1.2, 0.8 ] ])
        testdata_mask = numpy.array(
            [ [ False, True, False, False, False ],
              [ False, False, False, False, False ],
              [ True, False, False, False, False ] ])
        testdata = numpy.ma.masked_array(data=testdata_values, mask=testdata_mask)

        # Make file
        builder = FileBuilder()
        builder.create(testfile.name, title='A Test', institution='My Place', comment='What?', history='Very old', source='Ketchup')
        builder.add_dimension('dozen', 12)
        builder.add_dimension_and_variable('y_axis', OutputVariable('y_axis', numpy.float32, None), numpy.arange(3))
        builder.add_dimension_and_variable('x_axis', OutputVariable('x_axis', numpy.float32, None), numpy.arange(5))
        builder.add_variable(OutputVariable('stuff', numpy.float32, None), ('y_axis', 'x_axis'), testdata)
        builder.save_and_close()

        # Check results
        result = netCDF4.Dataset(testfile.name, 'r')
        self.assertEqual('A Test', result.title)
        self.assertEqual('My Place', result.institution)
        self.assertEqual('What?', result.comment)
        self.assertEqual('Very old', result.history)
        self.assertEqual('Ketchup', result.source)
        self.assertEqual(12, len(result.dimensions['dozen']))
        numpy.testing.assert_equal(numpy.arange(3), result.variables['y_axis'][:])
        numpy.testing.assert_equal(numpy.arange(5), result.variables['x_axis'][:])
        numpy.testing.assert_almost_equal(result.variables['stuff'][:], testdata, decimal=4)
        numpy.testing.assert_equal(result.variables['stuff'][:].mask, testdata_mask)
