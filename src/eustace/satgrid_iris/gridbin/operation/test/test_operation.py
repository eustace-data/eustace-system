"""Test operation API."""

import unittest
from ..operation import GridBinOperation
from ...gridfieldlist import GridFieldDescriptor
from ...gridfieldlist import GridFieldAxis
from ...gridbin import GridBins
from ...gridbin import GridBinMap
from iris.cube import Cube
from iris.coords import AuxCoord
import numpy

class TestGridBinOperation(unittest.TestCase):

    def test_init(self):
        op = GridBinOperation([ 'onestat', 'twostats' ], [ 'apples', 'pears', 'bananas' ])
        self.assertEqual({'onestat':None,'twostats':None}, op.outputdescriptors)
        self.assertEqual({'apples':None,'pears':None,'bananas':None}, op.inputdescriptors)


    def test_set_output_descriptor(self):

        op = GridBinOperation(['somestat'], ['oneinput', 'twoinputs'])

        # Check with wrong variable name
        op = GridBinOperation(['somestat'], ['oneinput', 'twoinputs'])
        self.assertRaises(ValueError, GridBinOperation.set_output_descriptor, op, 'anotherstat', GridFieldDescriptor('mystat'))

        # Test with the descriptor of wrong type
        op = GridBinOperation(['somestat'], ['oneinput', 'twoinputs'])
        self.assertRaises(TypeError, GridBinOperation.set_output_descriptor, op, 'somestat', 'mystat')

        # Test setting
        op.set_output_descriptor('somestat', GridFieldDescriptor('measurements', long_name='Important stuff'))
        self.assertEqual(op.outputdescriptors['somestat'].var_name, 'measurements')
        self.assertEqual(op.outputdescriptors['somestat'].long_name, 'Important stuff')


    def test_set_output_descriptor(self):

        op = GridBinOperation(['somestat'], ['oneinput', 'twoinputs'])

        # Check with wrong variable name
        op = GridBinOperation(['somestat'], ['oneinput', 'twoinputs'])
        self.assertRaises(ValueError, GridBinOperation.set_input_descriptor, op, 'threeinputs', GridFieldDescriptor('myinput'))

        # Test with the descriptor of wrong type
        op = GridBinOperation(['somestat'], ['oneinput', 'twoinputs'])
        self.assertRaises(TypeError, GridBinOperation.set_input_descriptor, op, 'twoinputs', 'myinput')

        # Test setting
        op.set_input_descriptor('twoinputs', GridFieldDescriptor('measurements', long_name='Important stuff'))
        self.assertEqual(op.inputdescriptors['twoinputs'].var_name, 'measurements')
        self.assertEqual(op.inputdescriptors['twoinputs'].long_name, 'Important stuff')

    def test_run(self):

        # Only checks that the NotImplementedError is thrown in base class
        # and that detection of possible gaps in descriptor assignments work ok.
        # For check of underlying functionality see integration test in test_sum.

        # Axes to use for output grid
        bin_latitude = GridFieldAxis(zeroth=-50.0, step=5.0, count=2, var_name='latitude', circular='false')
        bin_longitude = GridFieldAxis(zeroth=-6.0, step=3.0, count=3, var_name='longitude', circular='false')

        # Example to bin onto
        testbins = GridBins([bin_latitude, bin_longitude])

        # Make measurement values and flag final one as invalid
        measurements = numpy.ma.masked_array(
            data = [ 23.4 ], 
            mask = [ False ],
            dtype=numpy.float32)              
        coord_lat = numpy.array([-47.5 ], numpy.float32)  # noqa (allow multiple spaces)
        coord_lon = numpy.array([ 0.25 ], numpy.float32)  # noqa (allow multiple spaces)

        # cube with this info
        testcube = Cube(measurements, var_name='measurements')
        testcube.add_aux_coord(AuxCoord(coord_lat, var_name='latitude'), (0))
        testcube.add_aux_coord(AuxCoord(coord_lon, var_name='longitude'), (0))

        # Example mapping
        testmap = testbins.compute_input_map(testcube)

        # Operation with info partially set
        op = GridBinOperation(['somestat', 'anotherstat'], ['oneinput', 'twoinputs'])
        op.set_input_descriptor('oneinput', GridFieldDescriptor('measurements'))
        op.set_output_descriptor('somestat', GridFieldDescriptor('partialresults'))

        # Should raise a ValueError due to partial info
        self.assertRaisesRegexp(ValueError,
                                'Missing descriptors for identifiers: \[\'twoinputs\', \'anotherstat\'\]', 
                                GridBinOperation.run, op, testmap)

        # Add all inputs (should still fail on outputs)
        op.set_input_descriptor('twoinputs', GridFieldDescriptor('measurements'))
        self.assertRaisesRegexp(ValueError,
                                'Missing descriptors for identifiers: \[\'anotherstat\'\]', 
                                GridBinOperation.run, op, testmap)

        # Add all outputs
        op.set_output_descriptor('anotherstat', GridFieldDescriptor('moreresults'))

        # Should now just fail because the operation isn't implemented in base class
        self.assertRaises(NotImplementedError, GridBinOperation.run, op, testmap)
