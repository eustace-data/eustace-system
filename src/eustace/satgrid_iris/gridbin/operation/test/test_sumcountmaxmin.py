"""Test the sum and count operation."""

import unittest
import numpy
from iris.cube import Cube
from iris.coords import AuxCoord
from ..sumcountmaxmin import GridBinOperationSumCountMaxMin
from ...gridfieldlist import GridFieldDescriptor
from ...gridfieldlist import GridFieldAxis
from ...gridbin import GridBins

class TestGridBinOperationSumCount(unittest.TestCase):

    def test_run(self):

        # make a tiny test grid
        grid = GridBins(
            [ GridFieldAxis(zeroth=-3.0, step=2.0, count=2, var_name='latitude', circular=False),
              GridFieldAxis(zeroth=-3.0, step=1.0, count=5, var_name='longitude', circular=False) ] )

        # make measurement values and flag final one as invalid
        measurements = numpy.ma.masked_array(
            data = [ [ 1, 2 ], [ 3, 4 ], [ 5, 6 ] ], 
            mask = [ [ False, False ], [ False, False ], [ False, True ] ],
            dtype=numpy.float32)
              
        # make measurement locations (including edge case)
        # grid box (row,col):     (0,4) (0,3) (1,4) (1,4) (1,0) -invalid-
        # grid flat index   :       4     3     9     9     5   -invalid-
        coord_lat = numpy.array([ [ -0.2, -1.2 ], [ 1.8,  1.7 ], [  0.0, 0.5 ] ], numpy.float32)  # noqa (allow multiple spaces)
        coord_lon = numpy.array([ [  1.6,  1.4 ], [ 2.3,  1.5 ], [ -2.5, 0.5 ] ], numpy.float32)  # noqa (allow multiple spaces)

        # cube with this info
        testcube = Cube(measurements, var_name='measurements')
        testcube.add_aux_coord(AuxCoord(coord_lon, var_name='longitude'), (0,1))
        testcube.add_aux_coord(AuxCoord(coord_lat, var_name='latitude'), (0,1))
        
        # get map
        result = grid.compute_input_map(testcube)

        # new sum operation
        sum_operator = GridBinOperationSumCountMaxMin()
        sum_operator.set_input_descriptor(
            GridBinOperationSumCountMaxMin.INPUT, 
            GridFieldDescriptor('measurements'))
        sum_operator.set_output_descriptor(
            GridBinOperationSumCountMaxMin.OUTPUTTOTAL, 
            GridFieldDescriptor('total', numpy.float32))
        sum_operator.set_output_descriptor(
            GridBinOperationSumCountMaxMin.OUTPUTCOUNT, 
            GridFieldDescriptor('count', numpy.int32))
        sum_operator.set_output_descriptor(
            GridBinOperationSumCountMaxMin.OUTPUTMAX, 
            GridFieldDescriptor('maximum', numpy.float32),
            initialvalue=-999)
        sum_operator.set_output_descriptor(
            GridBinOperationSumCountMaxMin.OUTPUTMIN, 
            GridFieldDescriptor('minimum', numpy.float32),
            initialvalue=999)

        # run once and check
        sum_operator.run(result)
        numpy.testing.assert_array_equal(grid[0].data, [ [ 0, 0, 0, 1, 1 ], [ 1, 0, 0, 0, 2] ])
        numpy.testing.assert_array_equal(grid[1].data, [ [ 0, 0, 0, 2, 1 ], [ 5, 0, 0, 0, 7] ])
        numpy.testing.assert_array_equal(grid[2].data, [ [ 0, 0, 0, 2, 1 ], [ 5, 0, 0, 0, 3] ])
        numpy.testing.assert_array_equal(grid[3].data, [ [ 0, 0, 0, 2, 1 ], [ 5, 0, 0, 0, 4] ])
        numpy.testing.assert_array_equal(grid[0].data.mask, 
                                         [ [ True, True, True, False, False ], 
                                           [ False, True, True, True, False ] ])
        numpy.testing.assert_array_equal(grid[1].data.mask, 
                                         [ [ True, True, True, False, False ], 
                                           [ False, True, True, True, False ] ])
        numpy.testing.assert_array_equal(grid[2].data.mask, 
                                         [ [ True, True, True, False, False ], 
                                           [ False, True, True, True, False ] ])
        numpy.testing.assert_array_equal(grid[3].data.mask, 
                                         [ [ True, True, True, False, False ], 
                                           [ False, True, True, True, False ] ])

        # run again and check it has doubled
        sum_operator.run(result)
        numpy.testing.assert_array_equal(grid[0].data, [ [ 0, 0, 0, 2, 2 ], [ 2, 0, 0, 0, 4] ])
        numpy.testing.assert_array_equal(grid[1].data, [ [ 0, 0, 0, 4, 2 ], [10, 0, 0, 0,14] ])
        numpy.testing.assert_array_equal(grid[0].data.mask, 
                                         [ [ True, True, True, False, False ], 
                                           [ False, True, True, True, False ] ])
        numpy.testing.assert_array_equal(grid[1].data.mask, 
                                         [ [ True, True, True, False, False ], 
                                           [ False, True, True, True, False ] ])

