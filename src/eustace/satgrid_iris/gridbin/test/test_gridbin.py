"""Tests for gridbin."""

import unittest
import numpy
from iris.cube import Cube
from iris.coords import AuxCoord
from iris.coords import DimCoord
from ..gridfieldlist import GridFieldAxis
from ..gridfieldlist import GridFieldDescriptor
from ..gridbin import GridBinMap
from ..gridbin import GridBins

class TestGridBinMap(unittest.TestCase):

    def test_init(self):
        grid = GridBins([ 
                GridFieldAxis(zeroth=-0.5, step=1.0, count=8, var_name='a', circular=False),
                GridFieldAxis(zeroth=-0.5, step=1.0, count=8, var_name='b', circular=False) ])
        result = GridBinMap(grid, Cube(numpy.array([ 0 ]), var_name='What?'), [ 2, 3 ], [ 7, 5 ])
        self.assertEqual('a', result.outputgrid.axes[0].var_name)
        self.assertEqual('b', result.outputgrid.axes[1].var_name)
        self.assertEqual('What?', result.inputcube.var_name)
        self.assertEqual([ 2, 3 ], result.mapfrom)
        self.assertEqual([ 7, 5 ], result.mapto)

class TestGridBins(unittest.TestCase):

    def test_init(self):

        a = GridFieldAxis(zeroth=-0.5, step=1.0, count=8, var_name='a', circular=False)
        b = GridFieldAxis(zeroth=-0.5, step=1.0, count=8, var_name='b', circular=False)
        
        # check normal constructor works
        grid = GridBins([a,b])
        self.assertEqual('a', grid.axes[0].var_name)
        self.assertEqual('b', grid.axes[1].var_name)

        # check a wrong-sized array is rejected
        self.assertRaises(TypeError, GridBins, 1)
        self.assertRaises(ValueError, GridBins, [ a ])
        self.assertRaises(ValueError, GridBins, [ a, a, a ])

    def test_compute_input_map(self):

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
        
        # compute
        result = grid.compute_input_map(testcube)

        # check
        self.assertTrue(grid == result.outputgrid)
        self.assertTrue(testcube == result.inputcube)
        self.assertEqual((5,), result.mapfrom.shape)
        self.assertEqual(numpy.int32, result.mapfrom.dtype)
        self.assertEqual([0, 1, 2, 3, 4], result.mapfrom.tolist())
        self.assertEqual((5,), result.mapto.shape)
        self.assertEqual(numpy.int32, result.mapto.dtype)
        self.assertEqual([4, 3, 9, 9, 5], result.mapto.tolist())

    def test_get_or_create_field(self):

        # Axes to use for output grid
        bin_latitude = GridFieldAxis(zeroth=-50.0, step=5.0, count=2, var_name='latitude', circular='false')
        bin_longitude = GridFieldAxis(zeroth=-6.0, step=3.0, count=3, var_name='longitude', circular='false')
        
        # Example to bin onto
        testgrid = GridBins([bin_latitude, bin_longitude])

        # Example of existing data
        somefieldhere = numpy.ma.masked_array(
            data = [ [ 1, 2, 3 ], [ 4, 5, 6 ] ],
            mask = [ [ False, False, False], [ False, True, False ] ],
            dtype=numpy.float32)
        
        # Simulate one existing field
        testfield = Cube(somefieldhere, 
                        var_name='somefieldhere', 
                        dim_coords_and_dims=[(bin_latitude.axis_coordinates(),0),(bin_longitude.axis_coordinates(),1)],
                        standard_name='surface_temperature')
        testgrid.append(testfield)

        # Should have just this one field
        self.assertEqual(1, len(testgrid))

        # Retrieve it
        f = testgrid.get_or_create_field(GridFieldDescriptor('somefieldhere'), 999)
        
        # Should not have created anything new
        self.assertEqual(1, len(testgrid))

        # Retrieved field should be the same as the one we made
        self.assertEqual('surface_temperature', f.standard_name)
        numpy.testing.assert_equal([ [ 1, 2, 3 ], [ 4, 5, 6 ] ], f.data.data)
        numpy.testing.assert_equal([ [ False, False, False ], [ False, True, False ] ], f.data.mask)

        # Now describe something we don't yet have
        desc = GridFieldDescriptor(
            var_name='newfield',
            dtype=numpy.float64,
            standard_name='air_temperature',
            long_name="Some new field here (K)",
            units='K')

        # And retrieve it
        g = testgrid.get_or_create_field(desc, 876)

        # This time should have a new field
        self.assertEqual(2, len(testgrid))
        
        # And initialised to zeros with everything masked out
        numpy.testing.assert_equal([ [ 876, 876, 876 ], [ 876, 876, 876 ] ], g.data.data)
        numpy.testing.assert_equal([ [ True, True, True ], [ True, True, True ] ], g.data.mask)
        
        # re-requesting either should change nothing (and initial value is ignored)
        f = testgrid.get_or_create_field(GridFieldDescriptor('somefieldhere'), 23)
        f = testgrid.get_or_create_field(GridFieldDescriptor('somefieldhere'), 57)
        g = testgrid.get_or_create_field(desc, 82)
        g = testgrid.get_or_create_field(desc, 26)
        g = testgrid.get_or_create_field(desc, 10937)
        self.assertEqual(2, len(testgrid))
        numpy.testing.assert_equal([ [ 1, 2, 3 ], [ 4, 5, 6 ] ], f.data.data)
        numpy.testing.assert_equal([ [ False, False, False ], [ False, True, False ] ], f.data.mask)
        numpy.testing.assert_equal([ [ 876, 876, 876 ], [ 876, 876, 876 ] ], g.data.data)
        numpy.testing.assert_equal([ [ True, True, True ], [ True, True, True ] ], g.data.mask)
