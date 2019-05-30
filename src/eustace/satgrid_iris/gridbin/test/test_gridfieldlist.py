"""Tests for gridfieldlist."""

import unittest
import numpy
from iris.cube import Cube
from iris.coords import AuxCoord
from iris.coords import DimCoord
from ..gridfieldlist import GridFieldAxis
from ..gridfieldlist import GridFieldDescriptor
from ..gridfieldlist import GridFieldList

class TestGridFieldAxis(unittest.TestCase):

    def test_init(self):
        a = GridFieldAxis(zeroth=2.0, step=0.3, count=10, var_name='Splendid', circular=True)
        self.assertEqual(2.0, a.zeroth)
        self.assertEqual(0.3, a.step)
        self.assertEqual(10, a.count)
        self.assertEqual('Splendid', a.var_name)
        self.assertEqual(True, a.circular)

    def test_compute_indices(self):

        # test data
        testvalues = numpy.array([0.25, 0.65, 0.8, 1.6, 1.65], numpy.float32)

        # regrid onto grid boxes starting at -0.5, 0.0, 0.5, 1.0, 1.5
        axis = GridFieldAxis(zeroth=-0.75, step=0.5, count=5, var_name='test', circular=False)
        result = axis.compute_indices(testvalues)

        # check
        self.assertEqual([1, 2, 2, 4, 4], result.tolist())

    def test_compute_indices_error_check(self):

        # check for underflow (this would regrid onto negative indices)
        self.assertRaises(ValueError, GridFieldAxis.compute_indices,
                          GridFieldAxis(zeroth=-1.75, step=0.5, count=1000, var_name='test', circular=False),
                          numpy.array([-3.0, -1.0, 2.1, 5.9], numpy.float32))

        # check for overflow (this would regrid outside the axis length)
        self.assertRaises(ValueError, GridFieldAxis.compute_indices,
                          GridFieldAxis(zeroth=0.5, step=1.0, count=3, var_name='test', circular=False),
                          numpy.array([2.1, 5.9], numpy.float32))

    def test_compute_indices_circular(self):

        testvalues = numpy.array([0.0, -3.2, 4.0, 2.2, 3.99], numpy.float32)

        # without wrap this should have an error due to 4.0 being on the boundary
        self.assertRaises(ValueError, GridFieldAxis.compute_indices,
                          GridFieldAxis(zeroth=-4.5, step=1.0, count=8, var_name='test', circular=False), testvalues)

        # but with wrap we should get a result with 4.0 mapped to grid zero
        result = GridFieldAxis(zeroth=-4.5, step=1.0, count=8, var_name='test', circular=True).compute_indices(testvalues)
        numpy.testing.assert_equal(result, [4, 0, 0, 6, 7])

class TestGridFieldDescriptor(unittest.TestCase):

    def test_init(self):
        d = GridFieldDescriptor(
            var_name='Bob',
            dtype=numpy.float32,
            standard_name='bobs_stuff',
            long_name="The stuff of Bob (kg)",
            units='kg')
        self.assertEqual('Bob', d.var_name)
        self.assertEqual(numpy.float32, d.dtype)
        self.assertEqual('bobs_stuff', d.standard_name)
        self.assertEqual('The stuff of Bob (kg)', d.long_name)
        self.assertEqual('kg', d.units)


class TestGridFieldList(unittest.TestCase):

    def test_init(self):

        a = GridFieldAxis(zeroth=-0.5, step=1.0, count=8, var_name='a', circular=False)
        b = GridFieldAxis(zeroth=-0.5, step=1.0, count=8, var_name='b', circular=False)
        
        # check normal constructor works
        grid = GridFieldList([a,b])
        self.assertEqual('a', grid.axes[0].var_name)
        self.assertEqual('b', grid.axes[1].var_name)

        # check a wrong-sized array is rejected
        self.assertRaises(TypeError, GridFieldList, 1)

    def test_get_or_create_field(self):

        # Axes to use for output grid
        bin_latitude = GridFieldAxis(zeroth=-50.0, step=5.0, count=2, var_name='latitude', circular='false')
        bin_longitude = GridFieldAxis(zeroth=-6.0, step=3.0, count=3, var_name='longitude', circular='false')
        
        # Example to bin onto
        testgrid = GridFieldList([bin_latitude, bin_longitude])

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
