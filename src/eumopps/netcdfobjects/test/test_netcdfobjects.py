"""Tests for Catalogue input/output in NetCDF format."""

# pylint: disable=missing-docstring,invalid-name

import unittest
from netCDF4 import Dataset
from netCDF4 import chartostring
from netCDF4 import stringtoarr
from netCDF4 import num2date
import tempfile
import shutil
from eumopps.netcdfobjects import netcdfobjects
from datetime import datetime
import numpy
import os
import subprocess

class TestObject(object):
    def __init__(self, **kwargs):
        for key, value in kwargs.iteritems():
            setattr(self, key, value)

class TestNetCDFObjects(unittest.TestCase):
    """Tests for reading and writing NetCDF objects.
       These are tested further in catalogue integration tests."""
  
    def test_is_list_of_simple_python_objects(self):

        simple=[ TestObject(a='Hello'), TestObject(a='Bob') ]
        nested = [ TestObject(a=TestObject(a='Hello')), TestObject(a=TestObject(a='Bob')) ]
        self.assertTrue(netcdfobjects.is_list_of_simple_python_objects(simple))
        self.assertFalse(netcdfobjects.is_list_of_simple_python_objects(nested))
        
    def test_write_empty(self):

        # temporary file to use for output
        testfile = tempfile.NamedTemporaryFile(suffix='.nc')

        # instantiate writer and do output
        nc = netcdfobjects.nc_open_for_save(testfile.name)
        netcdfobjects.save_object(nc, 'emptytestobject', None)
        nc.close()

        # check we can load it with netCDF library
        netcdf = Dataset(testfile.name, 'r')   # pylint: disable=unused-variable

    def test_load_object_with_empty_list(self):
        """Read a file that has one object and an empty list."""

        # Empty file with unique temporary name
        testfile = tempfile.NamedTemporaryFile(suffix='.nc')

        # Build NetCDF manually
        nc_testdata = Dataset(testfile.name, 'w')
        nc_testdata.setncattr('Conventions', 'CF-1.6')
        nc_testdata.setncattr('python_class', 'eumopps.netcdfobjects.test.test_netcdfobjects.TestObject')
        group_a = nc_testdata.createGroup('a')
        group_a.createDimension('list_count', 0)
        nc_testdata.close()

        # Read
        nc_input = netcdfobjects.nc_open_for_load(testfile.name)
        stuff = netcdfobjects.load_object(nc_input)
        nc_input.close()

        # Check result (should be empty dataset list and have the specified
        # identifier)
        self.assertIsInstance(stuff, TestObject)
        self.assertTrue(isinstance(stuff.a, list))
        self.assertEqual(0, len(stuff.a))

    def test_save_simple_python_dict(self):

        testfile = tempfile.NamedTemporaryFile(suffix='.nc')
        nc_testdata = netcdfobjects.nc_open_for_save(testfile.name)
        testdata = { 
            'somename' : 'Bob',
            'anumber' : 529,
            'anothernumber' : -82.444, 
            'namelist' : [ 'some name', 'another name', 'something else' ],
            'oneday' : datetime(1922, 3, 3),
            'oneday_in_array' : [ datetime(1922, 3, 3) ],
            'datelist' : [ datetime(2017, 10, 10), datetime(1853, 6, 7) ],
            'myintegerlist' : [ -39, 3, 888, 19827 ],
            'somemorenumbers' : [ 89.999, -1.234, 3.33333 ],
            }
        netcdfobjects.save_simple_python_object( nc_testdata, "mytest", testdata)
        nc_testdata.close()
    

        nc_input = netcdfobjects.nc_open_for_load(testfile.name)

        self.assertEqual('CF-1.6', nc_input.Conventions)

        mytest = nc_input.groups['mytest']

        self.assertEqual('Bob', mytest.somename)
        self.assertEqual(529, mytest.anumber)
        self.assertEqual(-82.444, mytest.anothernumber)
        names = [ str(item) for item in chartostring(mytest.variables['namelist'][:]) ]
        self.assertEqual(3, len(names))
        self.assertEqual([ 'some name', 'another name', 'something else' ], names)
        self.assertEqual(datetime(1922, 3, 3), num2date(mytest.variables['oneday'][:], mytest.variables['oneday'].units))
        self.assertEqual([ datetime(1922, 3, 3) ],  num2date(mytest.variables['oneday_in_array'][:], mytest.variables['oneday_in_array'].units))
        datelist = num2date(mytest.variables['datelist'][:], mytest.variables['datelist'].units)
        self.assertEqual(2, len(datelist))
        self.assertEqual(datetime(2017, 10, 10), datelist[0])
        self.assertEqual(datetime(1853, 6, 7), datelist[1])
        numpy.testing.assert_equal([ -39, 3, 888, 19827 ], mytest.variables['myintegerlist'][:])
        numpy.testing.assert_equal([  89.999, -1.234, 3.33333 ], mytest.variables['somemorenumbers'][:])

    def test_save_simple_python_object(self):

        testfile = tempfile.NamedTemporaryFile(suffix='.nc')
        nc_testdata = netcdfobjects.nc_open_for_save(testfile.name)
        testdata = TestObject(
            somename='Bob',
            anumber=529,
            anothernumber=-82.444, 
            namelist=[ 'some name', 'another name', 'something else' ],
            oneday=datetime(1922, 3, 3),
            oneday_in_array=[ datetime(1922, 3, 3) ],
            datelist=[ datetime(2017, 10, 10), datetime(1853, 6, 7) ],
            myintegerlist=[ -39, 3, 888, 19827 ],
            somemorenumbers=[ 89.999, -1.234, 3.33333 ]
            )
        netcdfobjects.save_simple_python_object( nc_testdata, "mytest", testdata)
        nc_testdata.close()
    

        nc_input = netcdfobjects.nc_open_for_load(testfile.name)
        self.assertEqual('CF-1.6', nc_input.Conventions)
        mytest = nc_input.groups['mytest']
        self.assertEqual('eumopps.netcdfobjects.test.test_netcdfobjects.TestObject', mytest.python_class)
        self.assertEqual('Bob', mytest.somename)
        self.assertEqual(529, mytest.anumber)
        self.assertEqual(-82.444, mytest.anothernumber)
        names = [ str(item) for item in chartostring(mytest.variables['namelist'][:]) ]
        self.assertEqual(3, len(names))
        self.assertEqual([ 'some name', 'another name', 'something else' ], names)
        self.assertEqual(datetime(1922, 3, 3), num2date(mytest.variables['oneday'][:], mytest.variables['oneday'].units))
        self.assertEqual([ datetime(1922, 3, 3) ],  num2date(mytest.variables['oneday_in_array'][:], mytest.variables['oneday_in_array'].units))
        datelist = num2date(mytest.variables['datelist'][:], mytest.variables['datelist'].units)
        self.assertEqual(2, len(datelist))
        self.assertEqual(datetime(2017, 10, 10), datelist[0])
        self.assertEqual(datetime(1853, 6, 7), datelist[1])
        numpy.testing.assert_equal([ -39, 3, 888, 19827 ], mytest.variables['myintegerlist'][:])
        numpy.testing.assert_equal([  89.999, -1.234, 3.33333 ], mytest.variables['somemorenumbers'][:])


    def test_save_load_simple_python_object(self):

        # Make example data as before
        testfile = tempfile.NamedTemporaryFile(suffix='.nc')
        nc_testdata = netcdfobjects.nc_open_for_save(testfile.name)
        testdata = TestObject(
            somename='Bob',
            anumber=529,
            anothernumber=-82.444, 
            namelist=[ 'some name', 'another name', 'something else' ],
            oneday=datetime(1922, 3, 3),
            oneday_in_array=[ datetime(1922, 3, 3) ],
            datelist=[ datetime(2017, 10, 10), datetime(1853, 6, 7) ],
            myintegerlist=[ -39, 3, 888, 19827 ],
            somemorenumbers=[ 89.999, -1.234, 3.33333 ]
            )
        netcdfobjects.save_simple_python_object( nc_testdata, "mytest", testdata)
        nc_testdata.close()

        # Delete just to be sure we really load new things
        del testdata
        del nc_testdata

        # Now load it (should construct object)
        nc_input = netcdfobjects.nc_open_for_load(testfile.name)
        result = netcdfobjects.load_object(nc_input.groups['mytest'])
        nc_input.close()

        # Should be instance of original
        self.assertIsInstance(result, TestObject)

        # And values should be the same
        self.assertEqual('Bob', result.somename)
        self.assertEqual(529, result.anumber)
        self.assertEqual(-82.444, result.anothernumber)
        self.assertEqual([ 'some name', 'another name', 'something else' ], result.namelist)
        self.assertEqual(datetime(1922, 3, 3), result.oneday)
        self.assertEqual([ datetime(1922, 3, 3) ],  result.oneday_in_array)
        self.assertEqual(2, len(result.datelist))
        self.assertEqual(datetime(2017, 10, 10), result.datelist[0])
        self.assertEqual(datetime(1853, 6, 7), result.datelist[1])
        numpy.testing.assert_equal([ -39, 3, 888, 19827 ], result.myintegerlist)
        numpy.testing.assert_equal([  89.999, -1.234, 3.33333 ], result.somemorenumbers)
