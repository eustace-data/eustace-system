"""Tests for definition to ensure full range of EUSTACE temperatures can be represented."""

import unittest
import tempfile
import numpy
from eustace.outputformats import definitions
from eustace.outputformats.globalfield_filebuilder import DatasetAttributesGlobalField
from eustace.timeutils.epoch import days_since_epoch
from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eustace.outputformats.outputvariable import OutputVariable
import netCDF4
import time
from datetime import datetime

class TestFileBuilder(unittest.TestCase):

    def test_all_methods(self):

        # Test data
        testfile = tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.test_filebuilder.', suffix='.nc')

    # Build attributes for example
        attributes = DatasetAttributesGlobalField(
            dataset='Example',
            version='A',
            mainvariable='tas',
            source='B',
            institution='MO',
            comment='EUSTACE project example file format for global field',
            history='Created ' + time.strftime('%c'))

# Day number for the given date
        daynumber = int(days_since_epoch(datetime(2015, 11, 5)))

# object to build global field file at current time
        builder = FileBuilderGlobalField(testfile.name, daynumber, **attributes.__dict__)

#fill field with -95 degree C temperatures which is representative of the lowest temperature likely
#to be found in  EUSTACE
        shape = definitions.GLOBAL_FIELD_SHAPE

        testdata_values = numpy.full(shape,-95.+273.15)
        testdata_mask   = numpy.full(shape, False)

        testdata = numpy.ma.masked_array(data=testdata_values, mask=testdata_mask)

        builder.add_global_field(definitions.TAS, testdata)

        builder.save_and_close()

        # Check results
        result = netCDF4.Dataset(testfile.name, 'r')

        #check that the data haven't wrapped
        numpy.testing.assert_almost_equal(result.variables['tas'][:], testdata, decimal=4)

        #check that the time offsets match the longitudes
        numpy.testing.assert_almost_equal(result.variables['longitude'][:]/360., result.variables['timeoffset'][:], decimal=6)
