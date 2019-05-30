"""Tests of file format verification."""

# pylint: disable=missing-docstring,attribute-defined-outside-init,invalid-name

__version__ = "$Revision: 562 $"
__author__ = "Joel R. Mitchelson"

from ..verify_globalfield import VerifyGlobalField
from ..verify_stations import VerifyStationData
import unittest
from netCDF4 import Dataset
import tempfile
import numpy
import os
import StringIO
import sys # uncomment if printed diagnostics are needed

class TestArgs(object):

    def __init__(self):
        pass


class TestVerifyGlobalField(unittest.TestCase):

    def test_verify_global_field_none(self):
        testargs = TestArgs()
        testargs.pathname = None
        self.assertFalse(VerifyGlobalField().run(
            testargs, StringIO.StringIO()))

    def test_verify_global_field_empty(self):
        outputfile = tempfile.NamedTemporaryFile(
            prefix='eustace.test_verify.', suffix='.nc')
        nc = Dataset(outputfile.name, 'w')
        nc.close()
        testargs = TestArgs()
        testargs.pathname = outputfile.name
        self.assertFalse(VerifyGlobalField().run(
            testargs, StringIO.StringIO()))

    def test_verify_global_field_tasmax_only(self):

        # temporary file for testing
        outputfile = tempfile.NamedTemporaryFile(prefix='eustace.test_verify.', suffix='.nc')

        # Create
        nc = Dataset(outputfile.name, 'w')

        # Set header info
        nc.title = 'EUSTACE Surface Air Temperature Estimates Global'
        nc.institution = 'Test Institution Name'
        nc.source = __file__
        nc.history = ''
        nc.Conventions = 'CF-1.6'

        # data to use
        data_lat = numpy.arange(-90.0, 89.75001, 0.25, numpy.float32)
        data_lon = numpy.arange(-180.0, 179.75001, 0.25, numpy.float32)
        data_time = numpy.array([50000.5], numpy.float32)
        data_timebounds = numpy.array([ [50000.0, 50001.0] ], numpy.float32)
        data_timeoffset = numpy.arange(-0.5, 0.499999, (1.0/data_lon.size), numpy.float32)
        data_tasmax = numpy.random.normal(288.0, 50.0, (1, data_lat.size, data_lon.size)).astype(numpy.float32)

        # dimensions
        nc.createDimension('latitude', data_lat.size)
        nc.createDimension('longitude', data_lon.size)
        nc.createDimension('time', None)
        nc.createDimension('bounds', 2)

        # add data
        latitude = nc.createVariable('latitude', numpy.float32, ('latitude',))
        latitude.units = 'degrees_north'
        latitude.standard_name = 'latitude'
        latitude.long_name = 'Latitude (deg)'
        latitude.axis = 'Y'
        latitude[:] = data_lat
        longitude = nc.createVariable('longitude', numpy.float32, ('longitude'))
        longitude.units = 'degrees_east'
        longitude.standard_name = 'longitude'
        longitude.long_name = 'Longitude (deg)'
        longitude.axis = 'X'
        longitude[:] = data_lon
        time = nc.createVariable('time', numpy.float32, ('time',))
        time.units = 'days since 1850-01-01T00:00:00Z'
        time.standard_name = 'time'
        time.long_name = 'Time at zero longitude'
        time.calendar = 'gregorian'
        time.bounds='timebounds'
        time.axis='T'
        time.ancillary_variables='timeoffset'
        time[:] = data_time
        timebounds = nc.createVariable('timebounds', numpy.float32, ('time', 'bounds',))
        timebounds[:] = data_timebounds
        timeoffset = nc.createVariable('timeoffset', numpy.float32, ('longitude',))
        timeoffset.long_name = 'Local time offset from UTC (days)'
        timeoffset.units = 'days'
        timeoffset[:] = data_timeoffset
        tasmax = nc.createVariable('tasmax', numpy.int16, ('time', 'latitude', 'longitude'))
        tasmax.scale_factor = 0.001
        tasmax.add_offset = 273.15
        tasmax.units = 'K'
        tasmax.standard_name = 'air_temperature'
        tasmax.long_name = 'Maximum daily surface air temperature'
        tasmax.cell_methods = 'time: maximum'
        tasmax[:] = data_tasmax

        # save and close
        nc.close()

        # verification params
        testargs = TestArgs()
        testargs.pathname = outputfile.name
        outputstream = StringIO.StringIO()
        # outputstream = sys.stderr # use stderr if diagnostics needed

        # verify - should be ok
        self.assertTrue(VerifyGlobalField().run(testargs, outputstream))


class TestVerifyStationData(unittest.TestCase):

    def test_verify_station_data_none(self):
        testargs = TestArgs()
        testargs.pathnameprefix = None
        self.assertFalse(VerifyStationData().run(testargs, StringIO.StringIO()))

    def test_verify_station_data_empty(self):

        # this placeholder remains empty but we use the name to make temporary
        # temperature and status files
        outputfile = tempfile.NamedTemporaryFile(prefix='eustace_stations_testverify_empty', suffix='.placeholder')

        # generate files using prefix
        pathname_prefix = outputfile.name.split('.placeholder')[0]
        pathname_temperature = pathname_prefix + '_temperature.nc'
        pathname_status = pathname_prefix + '_status.nc'
        temperature_nc = Dataset(pathname_temperature, 'w')
        status_nc = Dataset(pathname_status, 'w')
        temperature_nc.close()
        status_nc.close()

        # args for test
        testargs = TestArgs()
        testargs.pathnameprefix = pathname_prefix

        # Use try-finally to ensure file removed even if test fails
        try:
            self.assertFalse(VerifyStationData().run(testargs, StringIO.StringIO()))
        finally:
            # print 'Unlinking: ', pathname_temperature, ',', pathname_status
            os.unlink(pathname_temperature)
            os.unlink(pathname_status)

    def test_verify_station_data(self):

        # this placeholder remains empty but we use the name to make temporary
        # temperature and status files
        outputfile = tempfile.NamedTemporaryFile(prefix='eustace_stations_testverify_example', suffix='.placeholder')

        # generate files using prefix
        pathname_prefix = outputfile.name.split('.placeholder')[0]
        pathname_temperature = pathname_prefix + '_temperature.nc'
        pathname_status = pathname_prefix + '_status.nc'

        # make files here
        temperature_nc = Dataset(pathname_temperature, 'w')
        status_nc = Dataset(pathname_status, 'w')

        # Set header info (temperature)
        temperature_nc.title = 'EUSTACE Surface Air Temperature Station Temperatures'
        temperature_nc.institution = 'Test Institution Name'
        temperature_nc.source = __file__
        temperature_nc.history = ''
        temperature_nc.Conventions = 'CF-1.6'

        # Set header info (status)
        status_nc.title = 'EUSTACE Surface Air Temperature Station Status'
        status_nc.institution = 'Test Institution Name'
        status_nc.source = __file__
        status_nc.history = ''
        status_nc.Conventions = 'CF-1.6'

        # Time UNLIMITED, 3 stations
        temperature_nc.createDimension('time', None)
        temperature_nc.createDimension('station', 3)
        temperature_nc.createDimension('name_strlen', 32)

        # set time range
        measurement_time = temperature_nc.createVariable('time', numpy.float32, ('time',))
        measurement_time[0:1095] = numpy.arange(59900, 60995, dtype=numpy.float32)
        measurement_time.units = 'days since 1850-01-01T00:00:00Z'
        measurement_time.standard_name = 'time'
        measurement_time.long_name = 'Time at zero longitude'
        measurement_time.calendar = 'gregorian'

        # create variables depending on time
        tasmin = temperature_nc.createVariable(
            'tasmin', numpy.int16, ('time', 'station'))
        tasmax = temperature_nc.createVariable(
            'tasmax', numpy.int16, ('time', 'station'))

        tasmin.units = 'K'
        tasmin.scale_factor = 0.001
        tasmin.add_offset = 273.15
        tasmin.standard_name = 'air_temperature'
        tasmin.long_name = 'Minimum daily surface air temperature'
        tasmin.cell_methods = 'time: minimum'
        tasmin[:] = numpy.random.normal(288.15, 30.0, (1095, 3))

        tasmax.units = 'K'
        tasmax.scale_factor = 0.001
        tasmax.add_offset = 273.15
        tasmax.standard_name = 'air_temperature'
        tasmax.long_name = 'Maximum daily surface air temperature'
        tasmax.cell_methods = 'time: maximum'
        tasmax[:] = numpy.random.normal(288.15, 30.0, (1095, 3))

        status_nc.createDimension('detection_time', None)
        status_nc.createDimension('tasmin_break', None)
        status_nc.createDimension('tasmax_break', None)
        status_nc.createDimension('station', 2)
        status_nc.createDimension('bounds', 2)

        status_time = status_nc.createVariable(
            'detection_time', numpy.float32, ('detection_time',))
        status_time.units = 'days since 1850-01-01T00:00:00Z'
        status_time.standard_name = 'time'
        status_time.long_name = 'Start time of period for break detection status report (days)'
        status_time[0:3] = numpy.array([59900, 60265, 60630], numpy.float32)

        status_nc.close()

        # args for test
        testargs = TestArgs()
        testargs.pathnameprefix = pathname_prefix

        # Use try-finally to ensure file removed even if test fails
        try:
            outputstream = StringIO.StringIO()
            # outputstream = sys.stderr # use stderr if diagnostics needed
            self.assertTrue(VerifyStationData().run(testargs, outputstream))
        finally:
            # print 'Unlinking: ', pathname_temperature, ',', pathname_status
            os.unlink(pathname_temperature)
            os.unlink(pathname_status)
