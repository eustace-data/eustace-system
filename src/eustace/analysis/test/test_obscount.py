"""Test the obscount class."""

import unittest
from datetime import datetime
import os
import numpy
import uuid
from ..obscount import ObsCounter
from ..observationsource import Observations
from ..fileio.observationsource_rawbinary import ObservationRawBinaryWriter
from ..fileio.observationsource_rawbinary import LocalCorrelationRangeRawBinaryWriter

from tempfile import NamedTemporaryFile
# There happens to be a useful temp directory manager in this module
# - could be placed somewhere more generic perhaps
from eumopps.catalogue.fileio.test.test_formatnetcdf import NamedTemporaryDirectory

class TestObsCounter(unittest.TestCase):

    def test_init(self):

        counter = ObsCounter(path='/some/where', source='myinfo', observable='what', startdate='18500104', enddate='18510105')
        self.assertEqual('/some/where', counter.path)
        self.assertEqual('myinfo', counter.source)
        self.assertEqual('what', counter.observable)
        self.assertEqual(3, counter.startday) # these use the 01/01/1850 epoch
        self.assertEqual(369, counter.endday)

    def test_count_singleday(self):

        # Make directory
        tempdir = NamedTemporaryDirectory()

        # Invent some data
        testdata = Observations(mask=numpy.array([False, False, False], numpy.bool),
                                time=numpy.float32(5),
                                location=numpy.array([78,22,5999999], numpy.int32), 
                                measurement=numpy.array([23.2,88.8,999.0], numpy.float64), 
                                uncorrelatederror=numpy.array([3.5,2.5,1.111], numpy.float64),
                                locallycorrelatederror=numpy.array([ [ 2.3, 3.4, 5.6 ], [ 9.9, 8.8, 7.7 ] ], numpy.float64))
        
        # Make writer instance and store test data
        ObservationRawBinaryWriter().write_day(os.path.join(*[ tempdir.pathname, 'mysurface', '1850', 'mysurface_measurements_18500106.bin' ]),
                                               uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), testdata, 5)

        # Instantiate counter
        counter = ObsCounter(tempdir.pathname, source='mysurface', observable='measurements', startdate='18500101', enddate='20170602')

        # Check this one file
        self.assertEqual(3, counter.count_single_day(num_local_correlation_ranges=2, daynumber=5))

        # Any other day should return None
        self.assertIsNone(counter.count_single_day(num_local_correlation_ranges=2, daynumber=12345))


    def test_read_num_local_correlation_ranges(self):

        # Make directory
        tempdir = NamedTemporaryDirectory()

        # Make test data
        LocalCorrelationRangeRawBinaryWriter().write(
            os.path.join(tempdir.pathname, 'mysurface_measurements_localcorrelationranges.bin'),
            [ 0.99, 0.88, 0.77, 0.66, 0.55 ])

        # Instantiate counter
        counter = ObsCounter(tempdir.pathname, source='mysurface', observable='measurements', startdate='18500101', enddate='20170602')

        # Check ranges
        self.assertEqual(5, counter.read_num_local_correlation_ranges())


