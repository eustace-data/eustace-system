import unittest
import numpy
import uuid
from tempfile import NamedTemporaryFile
import struct

from ..locationpointers import LocationPointersRawBinaryWriter
from ..observationsource_rawbinary import ObservationRawBinaryWriter
from ..observationsource_rawbinary import LocalCorrelationRangeRawBinaryWriter
from ..observationsource_rawbinary import LocationLookupRawBinaryWriter
from ..observationsource_rawbinary import LocationLookupWithID
from ...observationsource import Observations

class TestLocationPointerRawBinaryWriter(unittest.TestCase):

    def test_init(self):

        # Test data: correlation ranges
        testcorfile = NamedTemporaryFile(prefix='TestLocationPointerRawBinaryWriter_cor_', suffix='.bin')
        LocalCorrelationRangeRawBinaryWriter().write(testcorfile.name, [ 0.2345 ])

        # Test data: location lookup
        testlocfile = NamedTemporaryFile(prefix='TestLocationPointerRawBinaryWriter_loc_', suffix='.bin')
        testlocations = numpy.array( [ [   52.0,  -9.0, 22.0 ],
                                       [ -179.0, 133.0,  8.2 ] ], numpy.float64)
        testuuid =  uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6')
        LocationLookupRawBinaryWriter().write(testlocfile.name, LocationLookupWithID(testuuid, testlocations))

        # instantiate
        writer = LocationPointersRawBinaryWriter(pathname_locationlookup=testlocfile.name, pathname_correlationranges=testcorfile.name)

        # check values
        self.assertEqual(uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), writer.location_uuid)
        self.assertEqual(1, writer.numcorrelationranges)
        self.assertEqual(3, len(writer.pointers))
        self.assertTrue(isinstance(writer.pointers[0], list))
        self.assertTrue(isinstance(writer.pointers[1], list))
        self.assertTrue(isinstance( writer.pointers[2], list))
        self.assertEqual(0, len(writer.pointers[0]))
        self.assertEqual(0, len(writer.pointers[1]))
        self.assertEqual(0, len(writer.pointers[2]))

    def test_append_days_and_write(self):

        # Test data: UUID
        testuuid = uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6')

        # Test data: location lookup
        testlocfile = NamedTemporaryFile(prefix='TestLocationPointerRawBinaryWriter_loc_', suffix='.bin')
        testlocations = numpy.array( [ [   52.0,  -9.0, 22.0 ],
                                       [ -179.0, 133.0,  8.2 ] ], numpy.float64)
        LocationLookupRawBinaryWriter().write(testlocfile.name, LocationLookupWithID(testuuid, testlocations))

        # Test data: correlation ranges
        testcorfile = NamedTemporaryFile(prefix='TestLocationPointerRawBinaryWriter_cor_', suffix='.bin')
        LocalCorrelationRangeRawBinaryWriter().write(testcorfile.name, [ 0.2345 ])

        # Invent two days of data
        testdataA = Observations(mask=numpy.array([False, False, False], numpy.bool),
                                time=numpy.float32(38888),
                                location=numpy.array([0,1,2], numpy.int32), 
                                measurement=numpy.array([23.2,88.8,999.0], numpy.float64), 
                                uncorrelatederror=numpy.array([3.5,2.5,1.111], numpy.float64),
                                locallycorrelatederror=numpy.array([ [ 2.3, 3.4, 5.6 ] ], numpy.float64))

        testdataB = Observations(mask=numpy.array([False, False], numpy.bool),
                                time=numpy.float32(40981),
                                location=numpy.array([2,1], numpy.int32), 
                                measurement=numpy.array([54.2, 67.22], numpy.float64), 
                                uncorrelatederror=numpy.array([1.111,2.222], numpy.float64),
                                locallycorrelatederror=numpy.array([ [ 9.0, 7.8 ] ], numpy.float64))

        # Write
        testfileA = NamedTemporaryFile(prefix='TestLocationPointerRawBinaryWriter_obs_', suffix='.bin')
        testfileB = NamedTemporaryFile(prefix='TestLocationPointerRawBinaryWriter_obs_', suffix='.bin')
        ObservationRawBinaryWriter().write_day(testfileA.name, testuuid, testdataA, 38888)
        ObservationRawBinaryWriter().write_day(testfileB.name, testuuid, testdataB, 40981)

        # instantiate
        writer = LocationPointersRawBinaryWriter(pathname_locationlookup=testlocfile.name, pathname_correlationranges=testcorfile.name)

        # append days
        writer.append_day(testfileA.name, 38888)
        writer.append_day(testfileB.name, 40981)

        # header size is 48 bytes and each observation record is:
        #   * location: 8 bytes
        #   * measurement: 8 bytes
        #   * uncorrelated error: 8 bytes
        #   * locally correlated error: 8 bytes
        #   ---------
        #    32 bytes
        #
        # offsets add 8 bytes to this to skip location id
        #

        self.assertEqual(3, len(writer.pointers))
        self.assertEqual([ (38888,  56) ], writer.pointers[0])
        self.assertEqual([ (38888,  88), (40981, 88) ], writer.pointers[1])
        self.assertEqual([ (38888, 120), (40981, 56) ], writer.pointers[2])

        # Write output
        testoutput = NamedTemporaryFile(prefix='TestLocationPointerRawBinaryWriter_pointers_', suffix='.bin')
        writer.write(testoutput.name)

        # Read output contents
        file_handle = open(testoutput.name, 'rb')
        formatid = str(file_handle.read(16).decode('utf-8'))
        codeid = str(file_handle.read(16).decode('utf-8'))
        locuuid = uuid.UUID(bytes=file_handle.read(16))
        loccount =  struct.unpack('q', file_handle.read(8))[0]
        table = struct.unpack('qqqqqq', file_handle.read(6*8))
        expected_offset0 = file_handle.tell()
        records0 = struct.unpack('qq', file_handle.read(2*8))
        expected_offset1 = file_handle.tell()
        records1 = struct.unpack('qqqq', file_handle.read(4*8))
        expected_offset2 = file_handle.tell()
        records2 = struct.unpack('qqqq', file_handle.read(4*8))
        file_handle.close()

        # Check values
        self.assertEqual('EUSTACELOCP00001', formatid)
        self.assertEqual('NOCODEVERSION   ', codeid)
        self.assertEqual(uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), locuuid)
        self.assertEqual(3, loccount)
        self.assertEqual((expected_offset0, 1, expected_offset1, 2, expected_offset2, 2), table)
        self.assertEqual((38888,  56), records0)
        self.assertEqual((38888,  88, 40981, 88), records1)
        self.assertEqual((38888, 120, 40981, 56), records2)
        
        # Check that reading an observation record also works and corresponds to the offsets
        file_handle = open(testfileA.name, 'rb')
        file_handle.seek(120)
        obscheck = struct.unpack('ddd', file_handle.read(24))
        file_handle.close()
        self.assertAlmostEqual((999.0, 1.111, 5.6), obscheck)
