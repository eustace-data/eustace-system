"""Test for file I/O as raw binary files."""

import unittest
import numpy
import os
import struct
import uuid
from tempfile import NamedTemporaryFile
from ...observationsource import ObservationSource
from ...observationsource import Observations
from ..observationsource_rawbinary import ObservationSourceBinaryFilenameGenerator
from ..observationsource_rawbinary import ObservationRawBinaryWriter
from ..observationsource_rawbinary import ObservationRawBinaryReader
from ..observationsource_rawbinary import LocationLookupWithID
from ..observationsource_rawbinary import LocationLookupRawBinaryWriter
from ..observationsource_rawbinary import LocationLookupRawBinaryReader
from ..observationsource_rawbinary import LocalCorrelationRangeRawBinaryWriter
from ..observationsource_rawbinary import LocalCorrelationRangeRawBinaryReader
from ...test.test_diagnosticgrid import SimulatedObservationSource
from eustace.satgrid_iris.gridbin.gridfieldlist import GridFieldAxis
from eustace.outputformats import ensuredirectory

# There happens to be a useful temp directory manager in this module
# - could be placed somewhere more generic perhaps
from eumopps.catalogue.fileio.test.test_formatnetcdf import NamedTemporaryDirectory

class TestObservationSourceRawBinaryFilenameGenerator(unittest.TestCase):
    """Test the utility for filename generation."""

    def test_init(self):
        result = ObservationSourceBinaryFilenameGenerator(sourcename='garden', basepath='/bob/gate')
        self.assertEqual(result.basepath, '/bob/gate')
        self.assertEqual(result.sourcename, 'garden')

    def test_filename_observations(self):
        generator = ObservationSourceBinaryFilenameGenerator(sourcename='gate', basepath='/garden/path')
        result = generator.filename_observations(observable='geraniums', daynumber=57500)
        self.assertEqual('/garden/path/gate/2007/gate_geraniums_20070607.bin', result)

    def test_filename_patterns_observations(self):
        generator = ObservationSourceBinaryFilenameGenerator(sourcename='gate', basepath='/garden/path')
        result = generator.filename_patterns_observations(observable='geraniums')
        self.assertEqual(['/garden/path/gate/%Y', 'gate_geraniums_%Y%m%d.bin' ], result)

    def test_filename_local_correlation_ranges(self):
        generator = ObservationSourceBinaryFilenameGenerator(sourcename='gate', basepath='/garden/path')
        result = generator.filename_local_correlation_ranges(observable='daffodils')
        self.assertEqual('/garden/path/gate_daffodils_localcorrelationranges.bin', result)

    def test_filename_patterns_mobile_locations(self):
        generator = ObservationSourceBinaryFilenameGenerator(sourcename='gate', basepath='/garden/path')
        result = generator.filename_patterns_mobile_locations(observable='geraniums')
        self.assertEqual(['/garden/path/gate/%Y', 'gate_geraniums_mobilelocations_%Y%m%d.bin' ], result)

    def test_filename_mobile_locations(self):
        generator = ObservationSourceBinaryFilenameGenerator(sourcename='gate', basepath='/garden/path')
        result = generator.filename_mobile_locations(observable='geraniums', daynumber=57500)
        self.assertEqual('/garden/path/gate/2007/gate_geraniums_mobilelocations_20070607.bin', result)


class TestObservationRawBinaryWriter(unittest.TestCase):
    """Test storage."""

    def test_init(self):

        # make item
        generator = ObservationRawBinaryWriter()

    def test_write_day(self):

        # Make temp directory (gets deleted automatically when out of scope)
        tempdirectory = NamedTemporaryDirectory()
        # print 'Saving to: ', tempdirectory.pathname
        
        # Make writer instance
        generator = ObservationSourceBinaryFilenameGenerator(sourcename='playing', basepath=tempdirectory.pathname)
        writer =  ObservationRawBinaryWriter()

        # Simulated data
        source = SimulatedObservationSource()
        obs = source.observations('pretend')

        # Store two days (should only have non-empty for the first)
        writer.write_day(generator.filename_observations('pretend', 57500), uuid.UUID('4996797e-f0dc-4495-b91d-4345189fcb71'), obs, 57500)
        writer.write_day(generator.filename_observations('pretend', 57501), uuid.UUID('cc539a51-7f78-4fb8-9caa-d57229994a7c'), obs, 57501)

        # Should have made files for each day (even though second one has no data)
        fileA = open(os.path.join(*[tempdirectory.pathname, 'playing', '2007', 'playing_pretend_20070607.bin']), 'rb')
        fileB = open(os.path.join(*[tempdirectory.pathname, 'playing', '2007', 'playing_pretend_20070608.bin']), 'rb')

        # Check headers
        formatidA = str(fileA.read(16).decode('utf-8'))
        formatidB = str(fileB.read(16).decode('utf-8'))
        codeA = str(fileA.read(16).decode('utf-8'))
        codeB = str(fileB.read(16).decode('utf-8'))
        locuuidA = uuid.UUID(bytes=fileA.read(16))
        locuuidB = uuid.UUID(bytes=fileB.read(16))
        self.assertEqual('EUSTACEOBSN00001', formatidA)
        self.assertEqual('EUSTACEOBSN00001', formatidB)
        self.assertEqual('NOCODEVERSION   ', codeA)
        self.assertEqual('NOCODEVERSION   ', codeB)
        self.assertEqual(uuid.UUID('4996797e-f0dc-4495-b91d-4345189fcb71'), locuuidA)
        self.assertEqual(uuid.UUID('cc539a51-7f78-4fb8-9caa-d57229994a7c'), locuuidB)

        # Read data
        a = numpy.fromfile(
            fileA,
            dtype=[('location', numpy.uint64), ('measurement', numpy.float64), ('uncorrelatederror', numpy.float64), ('local_0', numpy.float64), ('local_1', numpy.float64)])
        b = numpy.fromfile(
            fileB,
            dtype=[('location', numpy.uint64), ('measurement', numpy.float64), ('uncorrelatederror', numpy.float64), ('local_0', numpy.float64), ('local_1', numpy.float64)])

        # Check sizes as expected
        self.assertEqual((3,), a.shape)
        self.assertEqual((0,), b.shape)

        # Check values
        numpy.testing.assert_equal([0, 3, 4], a['location'])
        numpy.testing.assert_equal([275, 204, 310], a['measurement'])
        numpy.testing.assert_equal([5, 7, 8], a['uncorrelatederror'])
        numpy.testing.assert_almost_equal([2.3, 1.2, 8.9], a['local_0'], decimal=6)
        numpy.testing.assert_almost_equal([0.222, 8.888, 9.999], a['local_1'], decimal=6)

class TestObservationRawBinaryReader(unittest.TestCase):
    """Test retrieval."""

    def test_read(self):
        """Test using a round-trip from storage."""        

        # Test file to write (deleted automatically on exit)
        testfile = NamedTemporaryFile(prefix='TestObservationRawBinaryReader', suffix='.bin')
        
        # Invent some data
        testdata = Observations(mask=numpy.array([False, False, False], numpy.bool),
                                time=numpy.float32(38888),
                                location=numpy.array([78,22,5999999], numpy.int32), 
                                measurement=numpy.array([23.2,88.8,999.0], numpy.float64), 
                                uncorrelatederror=numpy.array([3.5,2.5,1.111], numpy.float64),
                                locallycorrelatederror=numpy.array([ [ 2.3, 3.4, 5.6 ], [ 9.9, 8.8, 7.7 ] ], numpy.float64))
        
        # Make writer instance and store test data
        ObservationRawBinaryWriter().write_day(testfile.name, uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), testdata, 38888)

        # Load
        result = ObservationRawBinaryReader(2).read_day(testfile.name, uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), 38888)

        # Check results
        self.assertEqual(38888, result.time)
        numpy.testing.assert_equal([78,22,5999999], result.location)
        numpy.testing.assert_equal([23.2,88.8,999.0], result.measurement)
        numpy.testing.assert_equal([3.5,2.5,1.111], result.uncorrelatederror)
        numpy.testing.assert_equal([ [ 2.3, 3.4, 5.6 ], [ 9.9, 8.8, 7.7 ] ], result.locallycorrelatederror)

	# Check it correctly handles cases where file name is None
        result = ObservationRawBinaryReader(2).read_day(None, uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), 38888)
        numpy.testing.assert_equal([0], result.location)
        numpy.testing.assert_equal([0.], result.measurement)
        numpy.testing.assert_equal([0.], result.uncorrelatederror)
        
    def test_read_location_ids(self):
        """Test using a round-trip from storage."""        

        # Test file to write (deleted automatically on exit)
        testfile = NamedTemporaryFile(prefix='TestObservationRawBinaryReader', suffix='.bin')
        
        # Invent some data
        testdata = Observations(mask=numpy.array([False, False, False], numpy.bool),
                                time=numpy.float32(38888),
                                location=numpy.array([78,22,5999999], numpy.int32), 
                                measurement=numpy.array([23.2,88.8,999.0], numpy.float64), 
                                uncorrelatederror=numpy.array([3.5,2.5,1.111], numpy.float64),
                                locallycorrelatederror=numpy.array([ [ 2.3, 3.4, 5.6 ], [ 9.9, 8.8, 7.7 ] ], numpy.float64))
        
        # Make writer instance and store test data
        ObservationRawBinaryWriter().write_day(testfile.name, uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), testdata, 38888)

        # Load
        result = ObservationRawBinaryReader(2).read_location_ids(testfile.name, uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'))

        # Check results
        numpy.testing.assert_equal([78,22,5999999], result)

    def test_read_observation_count(self):

        # Test file to write (deleted automatically on exit)
        testfile = NamedTemporaryFile(prefix='TestObservationRawBinaryReader', suffix='.bin')
        
        # Invent some data
        testdata = Observations(mask=numpy.array([False, False, False], numpy.bool),
                                time=numpy.float32(38888),
                                location=numpy.array([78,22,5999999], numpy.int32), 
                                measurement=numpy.array([23.2,88.8,999.0], numpy.float64), 
                                uncorrelatederror=numpy.array([3.5,2.5,1.111], numpy.float64),
                                locallycorrelatederror=numpy.array([ [ 2.3, 3.4, 5.6 ], [ 9.9, 8.8, 7.7 ] ], numpy.float64))
        
        # Make writer instance and store test data
        ObservationRawBinaryWriter().write_day(testfile.name, uuid.UUID('cb456b34-e723-4053-b755-0ab311c58ae6'), testdata, 38888)

        # Load
        result = ObservationRawBinaryReader(2).read_observation_count(testfile.name)

        # Check results
        self.assertEqual(3, result)

class TestLocalCorrelationRangeRawBinaryReader(unittest.TestCase):

    def test_read(self):

        # temporary file of test data
        testfile = NamedTemporaryFile(prefix='TestLocalCorrelationRangeRawBinaryReader', suffix='.bin')

        # make test data by writing binary bytes directly
        simwriter = open(testfile.name, 'wb')
        simwriter.write('EUSTACECORN00001')
        simwriter.write('NOCODEVERSION   ')
        simwriter.write(struct.pack('ddd', 5.55, 2.333, 8.99999))
        simwriter.close()

        # do the read
        result = LocalCorrelationRangeRawBinaryReader().read(testfile.name)
        self.assertEqual((3,), result.shape)
        numpy.testing.assert_equal([ 5.55, 2.333, 8.99999 ],  result)

class TestLocalCorrelationRangeRawBinaryWriter(unittest.TestCase):

    def test_write(self):

        # temporary file of test data
        testfile = NamedTemporaryFile(prefix='TestLocalCorrelationRangeRawBinaryWriter', suffix='.bin')
        
        # make test data
        testdata = [ 0.99, 0.88 ]
        LocalCorrelationRangeRawBinaryWriter().write(testfile.name, testdata)

        # load it
        file_handle = open(testfile.name, 'rb')
        formatid = file_handle.read(16)
        codeid = file_handle.read(16)
        result = struct.unpack('dd', file_handle.read(2*8))
        file_handle.close()

        # check contents
        self.assertEqual('EUSTACECORN00001', formatid)
        self.assertEqual('NOCODEVERSION   ', codeid)
        self.assertEqual(( 0.99, 0.88 ), result)


class TestLocationLookupRawBinaryReader(unittest.TestCase):

    def test_read(self):

        # temporary file of test data
        testfile = NamedTemporaryFile(prefix='TestLocationLookupRawBinaryReader', suffix='.bin')

        # make test data by writing binary bytes directly
        simwriter = open(testfile.name, 'wb')
        simwriter.write('EUSTACELOCN00001')
        simwriter.write('NOCODEVERSION   ')
        simwriter.write(uuid.UUID('6d6c833d-064c-468a-92e6-c20e96675661').bytes)
        simwriter.write(struct.pack('dddddd', 45.0, 0.5, 43.0, 20.0, -20.0, -93.0))
        simwriter.close()

        # do the read
        result = LocationLookupRawBinaryReader().read(testfile.name)
        self.assertEqual(uuid.UUID('6d6c833d-064c-468a-92e6-c20e96675661'), result.uuid)
        self.assertEqual((2,3), result.lookuptable.shape)
        numpy.testing.assert_equal([ [ 45.0, 43.0, -20.0 ],
                                     [  0.5, 20.0, -93.0 ] ],  result.lookuptable)

    def test_read_uuid_and_count(self):

        # temporary file of test data
        testfile = NamedTemporaryFile(prefix='TestLocationLookupRawBinaryReader', suffix='.bin')

        # make test data by writing binary bytes directly
        simwriter = open(testfile.name, 'wb')
        simwriter.write('EUSTACELOCN00001')
        simwriter.write('NOCODEVERSION   ')
        simwriter.write(uuid.UUID('6d6c833d-064c-468a-92e6-c20e96675661').bytes)
        simwriter.write(struct.pack('dddddd', 45.0, 0.5, 43.0, 20.0, -20.0, -93.0))
        simwriter.close()

        # do the read
        location_uuid, count = LocationLookupRawBinaryReader().read_uuid_and_count(testfile.name)
        self.assertEqual(uuid.UUID('6d6c833d-064c-468a-92e6-c20e96675661'), location_uuid)
        self.assertEqual(3, count)

class TestLocationLookupRawBinaryWriter(unittest.TestCase):

    def test_write(self):

        # temporary file of test data
        testfile = NamedTemporaryFile(prefix='TestLocationLookupRawBinaryWriter', suffix='.bin')
        
        # make test data
        testuuid = uuid.UUID('00000000-0000-0000-0000-000000244c88')
        testdata = numpy.array( [ [   52.0,  -9.0, 22.0 ],
                                  [ -179.0, 133.0,  8.2 ] ], numpy.float64)
        LocationLookupRawBinaryWriter().write(testfile.name, LocationLookupWithID(testuuid, testdata))

        # load it
        file_handle = open(testfile.name, 'rb')
        formatid = str(file_handle.read(16).decode('utf-8'))
        codeid = str(file_handle.read(16).decode('utf-8'))
        locuuid = uuid.UUID(bytes=file_handle.read(16))
        table = struct.unpack('dddddd', file_handle.read(6*8))
        file_handle.close()
        self.assertEqual('EUSTACELOCN00001', formatid)
        self.assertEqual('NOCODEVERSION   ', codeid)
        self.assertEqual(uuid.UUID('00000000-0000-0000-0000-000000244c88'), locuuid)
        self.assertEqual(( 52.0, -179.0, -9.0, 133.0, 22.0, 8.2 ), table)
