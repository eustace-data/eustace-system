"""Compute location pointers."""

from observationsource_rawbinary import LocationLookupRawBinaryReader
from observationsource_rawbinary import ObservationRawBinaryReader
from observationsource_rawbinary import ObservationRawBinaryWriter
from observationsource_rawbinary import LocalCorrelationRangeRawBinaryReader
from observationsource_rawbinary import CATEGORY_OBSERVATION
import numpy

FORMATID_LOCATION_POINTERS = 'EUSTACELOCP00001'
"""Format identifier placed at start of correlated uncertainty binary files."""

CODEID =                           'NOCODEVERSION   '
"""Placeholder for code version information should we wish to use that."""

DTYPE = [('daynumber', numpy.uint64), ('fileoffset', numpy.uint64)]

class LocationPointersRawBinaryWriter(object):

    def __init__(self, pathname_locationlookup, pathname_correlationranges):
        """Use specified location lookup just to find total IDs and global UUID."""
        
        # Use location lookup to get UUID and count
        self.location_uuid, location_count =  LocationLookupRawBinaryReader().read_uuid_and_count(pathname_locationlookup)

        # Also count correlation ranges
        self.numcorrelationranges = LocalCorrelationRangeRawBinaryReader().read(pathname_correlationranges).shape[0]

        # Make list of N empty lists
        self.pointers = [ [ ] for index in range(location_count) ]

    def append_day(self, pathname_dailyobservations, daynumber):
        """Append location pointers for the specified daily data."""

        # get location ids
        locationids = ObservationRawBinaryReader(self.numcorrelationranges).read_location_ids(pathname_dailyobservations, self.location_uuid)

        # header size common to all files is 48 and add 8 to skip location id in each observation
        offset_size_bytes = 56

        # record size for this number of correlation ranges
        record_size_bytes = numpy.dtype(ObservationRawBinaryWriter.build_dtype(CATEGORY_OBSERVATION, self.numcorrelationranges)).itemsize

        # map from location id to offset in file
        for index, locationid in enumerate(locationids):
            self.pointers[locationid].append( (daynumber, offset_size_bytes + record_size_bytes*index) )

    def write(self, pathname):
        """Write to file."""

        # open it
        file_handle = open(pathname, 'wb')

        # header including lookup table UUID
        file_handle.write(FORMATID_LOCATION_POINTERS)
        file_handle.write(CODEID)
        file_handle.write(self.location_uuid.bytes)

        # count of location IDs
        location_count = numpy.uint64(len(self.pointers))
        file_handle.write(location_count)

        # the size of the above writes should add up to this
        header_size_bytes = numpy.uint64(file_handle.tell())

        # allocate offset table giving offsets to records for each location
        offset_table = numpy.zeros((location_count, 2), dtype=numpy.uint64)

        # the size of the offset table
        offset_table_size_bytes = offset_table.nbytes

        # cumulative offsets to pointer records
        pointer_record_offset = header_size_bytes + offset_table_size_bytes
        pointer_record_size_bytes = numpy.dtype(DTYPE).itemsize

        if location_count > 0:
        
            # populate offset table
            record_counts = numpy.array([ len(pointer_records) for pointer_records in self.pointers ], dtype=numpy.uint64)
            cumulative_record_counts = numpy.hstack((0, numpy.cumsum(record_counts)[0:-1]))
            record_offsets = (cumulative_record_counts * pointer_record_size_bytes) + pointer_record_offset
            offset_table[:,0] = record_offsets
            offset_table[:,1] = record_counts

            # write offset table
            offset_table.tofile(file_handle)

            # write offsets for each location id
            for pointer_records in self.pointers:
                numpy.array(pointer_records, dtype=DTYPE).tofile(file_handle)

        # finished
        file_handle.close()
        del file_handle
        
class LocationPointersRawBinaryReader(object):

    def __init__(self):
        self.formatid = FORMATID_LOCATION_POINTERS

    def read_pointer_records(self, pathname, verify_location_uuid, location_id):

        # open it
        file_handle = open(pathname, 'rb')

        # read and check format ID
        formatid = file_handle.read(16)
        if formatid != self.formatid:
            raise Exception('Expected format ID not found: {0}'.format(self.formatid))

        # read and ignore version
        codeversionid = file_handle.read(16)

        # read uuid
        location_uuid = uuid.UUID(bytes=file_handle.read(16))

        # check it matches the one given for verification
        if location_uuid != verify_location_uuid:
            raise Exception('Location UUID mismatch')

        # read count
        location_count = numpy.fromfile(file_handle, dtype=numpy.uint64, count=1)

        # check locaton_id in range
        if (location_id < 0) or (location_id >= location_count):
            raise Exception('Location ID {0} is outside of expected range {1}'.format(location_id, location_count))

        # go to offset table entry for this location
        file_handle.seek(16 * location_id, 1)   # The 1 here indicates seek from current position

        # read offset table info
        record_offset = numpy.fromfile(file_handle, dtype=numpy.uint64, count=1)
        record_count = numpy.fromfile(file_handle, dtype=numpy.uint64, count=1)

        # seek to offset
        file_handle.seek(record_offset)

        # read into NumPy structure
        records = numpy.fromfile(file_handle, dtype=DTYPE, count=record_count)

        # done
        return records
