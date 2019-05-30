"""File I/O as raw binary files, suitable for input to analysis processes and diagnostics."""

import numpy
import os
import uuid
from ..observationsource import Observations
from ..observationsource import ObservationSource
from eustace.timeutils.epoch import epoch_plus_days
from eustace.outputformats.ensuredirectory import ensuredirectory
from eumopps.timeutils import datetime_numeric
from eustace.analysis.fileio.observationsource_missing import ObservationsMissing

CATEGORY_OBSERVATION = 'obs'
CATEGORY_LOCATION_LOOKUP = 'loc'
CATEGORY_LOCAL_CORRELATION_RANGE = 'cor'
CATEGORIES = [ CATEGORY_OBSERVATION, CATEGORY_LOCATION_LOOKUP, CATEGORY_LOCAL_CORRELATION_RANGE ]
DTYPE = {
    CATEGORY_OBSERVATION: [('location', numpy.uint64), ('measurement', numpy.float64), ('uncorrelatederror', numpy.float64) ],
    CATEGORY_LOCATION_LOOKUP: [('latitude', numpy.float64), ('longitude', numpy.float64) ],
    CATEGORY_LOCAL_CORRELATION_RANGE: [('range', numpy.float64) ]
}

FORMATID_OBSERVATION =             'EUSTACEOBSN00001'
"""Format identifier placed at start of observation binary files."""

FORMATID_LOCATION_LOOKUP =         'EUSTACELOCN00001'
"""Format identifier placed at start of location binary files."""

FORMATID_LOCAL_CORRELATION_RANGE = 'EUSTACECORN00001'
"""Format identifier placed at start of correlated uncertainty binary files."""

CODEID =                           'NOCODEVERSION   '
"""Placeholder for code version information should we wish to use that."""

class RawBinaryReader(object):
    """Generic reader object used as base class."""

    def __init__(self, formatid):
        """Build this reader."""

        self.formatid = formatid

    def read_header_from_stream(self, file_handle):
        """Read header from open file stream and check format is the one expected."""

        # read and check format ID
        formatid = file_handle.read(16)
        if formatid != self.formatid:
            raise Exception('Expected format ID not found: {0}'.format(self.formatid))

        # read and ignore version
        codeversionid = file_handle.read(16)

    @staticmethod
    def count_remaining(file_handle, record_size_bytes):
        """Count total records between current position in stream and end of stream and seek to end of stream."""

        # measure file pointer at start
        start_offset = file_handle.tell()

        # compute data size by seeking to end
        file_handle.seek(0, 2)
        total_size_bytes = file_handle.tell() - start_offset

        # count is total size divided by record size
        count = total_size_bytes / record_size_bytes

        return count

class RawBinaryLocationCategoryReader(RawBinaryReader):

    def __init__(self, formatid, basecategory):
        """Build reader for files relating to location."""

        super(RawBinaryLocationCategoryReader, self).__init__(formatid)
        self.basecategory = basecategory

    def read_header_from_stream(self, file_handle):
        """Read header including UUID of location lookup."""

        # read header as usual
        super(RawBinaryLocationCategoryReader, self).read_header_from_stream(file_handle)

        # also read and return location uuid at end of header
        location_uuid = uuid.UUID(bytes=file_handle.read(16))
        return location_uuid

    def count_remaining(self, file_handle):
        """Count total records between current position in stream and end of stream and seek to end of stream."""

        # find record size
        record_size_bytes = numpy.dtype(DTYPE[self.basecategory]).itemsize

        # compute remaining records
        return RawBinaryReader.count_remaining(file_handle, record_size_bytes)

    def read_uuid_and_count(self, pathname):
        """Quick retrieval of UUID and total number of location IDs"""

        # open it
        file_handle = open(pathname, 'rb')

        # get uuid
        location_uuid = self.read_header_from_stream(file_handle)

        # get count
        count = self.count_remaining(file_handle)

        # done with file
        file_handle.close()
        del file_handle

        return location_uuid, count

class ObservationRawBinaryReader(RawBinaryLocationCategoryReader):
    """Read observation objects from binary."""

    def __init__(self, num_local_correlation_ranges):
        """Build reader to read days of data on demand."""

        super(ObservationRawBinaryReader, self).__init__(FORMATID_OBSERVATION, CATEGORY_OBSERVATION)
        self.num_local_correlation_ranges = num_local_correlation_ranges

    def read_contents(self, pathname, verify_location_uuid):
        """Read contents into NumPy structure."""
        
        # open file
        file_handle = open(pathname, 'rb')

        # get uuid
        location_uuid = self.read_header_from_stream(file_handle)

        # read data
        contents = numpy.fromfile(file_handle, dtype=ObservationRawBinaryWriter.build_dtype(self.basecategory, self.num_local_correlation_ranges))

	# close it
        file_handle.close()
        del file_handle

        # return NumPy result
        return contents

    def read_day(self, pathname,  verify_location_uuid, daynumber):
        """Read data and populate time member with given parameter.
           Verify location UUID."""

        # read into NumPy structure
        if pathname is not None:
            contents = self.read_contents(pathname, verify_location_uuid)
            # mask (allow everything)
            mask = numpy.zeros((contents.shape[0],), numpy.bool)
        else:
            none_type = ObservationRawBinaryWriter.build_dtype(self.basecategory, self.num_local_correlation_ranges)
            contents = numpy.array([tuple([0.]*len(none_type))], dtype=none_type)
            # mask (allow everything)
            mask = numpy.array([1], numpy.bool)
        
        
        # build locally correlated error list
        locallycorrelatederror = [ ]
        for index in range(self.num_local_correlation_ranges):
            locallycorrelatederror.append( contents[ObservationRawBinaryWriter.locally_correlated_uncertainty_field_name(index)] )

        # Build observations class
        return self.build_observations_instance(mask, daynumber, contents, locallycorrelatederror)

    def build_observations_instance(self, mask, daynumber, contents, locallycorrelatederror):
        """Create Observations object using specified data."""
       
        # return as Observations object
        return Observations(mask, daynumber, contents['location'], contents['measurement'], contents['uncorrelatederror'], locallycorrelatederror)

    def read_location_ids(self, pathname, verify_location_uuid):
        """Read location IDs only."""

        # read into NumPy structure
        contents = self.read_contents(pathname, verify_location_uuid)

        # extract location ids only
        return contents['location']

    def read_observation_count(self, pathname):
        """Read total number of observations."""

        # Compute record size
        record_size_bytes = numpy.dtype(ObservationRawBinaryWriter.build_dtype(self.basecategory, self.num_local_correlation_ranges)).itemsize

        # Open file and read past header
        file_handle = open(pathname, 'rb')
        self.read_header_from_stream(file_handle)

        # Compute count
        count = RawBinaryReader.count_remaining(file_handle, record_size_bytes)

        # Close file
        file_handle.close()
        del file_handle

        # Return result
        return count

class ObservationRawBinaryWriter(object):
    """Write portions of observation objects to binary."""
    
    def __init__(self):
        """Initialise."""
        self.formatid = FORMATID_OBSERVATION
        self.basecategory = CATEGORY_OBSERVATION

    @staticmethod
    def locally_correlated_uncertainty_field_name(index):
         return 'cor{0}'.format(index)

    @staticmethod
    def build_dtype(basecategory, num_local_correlation_ranges):
        """Build the dtype assuming a given number of locally correlated uncertainty components."""

        dtype = [ ]
        dtype.extend( DTYPE[ basecategory ] )
        for index in range(num_local_correlation_ranges):
            fieldname = ObservationRawBinaryWriter.locally_correlated_uncertainty_field_name(index)
            dtype.append( (fieldname, numpy.float64) )
        return dtype

    def write_day(self, pathname, location_uuid, obs, daynumber):
        """Write the specified object (one day of data)."""

        #  create directory if we don't already have it
        ensuredirectory(pathname)

        # open it
        file_handle = open(pathname, 'wb')
        
        # header
        file_handle.write(self.formatid)
        file_handle.write(CODEID)
        file_handle.write(location_uuid.bytes)

        # data
        valid = numpy.nonzero(numpy.logical_and(obs.mask == False, obs.time.astype(numpy.int32) == daynumber))[0]
        obscount = valid.shape[0]
        corcount = len(obs.locallycorrelatederror)

        # allocate output buffer for writing
        result = numpy.empty(obscount, dtype=ObservationRawBinaryWriter.build_dtype(self.basecategory, corcount))

        # populate output structure from observation info
        self.populate_output_structure(result, obs, valid, corcount)

        # store
        result.tofile(file_handle)

        # done
        file_handle.close()
        del file_handle

    def populate_output_structure(self, result, obs, valid, corcount):
        """Convert observations instance to NumPy structure ready for output."""

        result['location'] = obs.location[valid]
        result['measurement'] = obs.measurement[valid]
        result['uncorrelatederror'] = obs.uncorrelatederror[valid]
        for index in range(corcount):
            result[ObservationRawBinaryWriter.locally_correlated_uncertainty_field_name(index)] = obs.locallycorrelatederror[index][valid]
        return result

class LocationLookupWithID(object):
    """Location lookup table together with unique identifier."""

    def __init__(self, uuid, lookuptable):
        """Construct from given info."""

        self.uuid = uuid
        self.lookuptable = lookuptable
      
class LocationLookupRawBinaryReader(RawBinaryLocationCategoryReader):
    """Read observation locations."""

    def __init__(self):
        """Initialise."""

        super(LocationLookupRawBinaryReader, self).__init__(FORMATID_LOCATION_LOOKUP, CATEGORY_LOCATION_LOOKUP)

    
    def read(self, pathname):
        """Read locations from file handle."""

        # open it
        file_handle = open(pathname, 'rb')

        # read location uuid
        location_uuid = self.read_header_from_stream(file_handle)

        # read data
        data = numpy.fromfile(file_handle, dtype=DTYPE[ CATEGORY_LOCATION_LOOKUP ])
        lookuptable = numpy.vstack([ data['latitude'], data['longitude'] ])

        # done
        file_handle.close()
        del file_handle

        return LocationLookupWithID(location_uuid, lookuptable)

class LocationLookupRawBinaryWriter(object):
    """Write observation locations."""

    def __init__(self):
        """Initialise."""
        pass

    def write(self, pathname, data):
        """Write locations to file handle."""

        ensuredirectory(pathname)

        # open it
        file_handle = open(pathname, 'wb')

        # header including lookup table UUID
        file_handle.write(FORMATID_LOCATION_LOOKUP)
        file_handle.write(CODEID)
        file_handle.write(data.uuid.bytes)

        # lookup table
        result = numpy.empty(data.lookuptable.shape[1], dtype=DTYPE[ CATEGORY_LOCATION_LOOKUP ])
        result['latitude'] = data.lookuptable[0,:]
        result['longitude'] = data.lookuptable[1,:]
        result.tofile(file_handle)

        # done
        file_handle.close()
        del file_handle        

class LocalCorrelationRangeRawBinaryReader(RawBinaryReader):
    """Read observation locations."""

    def __init__(self):
        """Initialise."""

        super(LocalCorrelationRangeRawBinaryReader, self).__init__(FORMATID_LOCAL_CORRELATION_RANGE)

    def read(self, pathname):
        """Read locations from file handle."""

        # open it
        file_handle = open(pathname, 'rb')

        # read and check format ID
        self.read_header_from_stream(file_handle)

        # read data
        result = numpy.fromfile(file_handle, dtype=DTYPE[ CATEGORY_LOCAL_CORRELATION_RANGE ])

        # done with file
        file_handle.close()
        del file_handle
        
        return result['range']

class LocalCorrelationRangeRawBinaryWriter(object):
    """Write observation locations."""

    def __init__(self):
        """Initialise."""
        pass

    def write(self, pathname, ranges):
        """Write locations to file handle."""
        
        ensuredirectory(pathname)
        
        file_handle = open(pathname, 'wb')
        file_handle.write(FORMATID_LOCAL_CORRELATION_RANGE)
        file_handle.write(CODEID)
        result = numpy.empty(len(ranges), dtype=DTYPE[ CATEGORY_LOCAL_CORRELATION_RANGE ])
        result['range'] = ranges
        result.tofile(file_handle)
        file_handle.close()

class ObservationSourceBinaryFilenameGenerator(object):
    """Generate filenames suitable for storage of observations."""

    SUBDIRFORMAT = os.path.join('{sourcename}', '%Y')
    DATEFORMAT = '%Y%m%d'
    OBSFORMAT = '{sourcename}_{observable}_{datestring}.bin'
    CORFORMAT = '{sourcename}_{observable}_localcorrelationranges.bin'
    MOVFORMAT = '{sourcename}_{observable}_mobilelocations_{datestring}.bin'

    def __init__(self, sourcename, basepath):
        """Generate based on name of source and basepath."""

        self.sourcename = sourcename
        self.basepath = basepath

    def filename_patterns_observations(self, observable):
        """Like generate_filename but makes a generic pattern suitable for EUMOPPS."""

        subdir = ObservationSourceBinaryFilenameGenerator.SUBDIRFORMAT.format(sourcename=self.sourcename)

        name = ObservationSourceBinaryFilenameGenerator.OBSFORMAT.format(
            sourcename=self.sourcename,
            observable=observable,
            datestring=ObservationSourceBinaryFilenameGenerator.DATEFORMAT)

        return [ os.path.join(self.basepath, subdir), name ]

    def filename_local_correlation_ranges(self, observable):
        """Generate filename for correlation ranges for given observable."""

        name = ObservationSourceBinaryFilenameGenerator.CORFORMAT.format(
            sourcename=self.sourcename,
            observable=observable)
        return os.path.join(self.basepath, name)

    def filename_patterns_mobile_locations(self, observable):
        """Like generate_filename but makes a generic pattern suitable for EUMOPPS."""

        subdir = ObservationSourceBinaryFilenameGenerator.SUBDIRFORMAT.format(sourcename=self.sourcename)

        name = ObservationSourceBinaryFilenameGenerator.MOVFORMAT.format(
            sourcename=self.sourcename,
            observable=observable,
            datestring=ObservationSourceBinaryFilenameGenerator.DATEFORMAT)

        return [ os.path.join(self.basepath, subdir), name ]

    def filename_observations(self, observable, daynumber):
        """Generate individual filename for data on specified day for the specified observable."""

        patterns = self.filename_patterns_observations(observable)
        return os.path.join(*[ datetime_numeric.build_from_pattern(pattern, epoch_plus_days(daynumber)) for pattern in patterns ])

    def filename_mobile_locations(self, observable, daynumber):
        """Generate individual filename for locations."""

        patterns = self.filename_patterns_mobile_locations(observable)
        return os.path.join(*[ datetime_numeric.build_from_pattern(pattern, epoch_plus_days(daynumber)) for pattern in patterns ])


class ObservableFileSpec(object):
    """Specs needed to provide observation source interface to file."""

    def __init__(self, filename, localcorrelationranges):
        self.filename = filename
        self.localcorrelationranges = localcorrelationranges

class ObservationSourceSingleDayRawBinary(ObservationSource):
    """Observation source interface to raw binary data."""

    def __init__(self, location_lookup_with_id, filespecs, daynumber):
        """Initialise with:
           * location_lookup - known locations
           * filespecs - dictionary of ObservableFileSpec classes
           * daynumber - load this day only."""
        self.location_lookup_with_id = location_lookup_with_id
        self.filespecs = filespecs
        self.daynumber = daynumber

    def observables(self):
        """The names of variables estimated from this source."""
        return self.filespecs.keys()

    def observation_location_lookup(self):
        """NumPy array in which column number corresponds to location id and rows are latitude and longitude."""
        
        if self.location_lookup_with_id is not None:
            return self.location_lookup_with_id.lookuptable
        else:
            return None

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""
        spec = self.filespecs[observable]
        reader = ObservationRawBinaryReader(len(spec.localcorrelationranges))
        if spec.filename is not None and self.location_lookup_with_id is not None:
            return reader.read_day(spec.filename, self.location_lookup_with_id.uuid, self.daynumber)
        else:
            return ObservationsMissing(self.daynumber)
        
    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty."""
        return self.filespecs[observable].localcorrelationranges
