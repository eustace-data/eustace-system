"""
In-situ ocean source
--------------------
Observation source interface to output of in-situ ocean preprocessing.
"""

import numpy
import gzip
import shutil
import tempfile
import os
import eustaceconfig
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from datetime import datetime
from eustace.timeutils.epoch import days_since_epoch

class HadNMAT2Format(object):
    """Map from required fields to column numbers in HadNMAT2 file."""

    # Column names
    ID = 'ID'
    LAT = 'Lat'
    LON = 'Lon'
    YEAR = 'year'
    MONTH = 'month'
    DAY = 'day'
    HOUR = 'hour'
    AIRT = 'Air'  # workaround for fact that it actually says Air T

    # All fields
    FIELDS = [ ID, LAT, LON, YEAR, MONTH, DAY, HOUR, AIRT ]

    # Corresponding dtypes
    DTYPES = [ (ID, 'S8'), (LAT, 'i4'), (LON, 'i4'), (YEAR, 'i4'), (MONTH, 'i4'), (DAY, 'i4'), (HOUR, 'i4'), (AIRT, 'i4') ]

    # widths of fields
    FIELDWIDTH = 8

    # Temperature converstion
    TEMPERATURE_SCALE = numpy.float32(0.01)
    TEMPERATURE_OFFSET = numpy.float32(273.15)

    # Location scale factor for lat and lon
    LOCATION_SCALE = numpy.float32(0.1)

    def __init__(self, filename):
        """Construct from format file."""

        formatfile = open(filename, 'r')
        header = formatfile.readline()
        blankline = formatfile.readline()
        titles = formatfile.readline().split()
        self.usecols = [ titles.index(field) for field in HadNMAT2Format.FIELDS ]

class ObservationSourceInsituOceanHadNMAT2(ObservationSource):
    """Provide observation interface to HadNMAT2 gzipped text file."""
   
    # Observation maps not used here but useful to retrieve list of observables without instantiating class
    OBSERVATIONMAPS = { ObservationSource.TMEAN: None }
    
    def __init__(self, format, filename):
        """Load and parse specified filename according to specified format instance."""

        super(ObservationSource, self).__init__()                

        # assume text-gzip
        textgzip = gzip.open(filename, 'rb')

        # temporary text file
        textfile = tempfile.NamedTemporaryFile(prefix='eustace.preprocess.insitu_ocean.', suffix='.txt')

        # unzip
        shutil.copyfileobj(textgzip, textfile)

        # parse tab-delimted format
        # comments=None is required because some stations have the hashtag character in their name
        textfile.seek(0)
        txtdata = numpy.genfromtxt(textfile, delimiter=HadNMAT2Format.FIELDWIDTH, usecols=format.usecols, names=HadNMAT2Format.FIELDS, dtype=HadNMAT2Format.DTYPES, comments=None)

        # build coordinates
        self.coords = numpy.vstack([ 
                txtdata[HadNMAT2Format.LAT] * HadNMAT2Format.LOCATION_SCALE, 
                txtdata[HadNMAT2Format.LON] * HadNMAT2Format.LOCATION_SCALE ])

        # build date values
        year = txtdata[HadNMAT2Format.YEAR]
        month = txtdata[HadNMAT2Format.MONTH]
        day = txtdata[HadNMAT2Format.DAY]
        self.time = numpy.array([ days_since_epoch(datetime(year[index], month[index], day[index])) for index in range(txtdata.shape[0]) ], numpy.float32)

        # compute mean in kelvin
        self.tmean = (txtdata[HadNMAT2Format.AIRT].astype(numpy.float32) * HadNMAT2Format.TEMPERATURE_SCALE)  + HadNMAT2Format.TEMPERATURE_OFFSET

        
    def observables(self):
        """Observable item is only Tmean in this case."""
        return ObservationSourceInsituOceanHadNMAT2.OBSERVATIONMAPS.keys()

    def number_of_observations(self):
        """Total number of observations."""
        return self.tmean.shape[0]

    def observation_location_lookup(self):
        """NumPy array of 3D coordinates of latitude, longitude, referencetime."""
        return self.coords

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""

        if observable == ObservationSource.TMEAN:

            # Anything that's read is valid
            mask = numpy.zeros(shape=self.tmean.shape, dtype=numpy.bool)

            # Placeholder for uncorrelated error
            uncorrelatederror = numpy.empty(shape=self.tmean.shape, dtype=numpy.float32)
            uncorrelatederror.fill(1.1)

            # result
            return Observations(mask, self.time, numpy.array(range(self.tmean.shape[0]), numpy.uint64), self.tmean, uncorrelatederror, [ ])
        
    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty (empty array in this case)."""
        
        return numpy.array([ ], numpy.float32)
        

class ObservationSourceInsituOceanHadNMAT2Default(ObservationSourceInsituOceanHadNMAT2):
    """Instantiate HadNMAT2 readers using a default format file."""

    FORMAT = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/HadNMAT2/Docs/format.txt')

    def __init__(self, filename):

        # load format
        format = HadNMAT2Format(ObservationSourceInsituOceanHadNMAT2Default.FORMAT)

        # load file
        super(ObservationSourceInsituOceanHadNMAT2Default, self).__init__(format, filename)

