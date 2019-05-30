"""
Surface-air lake source
------------------------
Read surface-air model output for lake surface.
"""

import numpy
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from datetime import datetime
from eustace.timeutils.epoch import days_since_epoch


class ObservationSourceLakeReading(ObservationSource):
    """Provide ObservationSource interface to lake surface data loaded from text file."""

    # Observation map is not used for lakes but this variable is useful for establishing observables without
    # instantiating a class
    OBSERVATIONMAPS = { ObservationSource.TMEAN: None }
   
    # Column names
    COLUMN_DATETIME = 'DateTime'
    COLUMN_LAT = 'lat'
    COLUMN_LON = 'lon'
    COLUMN_TMEAN = 'fTa'
    
    # Date format as used in files
    DATEFORMAT = '%Y-%m-%d'

    # Offset from celsius
    TEMPERATURE_OFFSET = numpy.float32(273.15)

    @staticmethod
    def parse_date(s):
        """Helper to parse date to integer since the epoch."""
        return days_since_epoch(datetime.strptime(s, ObservationSourceLakeReading.DATEFORMAT))

    def __init__(self, filename):
        """Load and parse specified filename (tab-delimited format)."""

        super(ObservationSource, self).__init__()                

        # parse tab-delimted format
        txtdata = numpy.genfromtxt(filename, names=True, converters={ObservationSourceLakeReading.COLUMN_DATETIME: ObservationSourceLakeReading.parse_date}, dtype=numpy.float32)

        # signatures for each location so we can spot duplicates
        location_signature = \
            numpy.int64(360 * 1000000 * txtdata[ObservationSourceLakeReading.COLUMN_LAT]) + \
            numpy.int64(      1000000 * txtdata[ObservationSourceLakeReading.COLUMN_LON])

        # find duplicate signatures
        unique_signature, unique_indices, self.location_ids = numpy.unique(location_signature, return_index=True, return_inverse=True)
 
        # coordinates of unique signatures only
        self.location_lookup = numpy.vstack([ 
                txtdata[ObservationSourceLakeReading.COLUMN_LAT][unique_indices], 
                txtdata[ObservationSourceLakeReading.COLUMN_LON][unique_indices] ])

        # get time
        self.time = txtdata[ObservationSourceLakeReading.COLUMN_DATETIME]

        # compute mean in kelvin
        self.tmean = txtdata[ObservationSourceLakeReading.COLUMN_TMEAN] + ObservationSourceLakeReading.TEMPERATURE_OFFSET

        
    def observables(self):
        """Observable item is only Tmean in this case."""
        return ObservationSourceLakeReading.OBSERVATIONMAPS.keys()

    def number_of_observations(self):
        """Total number of observations."""
        return self.tmean.shape[0]

    def observation_location_lookup(self):
        """NumPy array of 3D coordinates of latitude, longitude, referencetime."""
        return self.location_lookup

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""

        if observable == ObservationSource.TMEAN:

            # Anything that's read is valid
            mask = numpy.zeros(shape=self.tmean.shape, dtype=numpy.bool)

            # Placeholder for uncorrelated error
            uncorrelatederror = numpy.empty(shape=self.tmean.shape, dtype=numpy.float32)
            uncorrelatederror.fill(5.0)

            # result
            return Observations(mask, self.time, numpy.uint64(self.location_ids), self.tmean, uncorrelatederror, [ ])

    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty (empty array in this case)."""
        
        return numpy.array([ ], numpy.float32)
