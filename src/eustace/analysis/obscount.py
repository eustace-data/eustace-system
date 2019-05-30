"""Utility to count observations in raw binary files."""

from datetime import datetime
import argparse
import os
import json
from eumopps.timeutils import datetime_numeric
from eustace.timeutils.epoch import epoch_plus_days
from eustace.timeutils.epoch import days_since_epoch
from fileio.observationsource_rawbinary import ObservationRawBinaryReader
from fileio.observationsource_rawbinary import LocalCorrelationRangeRawBinaryReader
from fileio.observationsource_rawbinary import ObservationSourceBinaryFilenameGenerator

class ObsSummaryStatistics(object):

    def __init__(self, counts):
        """Generate summary statistics."""

        self.number_of_days = len(counts)
        self.total_observations = sum(x for x in counts if x)
        self.missing_days = sum(x is None for x in counts)
        days_not_missing = self.number_of_days - self.missing_days
        self.mean_observations_per_day = (float(self.total_observations) / days_not_missing) if days_not_missing else None

class ObsCounter(object):

    @staticmethod
    def todatetime(s):
        """If s is a string in format YYYYmmdd convert it to datetime object, otherwise return unchanged."""

        if isinstance(s, str):
            return datetime.strptime(s, ObservationSourceBinaryFilenameGenerator.DATEFORMAT)
        else:
            return s

    def __init__(self, path, source, observable, startdate, enddate):
        """Build counter instance."""

        # parse strings to datetime objects
        startdate = ObsCounter.todatetime(startdate)
        enddate = ObsCounter.todatetime(enddate)

        # set member variables
        self.path = path
        self.source = source
        self.observable = observable

        # express as day numbers
        self.startday = int( days_since_epoch(startdate) )
        self.endday = int ( days_since_epoch(enddate) )


    def count_single_day(self, num_local_correlation_ranges, daynumber):
        """Return observations on specified day or None if file is missing. Exceptions thrown if problems reading."""

        # Get pathname string
        pathname = ObservationSourceBinaryFilenameGenerator(self.source, self.path).filename_observations(self.observable, daynumber)

        # Check exists
        if not os.access(pathname, os.R_OK):
            return None

        # Count it
        return ObservationRawBinaryReader(num_local_correlation_ranges).read_observation_count(pathname)

    def read_num_local_correlation_ranges(self):
        """Read file containing number of correlation ranges."""
        
        # Get pathname string
        pathname = ObservationSourceBinaryFilenameGenerator(self.source, self.path).filename_local_correlation_ranges(self.observable)
        
        # Read it and return count
        return LocalCorrelationRangeRawBinaryReader().read(pathname).shape[0]

    def count_all_days(self):
        """Do count for all days"""

        # Need to know correlation count
        num_local_correlation_ranges = self.read_num_local_correlation_ranges()

        # Read
        return [ self.count_single_day(num_local_correlation_ranges, daynumber) for daynumber in range(self.startday, self.endday+1) ]

    def summary_statistics(self):
        """Count and get some summary statistics."""

        # Make the counts
        counts = self.count_all_days()

        # Generate stats
        return ObsSummaryStatistics(counts)
        
def main():

    # Parse arguments
    parser = argparse.ArgumentParser('obscount')
    parser.add_argument('--path', required=True, help='basepath')
    parser.add_argument('--source', required=True, help='name of source')
    parser.add_argument('--observable', required=True, help='name of observable')
    parser.add_argument('--startdate', required=True, help='start date')
    parser.add_argument('--enddate', required=True, help='end date')
    args = parser.parse_args()

    # Do the work
    print json.dumps( ObsCounter(**args.__dict__).summary_statistics().__dict__ )

if __name__ == '__main__':
    main()
