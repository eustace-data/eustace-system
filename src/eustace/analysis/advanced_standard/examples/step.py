"""Mid-month time step needed for HadCRUT4 processing."""

from eumopps.catalogue.step import Step
from eumopps.timeutils import datetime_numeric
from eumopps.timeutils.timebase import TimeBaseAnnual, TimeBaseMonthly
from eumopps.catalogue.dataset import CatalogueFileEntry
import os.path
from datetime import datetime
import numpy

class StepAnnual(Step):
    """One step per year between start and end dates."""

    @staticmethod
    def parsetime(t):
        """Utility to allow convenience of strings for indicating start and end rather than datetime objects."""

        if isinstance(t, datetime):

            return t

        elif isinstance(t, str):

            return datetime_numeric.parse(t)

        else:

            raise ValueError('Input parameter to StepAnnual is of unexpected type - should be datetime object or string')

    def __init__(self, start, end):
        """Initialise with text fields indicating start and end times - for convenience these may be strings, otherwise datetime objects expected."""

        self.start = StepAnnual.parsetime(start)
        self.end = StepAnnual.parsetime(end)

    def count(self):
        """Total operations."""

        return int(self.timebase().datetime_to_number(self.end)) + 1

    def timebase(self):
        """Get timebase as the start of the year containing self.start"""
        
        return TimeBaseAnnual( datetime( self.start.year, 1, 1, 0, 0 , 0 ) )

    def index_at_time(self, t):

        return int( numpy.floor(self.timebase().datetime_to_number(t)) )

    def time_at_index(self, operationindex):

        return self.timebase().number_to_datetime(operationindex)

    def filename_from_patterns(self, patterns, t):

        paths = [ datetime_numeric.build_from_pattern(pattern, t) for pattern in patterns ]
        
        #paths = [ datetime_numeric.build_from_pattern(pattern, t) if (t >= self.start or t <= self.end) for pattern in patterns ] # would require a list of times provided by time_at_index
        name = os.path.join(*paths)
        return name

    def create_output_entry(self, patterns, operation_index):
        """Use the patterns to create and return an output entry for the specified operation."""

        # Compute time at this index
        t = self.time_at_index(operation_index)

        # Build the entry
        return CatalogueFileEntry(self.filename_from_patterns(patterns, t), t)
    
    def is_uniquely_defined_by(self, patterns):
        """Return True if the patterns are sufficient to uniquely identify the step they correspond to, False otherwise"""

        concatpattern = ''.join(patterns)
        return (('%Y' in concatpattern))

class StepMonthly(Step):
    """One step per day between start and end days."""

    @staticmethod
    def parsetime(t):
        """Utility to allow convenience of strings for indicating start and end rather than datetime objects."""

        if isinstance(t, datetime):

            return t

        elif isinstance(t, str):

            return datetime_numeric.parse(t)

        else:

            raise ValueError('Input parameter to StepAnnual is of unexpected type - should be datetime object or string')

    def __init__(self, start, end):
        """Initialise with text fields indicating start and end times - for convenience these may be strings, otherwise datetime objects expected."""

        self.start = StepMonthly.parsetime(start)
        self.end = StepMonthly.parsetime(end)
        
        print self.start, self.end
        print self.timebase().datetime_to_number(self.start), self.timebase().datetime_to_number(self.end)
        

    def count(self):
        """Total operations."""

        return int(self.timebase().datetime_to_number(self.end)) + 1

    def timebase(self):
        """Get timebase as the start of the month containing self.start"""
        
        return TimeBaseMonthly( datetime( self.start.year, self.start.month, 1, 0, 0 , 0 ) )

    def index_at_time(self, t):

        return int( numpy.floor(self.timebase().datetime_to_number(t)) )

    def time_at_index(self, operationindex):

        return self.timebase().number_to_datetime(operationindex)

    def filename_from_patterns(self, patterns, t):

        paths = [ datetime_numeric.build_from_pattern(pattern, t) for pattern in patterns ]
        
        #paths = [ datetime_numeric.build_from_pattern(pattern, t) if (t >= self.start or t <= self.end) for pattern in patterns ] # would require a list of times provided by time_at_index
        name = os.path.join(*paths)
        return name

    def create_output_entry(self, patterns, operation_index):
        """Use the patterns to create and return an output entry for the specified operation."""

        # Compute time at this index
        t = self.time_at_index(operation_index)

        # Build the entry
        return CatalogueFileEntry(self.filename_from_patterns(patterns, t), t)
    
    def is_uniquely_defined_by(self, patterns):
        """Return True if the patterns are sufficient to uniquely identify the step they correspond to, False otherwise"""
        print patterns
        concatpattern = ''.join(patterns)
        return (('%Y' in concatpattern) and ('%m' in concatpattern))
