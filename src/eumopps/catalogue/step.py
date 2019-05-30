"""Classes to iterate over operations."""

from eumopps.timeutils import datetime_numeric
from eumopps.timeutils.timebase import TimeBaseDays
from eumopps.catalogue.dataset import CatalogueFileEntry
import os.path
from datetime import datetime

class Step(object):
    """Describes one step of an iterative operation."""

    def count(self):
        """Total operations."""

        raise NotImplementedError

    def timebase(self):

        raise NotImplementedError

    def time_at_index(self, operationindex):

        raise NotImplementedError

    def index_at_time(self, t):

        raise NotImplementedError

    def create_output_entry(self, patterns, operation_index):
        """Use the patterns to create and return an output entry for the specified operation."""

        raise NotImplementedError

    def is_uniquely_defined_by(self, patterns):
        """Return True if the patterns are sufficient to uniquely identify the step they correspond to, False otherwise"""

        raise NotImplementedError
        

class StepOnce(Step):

    def count(self):
        """Total operations."""

        return 1

    def index_at_time(self, t):
        """All times map to single step."""

        return 0

    def create_output_entry(self, patterns, operation_index):
        """Use the patterns to create and return an output entry for the specified operation."""

        return CatalogueFileEntry(os.path.join(*patterns))

    def is_uniquely_defined_by(self, patterns):
        """Return True if the patterns are sufficient to uniquely identify the step they correspond to.
           In our case it's always true as all patterns correspond to the single step."""

        return True

class StepDaily(Step):
    """One step per day between start and end days."""

    @staticmethod
    def parsetime(t):
        """Utility to allow convenience of strings for indicating start and end rather than datetime objects."""

        if isinstance(t, datetime):

            return t

        elif isinstance(t, str):

            return datetime_numeric.parse(t)

        else:

            raise ValueError('Input parameter to StepDaily is of unexpected type - should be datetime object or string')

    def __init__(self, start, end):
        """Initialise with text fields indicating start and end times - for convenience these may be strings, otherwise datetime objects expected."""

        self.start = StepDaily.parsetime(start)
        self.end = StepDaily.parsetime(end)

    def count(self):
        """Total operations."""

        return int(self.timebase().datetime_to_number(self.end)) + 1

    def timebase(self):

        return TimeBaseDays(self.start)

    def index_at_time(self, t):

        return int( self.timebase().datetime_to_number(t) )

    def time_at_index(self, operationindex):

        return self.timebase().number_to_datetime(operationindex)

    def filename_from_patterns(self, patterns, t):

        paths = [ datetime_numeric.build_from_pattern(pattern, t) for pattern in patterns ]
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
        return (('%Y' in concatpattern) and ('%m' in concatpattern) and ('%d' in concatpattern))
