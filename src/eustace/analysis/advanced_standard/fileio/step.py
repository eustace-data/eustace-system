"""Classes to iterate over with indices defined as days since the analysis epoch."""

from eumopps.timeutils import datetime_numeric
from eumopps.catalogue.step import StepDaily
from eumopps.timeutils.timebase import TimeBaseDays
from eumopps.catalogue.dataset import CatalogueFileEntry
from eustace.timeutils.epoch import EPOCH
import os.path
from datetime import datetime

class StepDailyEpoch(StepDaily):
    """One step per day between start and end days with indexes defined as days since the analysis reference epoch."""
    
    def __init__(self, start, end):
        """Initialise with text fields indicating start and end times - for convenience these may be strings, otherwise datetime objects expected."""
        
        super(StepDailyEpoch, self).__init__(start, end)
        
    def count(self):
        """Total operations."""
        print int(self.timebase().datetime_to_number(self.end) - self.timebase().datetime_to_number(self.start)) + 1
        return int(self.timebase().datetime_to_number(self.end) - self.timebase().datetime_to_number(self.start)) + 1

    def timebase(self):
        return TimeBaseDays(EPOCH)
        
    def create_output_entry(self, patterns, operation_index):
        """Use the patterns to create and return an output entry for the specified operation."""
        print patterns, operation_index

        # Compute time at this index
        t = self.time_at_index(operation_index)

        # Build the entry
        return CatalogueFileEntry(self.filename_from_patterns(patterns, t), t)
        
    def filename_from_patterns(self, patterns, t):
        
        paths = [ datetime_numeric.build_from_pattern(pattern, t) for pattern in patterns ]
        name = os.path.join(*paths)
        return name
        
    def index_at_time(self, t):
        return int( self.timebase().datetime_to_number(t) )

    def time_at_index(self, operationindex):
        return self.timebase().number_to_datetime(operationindex)