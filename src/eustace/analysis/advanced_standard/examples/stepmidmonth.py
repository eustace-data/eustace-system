"""Mid-month time step needed for HadCRUT4 processing."""

from datetime import datetime
from eumopps.timeutils import datetime_numeric
from eumopps.catalogue.dataset import CatalogueFileEntry
from eumopps.catalogue.step import Step
import os

class StepMidMonth(Step):
    """Time index at middle of each month, beginning with HadCRUT4 epoch."""

    def __init__(self, epoch, stepcount):
        """Initialise with text fields indicating start and end times."""

        self.epoch = epoch
        self.stepcount = stepcount

    def count(self):
        """Total operations."""

        return self.stepcount

    def time_at_index(self, operationindex):
        """Compute the mid-point of given month number since epoch and return as datetime object.
           This is needed because HadCRUT4 uses midmonth values on the time axis."""

        year = datetime_numeric.parse(self.epoch).year
        total = 12*year + operationindex
        t0 = datetime(total / 12, 1 + (total % 12), 1)
        t1 = datetime((total + 1)/12, 1 + ((total + 1) % 12), 1)
        halfmonth = (t1 - t0) / 2
        return t0 + halfmonth

    def operation_indices_for_input_entry(self, entry):
        """Returns a list of operation indices requiring the given input."""

        if entry.time:

            epochtime = datetime_numeric.parse(self.epoch)
            return [ (12*(entry.time.year - epochtime.year) + (entry.time.month - epochtime.month)) ]

        else:

            # Non-time references assumed to apply at all steps
            return range(self.count())

    def create_output_entry(self, patterns, operation_index):
        """Use the patterns to create and return an output entry for the specified operation."""

        # Compute time at this index
        t = self.time_at_index(operation_index)

        # Generate filename
        paths = [ datetime_numeric.build_from_pattern(pattern, t) for pattern in patterns ]
        name = os.path.join(*paths)

        # Build the entry
        return CatalogueFileEntry(name, t)
