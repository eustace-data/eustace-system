"""Convert datetimes to/from a time base as suitable for use in catalogue outputs."""

# Allow short names dt and t for time variables
# pylint: disable=invalid-name

__version__ = "$Revision: 1349 $"
__author__ = "Joel R. Mitchelson"

from numpy import int64
from numpy import float32
from datetime import timedelta, datetime
import numpy

from dateutil.relativedelta import relativedelta
from dateutil.rrule import rrule, MONTHLY

class TimeBase(object):
    """Base class describing interface to convert datetimes to/from a time base as suitable for use in catalogue outputs."""

    def __init__(self):
        """Initialise"""
        pass

    def datatype(self):
        """NumPy data type of timebase."""
        raise NotImplementedError

    def datetime_to_number(self, d):
        """Convert datetime object to a floating point or integer number."""
        raise NotImplementedError

    def number_to_datetime(self, n):
        """Convert number to datetime object."""
        raise NotImplementedError

    def units(self):
        """Units string to describe this timebase."""
        raise NotImplementedError


class TimeBaseDays(TimeBase):
    """Time base expressing days as a floating point number from a specified epoch."""

    def __init__(self, epoch):
        """Construct with reference to specified epoch."""
        super(TimeBaseDays, self).__init__()
        self.epoch = epoch

    def datatype(self):
        """NumPy data type of timebase."""
        return float32

    def datetime_to_number(self, d):
        """Convert datetime object to a floating point or integer number."""
        td = d - self.epoch
        return float(td.days) + (float(td.seconds) / 86400.0) + (float(td.microseconds) / 1000000.0)

    def number_to_datetime(self, n):
        """Convert number to datetime object."""
        return self.epoch + timedelta(days=n)

    def units(self):
        """Units string to describe this timebase."""

        return 'days since {0:02d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d} UTC'.format(
            self.epoch.year,
            self.epoch.month,
            self.epoch.day,
            self.epoch.hour,
            self.epoch.minute,
            self.epoch.second)


class TimeBaseSeconds(TimeBase):
    """Timebase as seconds from epoch, as 64-bit integer."""

    def __init__(self, epoch):
        """Construct with reference to specified epoch."""
        super(TimeBaseSeconds, self).__init__()
        self.epoch = epoch

    def datatype(self):
        """NumPy data type of timebase."""
        return int64

    def datetime_to_number(self, d):
        """Convert datetime object to 64-bit integer."""
        td = d - self.epoch
        return (int64(td.days) * int64(86400)) + int64(td.seconds)

    def number_to_datetime(self, n):
        """Convert integer number of seconds to datetime object."""
        return self.epoch + timedelta(seconds=long(n))

    def units(self):
        """Units string to describe this timebase."""

        return 'seconds since {0:02d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d} UTC'.format(
            self.epoch.year,
            self.epoch.month,
            self.epoch.day,
            self.epoch.hour,
            self.epoch.minute,
            self.epoch.second)


class TimeBaseAnnual(TimeBase):
    """Time base expressing days as a floating point number from a specified epoch."""

    def __init__(self, epoch):
        """Construct with reference to specified epoch."""
        super(TimeBaseAnnual, self).__init__()
        
        if TimeBaseAnnual.year_fraction(epoch) == 0.0:
            self.epoch = epoch
        else:
            raise  ValueError('epoch must be specified as the start of a year')

    def datatype(self):
        """NumPy data type of timebase."""
        return float32

    @staticmethod
    def year_fraction(d):
        """Convert a datetime to decimal year."""

        # The datetime object at start of year
        year_start = datetime(d.year, 1, 1)

        # The datetime object at start of next year
        year_end = datetime(d.year + 1, 1, 1)

        # Seconds since start of year
        t_seconds = float32( (d - year_start).total_seconds() )

        # Total seconds in this year
        total_seconds = float32( (year_end - year_start).total_seconds() )
        
        return t_seconds / total_seconds

    def datetime_to_number(self, d):
        """Convert datetime object to a floating point or integer number."""
        
        # Get year fraction for d
        d_location_in_year = TimeBaseAnnual.year_fraction(d)
        
        # Compute years complete years since epoch
        years_since_epoch = d.year - self.epoch.year
        
        # Return the total difference in decimal years
        return float(years_since_epoch + d_location_in_year)

    def number_to_datetime(self, n):
        """Convert number to datetime object accurate to second resolution."""
        
        # The datetime object at start of the year requiring that the epoch is the start of a year
        year_start = datetime(self.epoch.year + int( numpy.floor(n) ), 1, 1)
        
        # Get remaining fraction of year for n
        year_fraction = n - int( numpy.floor(n) )
        
        # Convert year_fraction to seconds
        year_end = datetime(self.epoch.year + int( numpy.floor(n) ) + 1, 1, 1)        
        total_seconds_in_year = float32( (year_end - year_start).total_seconds() )
        
        elapsed_seconds_this_year = year_fraction * total_seconds_in_year
        
        return year_start + timedelta(seconds=elapsed_seconds_this_year)

    def units(self):
        """Units string to describe this timebase."""

        return 'days since {0:02d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d} UTC'.format(
            self.epoch.year,
            self.epoch.month,
            self.epoch.day,
            self.epoch.hour,
            self.epoch.minute,
            self.epoch.second)
            
            
class TimeBaseMonthly(TimeBase):
    """Time base expressing days as a floating point number from a specified epoch.
    
    """

    def __init__(self, epoch):
        """Construct with reference to specified epoch."""
        super(TimeBaseMonthly, self).__init__()
        
        if TimeBaseMonthly.month_fraction(epoch) == 0.0:
            self.epoch = epoch
        else:
            raise  ValueError('epoch must be specified as the start of a month')

    def datatype(self):
        """NumPy data type of timebase."""
        return float32

    @staticmethod
    def year_fraction(d):
        """Convert a datetime to decimal year."""

        # The datetime object at start of year
        year_start = datetime(d.year, 1, 1)

        # The datetime object at start of next year
        year_end = datetime(d.year + 1, 1, 1)

        # Seconds since start of year
        t_seconds = float32( (d - year_start).total_seconds() )

        # Total seconds in this year
        total_seconds = float32( (year_end - year_start).total_seconds() )
        
        return t_seconds / total_seconds

    @staticmethod
    def month_fraction(d):
        """Convert a datetime to decimal month fraction."""

        # The datetime object at start of year
        month_start = datetime(d.year, d.month, 1)

        # The datetime object at start of next year
        month_end = month_start+relativedelta(months=1)

        # Seconds since start of month
        t_seconds = float32( (d - month_start).total_seconds() )

        # Total seconds in this month
        total_seconds = float32( (month_end - month_start).total_seconds() )
        
        return t_seconds / total_seconds

    def datetime_to_number(self, d):
        """Convert datetime object to a floating point or integer number number of months since epoch."""

        return (d.year * 12 + d.month) - (self.epoch.year * 12 + self.epoch.month) + TimeBaseMonthly.month_fraction(d)
        

    def number_to_datetime(self, n):
        #"""Convert number of months since epoch to datetime object accurate to second resolution."""
        
        # Get datetime for the month's start and the total number of seconds in the month
        month_start = self.epoch + relativedelta(months=int(numpy.floor(n)))
        month_end   = month_start + relativedelta(months=1)
        total_seconds_in_month = float32( (month_end - month_start).total_seconds() )
        
        # Get the time elapsed since the start of the month
        elapsed_fraction_of_month = n - numpy.floor(n)
        elapsed_seconds_this_month = elapsed_fraction_of_month * total_seconds_in_month
        
        return month_start + timedelta(seconds=elapsed_seconds_this_month)

    def units(self):
        """Units string to describe this timebase."""

        return 'days since {0:02d}-{1:02d}-{2:02d} {3:02d}:{4:02d}:{5:02d} UTC'.format(
            self.epoch.year,
            self.epoch.month,
            self.epoch.day,
            self.epoch.hour,
            self.epoch.minute,
            self.epoch.second)