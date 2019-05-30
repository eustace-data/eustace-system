"""Conversion of datetimes into decimal years"""

import numpy
from datetime import datetime

def datetime_to_decimal_year(t_datetime):
    """Convert a datetime to decimal year."""

    # The datetime object at start of year
    year_start = datetime(t_datetime.year, 1, 1)

    # The datetime object at start of next year
    year_end = datetime(t_datetime.year + 1, 1, 1)

    # Seconds since start of year
    t_seconds = numpy.float64( (t_datetime - year_start).total_seconds() )

    # Total seconds in this year
    total_seconds = numpy.float64( (year_end - year_start).total_seconds() )

    # Compute decimal
    return numpy.float64(t_datetime.year) + (t_seconds / total_seconds)