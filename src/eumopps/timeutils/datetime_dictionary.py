"""Allow datetime information to be read from dictionary of strings, for easy conversion from JSON."""

#pylint: disable=missing-docstring

from datetime import datetime

YEAR = 'year'
MONTH = 'month'
DAY = 'day'
HOUR = 'hour'
MINUTE = 'minute'
SECOND = 'second'

FIELDS = [YEAR, MONTH, DAY, HOUR, MINUTE, SECOND]


def parse(dictionary):
    """Convert a dictionary of date and time components into a datetime object."""

    # cast all to integers
    datevalues = {key: int(value) for key, value in dictionary.iteritems() if key in FIELDS}

    # allow null day
    if (YEAR in datevalues) and (MONTH in datevalues) and not (DAY in datevalues):
        datevalues[DAY] = 1
    elif (YEAR in datevalues) and not (MONTH in datevalues):
        datevalues[MONTH] = 1
        datevalues[DAY] = 1
    return datetime(**datevalues) if datevalues else None
