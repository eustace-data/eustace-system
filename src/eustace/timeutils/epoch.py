"""Define the EUSTACE epoch and compute times relative to it."""

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

from datetime import datetime
from eumopps.timeutils.timebase import TimeBaseDays
from eumopps.timeutils.timebase import TimeBaseSeconds

EPOCH = datetime(1850, 1, 1)
"""Represents time zero for all EUSTACE products."""

# Allow single-letter names here for time etc.
# pylint: disable=invalid-name

def days_since_epoch(t):
    """Compute days to datetime t from EUSTACE epoch 01/01/1850 00:00 UTC, as 32-bit floating point."""
    return TimeBaseDays(EPOCH).datetime_to_number(t)


def epoch_plus_days(d):
    """Add d days to EUSTACE epoch 01/01/1850 00:00 UTC and return as datetime object."""
    return TimeBaseDays(EPOCH).number_to_datetime(d)


def seconds_since_epoch(t):
    """Compute seconds to datetime t from EUSTACE epoch 01/01/1850 00:00 UTC, as NumPy 64-bit integer."""
    return TimeBaseSeconds(EPOCH).datetime_to_number(t)


def epoch_plus_seconds(s):
    """Add s seconds to EUSTACE epoch 01/01/1850 00:00 UTC and return as datetime object."""
    return TimeBaseSeconds(EPOCH).number_to_datetime(s)
