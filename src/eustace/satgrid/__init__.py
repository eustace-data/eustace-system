"""
Aggregation of satellite data onto a regular latitude-longitude grid.

Can be run from a command line.  Usage is:

  ``python2.7 -m eustace.satgrid [options] datetime``

For a full list of options type ``python2.7 -m eustace.satgrid --help``

To process a whole day of data, specify datetime in the form ``YYYYmmdd``

To process only a single time on one day, specify datetime in the form ``YYYYmmddHHMM``
"""

__version__ = "$Revision: 233 $"
__author__ = "Joel R. Mitchelson"
