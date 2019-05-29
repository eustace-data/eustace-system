#!/usr/bin/python2.7
"""Parse catalogue of data sets and show statistics."""

__version__ = "$Revision: 293 $"
__author__ = "Joel R. Mitchelson"

from eustace.catalogue import stats
import sys
import argparse
import os

# parse command line
parser = argparse.ArgumentParser(os.path.basename(__file__))
parser.add_argument('filename', help='name of json file')
args = parser.parse_args()

# process
stats.catalogue_stats(sys.stdout, open(args.filename, 'r'))
