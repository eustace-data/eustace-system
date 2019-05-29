#!/usr/bin/python2.7
"""Parse catalogue of data sets and show statistics."""

__version__ = "$Revision: 136 $"
__author__ = "Joel R. Mitchelson"

from eustace.catalogue import generate
import sys
import argparse
import os

# parse command line
parser = argparse.ArgumentParser(os.path.basename(__file__))
parser.add_argument('filename', help='name of json file')
args = parser.parse_args()

# process
generate.generate_checksum(sys.stdout, open(args.filename, 'r'))
