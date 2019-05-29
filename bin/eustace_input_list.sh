#!/usr/bin/env python2.7
"""Build catalogue of all data sets."""

__version__ = "$Revision: 293 $"
__author__ = "Joel R. Mitchelson"

from eustace.catalogue import generate
from eustace.catalogue.inputdescriptors import default_descriptors_filename
import sys
import argparse
import os

# parse command line
parser = argparse.ArgumentParser(os.path.basename(__file__))
parser.add_argument('-dataset', help='restrict to one data set', default=None)
parser.add_argument('-descriptors', help='descriptors file', default=default_descriptors_filename())
args = parser.parse_args()

# write to output
generate.generate_file_entries(
  outputstream=sys.stdout,
  inputstream=open(args.descriptors, 'r'),
  dataset_filter=args.dataset)

# terminate with newline
sys.stdout.write('\n')

