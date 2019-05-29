#!/usr/bin/python2.7
"""Build catalogue of specified data set and output by tag and sample times."""

__version__ = "$Revision: 136 $"
__author__ = "Joel R. Mitchelson"

from eustace import inputsources
import sys
import argparse
import os

# parse command line
parser = argparse.ArgumentParser(os.path.basename(__file__))
parser.add_argument('dataset', help='name of dataset')
args = parser.parse_args()

# write to output
inputsources.catalogue_generate_group_by_tag(
  outputstream=sys.stdout,
  dataset_list=inputsources.load_default_descriptors(),
  dataset_filter=args.dataset)

# terminate with newline
sys.stdout.write('\n')
