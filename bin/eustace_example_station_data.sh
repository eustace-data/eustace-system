#!/usr/bin/python2.7
"""Make example station data."""

import eustace.outputformats
from eustace.outputformats.examples.example_station_data import write_example_station_data
from eumopps.version.svn import get_revision_id_for_module
from eumopps.version.svn import set_allow_unversioned_code
import os
import argparse

# parse command line
parser = argparse.ArgumentParser(os.path.basename(__file__))
parser.add_argument('-i', '--institution', required=True, dest='institution')
parser.add_argument('-u', '--allow-unversioned-code', required=False, dest='allow_unversioned_code', action='store_true')
parser.add_argument('outputdirectory', help='name of output directory')
args = parser.parse_args()
print 'Output directory: ', args.outputdirectory

# option to disable versioning
set_allow_unversioned_code(args.allow_unversioned_code)

# write file
write_example_station_data(
  source=os.path.basename(__file__),
  version=get_revision_id_for_module(eustace.outputformats),
  outputdirectory=args.outputdirectory, 
  institution=args.institution)
