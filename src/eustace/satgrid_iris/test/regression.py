"""Test for compatibility with previously produced files.
   This is for a manual check - not part of test suite."""

import argparse
from netCDF4 import Dataset
import os
import numpy

# Get filenames from command line
parser = argparse.ArgumentParser('satgrid')
parser.add_argument('--showdiags', action='store_true')
parser.add_argument('filename_old')
parser.add_argument('filename_new')
args = parser.parse_args()

# List of the two filenames
filenames = [ args.filename_old, args.filename_new ]

# File sizes should only differ by a few characters
# due to different global attributes
# assert(abs(os.stat(filenames[0]).st_size - os.stat(filenames[1]).st_size) < 20000)

# Load NetCDF and get handles to variables
datasets = [ Dataset(filename, 'r') for filename in filenames ]
a = datasets[0].variables
b = datasets[1].variables

# Some diagnostics
if args.showdiags:

    mismatch_numobs = numpy.flatnonzero(a['total_number_of_observations'][:] != b['total_number_of_observations'][:])
    print mismatch_numobs
    print a['total_number_of_observations'][:].flat[mismatch_numobs]
    print b['total_number_of_observations'][:].flat[mismatch_numobs]

    mismatch_validobs = numpy.flatnonzero(a['ts_number_of_observations'][:] != b['ts_number_of_observations'][:])
    print mismatch_validobs
    print a['ts_number_of_observations'][:].flat[mismatch_validobs]
    print b['ts_number_of_observations'][:].flat[mismatch_validobs]

    mismatch_sfc = numpy.flatnonzero(a['tsmean_unc_loc_sfc'][:].mask != b['tsmean_unc_loc_sfc'][:].mask)
    print mismatch_sfc
    print a['tsmean_unc_loc_sfc'][:].mask.flat[mismatch_sfc]
    print b['tsmean_unc_loc_sfc'][:].mask.flat[mismatch_sfc]
    print a['tsmean_unc_loc_sfc'][:].data.flat[mismatch_sfc]
    print b['tsmean_unc_loc_sfc'][:].data.flat[mismatch_sfc]

    mismatch_atm = numpy.flatnonzero(a['tsmean_unc_loc_atm'][:].mask != b['tsmean_unc_loc_atm'][:].mask)
    print mismatch_atm
    print a['tsmean_unc_loc_atm'][:].mask.flat[mismatch_atm]
    print b['tsmean_unc_loc_atm'][:].mask.flat[mismatch_atm]
    print a['tsmean_unc_loc_atm'][:].data.flat[mismatch_atm]
    print b['tsmean_unc_loc_atm'][:].data.flat[mismatch_atm]

# Test values
numpy.testing.assert_equal(a['latitude'][:], b['latitude'][:])
numpy.testing.assert_equal(a['longitude'][:], b['longitude'][:])
numpy.testing.assert_equal(a['total_number_of_observations'][:], b['total_number_of_observations'][:])
numpy.testing.assert_equal(a['ts_number_of_observations'][:], b['ts_number_of_observations'][:])
numpy.testing.assert_equal(a['tsmean'][:].mask, b['tsmean'][:].mask)
numpy.testing.assert_equal(a['tsmean'][:].data, b['tsmean'][:].data)
numpy.testing.assert_equal(a['tsmax'][:].mask, b['tsmax'][:].mask)
numpy.testing.assert_equal(a['tsmax'][:].data, b['tsmax'][:].data)
numpy.testing.assert_equal(a['tsmin'][:].mask, b['tsmin'][:].mask)
numpy.testing.assert_equal(a['tsmin'][:].data, b['tsmin'][:].data)
numpy.testing.assert_equal(a['tsmean_unc_ran'][:].mask, b['tsmean_unc_ran'][:].mask)
numpy.testing.assert_equal(a['tsmean_unc_ran'][:].data, b['tsmean_unc_ran'][:].data)
numpy.testing.assert_equal(a['tsvariance'][:].mask, b['tsvariance'][:].mask)
numpy.testing.assert_equal(a['tsvariance'][:].data, b['tsvariance'][:].data)
numpy.testing.assert_equal(a['tsmean_unc_spl'][:].mask, b['tsmean_unc_spl'][:].mask)
numpy.testing.assert_equal(a['tsmean_unc_spl'][:].data, b['tsmean_unc_spl'][:].data)
numpy.testing.assert_equal(a['tsmean_unc_loc_sfc'][:].mask, b['tsmean_unc_loc_sfc'][:].mask)
numpy.testing.assert_equal(a['tsmean_unc_loc_sfc'][:].data, b['tsmean_unc_loc_sfc'][:].data)
numpy.testing.assert_equal(a['tsmean_unc_loc_atm'][:].mask, b['tsmean_unc_loc_atm'][:].mask)
numpy.testing.assert_equal(a['tsmean_unc_loc_atm'][:].data, b['tsmean_unc_loc_atm'][:].data)
numpy.testing.assert_equal(a['tsmean_unc_sys'][:].mask, b['tsmean_unc_sys'][:].mask)
numpy.testing.assert_equal(a['tsmean_unc_sys'][:].data, b['tsmean_unc_sys'][:].data)
