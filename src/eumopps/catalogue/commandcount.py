"""Quickly establish number of operations records in NetCDF dataset."""

import argparse
import sys
from fileio.formatnetcdf import CatalogueReaderNetCDF


def main(argvalues=None, outputstream=sys.stdout):
    """Entry point for command-count program to
       quickly establish number of operations records in NetCDF dataset."""

    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('catalogue', help='Catalogue input file (NetCDF)')
    parser.add_argument('runmodule', help='Operation module name')
    parser.add_argument('--batchsize', type=int, default=1, help='Count batches with this batch size')
    args = parser.parse_args(argvalues)

    # run
    count = CatalogueReaderNetCDF().operationcount(args.catalogue, args.runmodule, args.batchsize)
    outputstream.write('{0}\n'.format(count))

if __name__ == '__main__':
    main()
