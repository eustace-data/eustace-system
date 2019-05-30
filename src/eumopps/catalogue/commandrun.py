"""Extract single command line from catalogue."""

import argparse
import importlib
import os
import numpy
import sys
from datetime import datetime
from eumopps.version.svn import set_allow_unversioned_code
from eumopps.catalogue.fileio.formatnetcdf import CatalogueReaderNetCDF
from eumopps.jsonobjects import jsonobjects

def commandrun(timestamp, cataloguepathname, operationname, index, use_batch_size=None, batch=None, justlist=False, allow_toolchain_mismatch=False):
    """Run the command."""

    # compute operation index from batch index if requested
    if use_batch_size:
        if batch is None:
            raise ValueError('eumopps.catalogue.commandrun: use_batch_size is specified but no batch number is given')
        operationindex = (batch * use_batch_size) + index
    else:
        if batch:
            raise ValueError('eumopps.catalogue.commandrun: batch number is given but use_batch_size not specified')
        operationindex = index

    # get catalogue identifier
    identifier = CatalogueReaderNetCDF().load_identifier(cataloguepathname)

    # load the operation
    # this lists all inputs - must filter down to just this command
    operation = CatalogueReaderNetCDF().load_operation(cataloguepathname, operationname)

    if operation is not None:

        # The toolchain member of operation object is the one just loaded from catalogue
        catalogue_toolchain = operation.toolchain

        # This is the runtime toolchain evaluated from latest module info
        runtime_toolchain = operation.current_runtime_toolchain()

        try:

            # resolve pathnames (and use timestamp if required)
            operation.resolve_single_operation(cataloguepathname, operationindex, timestamp)

        except IndexError:

            raise ValueError('operation index out of range')

        if justlist:

            # Just print everything, don't actually run
            print 'Timestamp: ', timestamp
            print 'Catalogue identifier: ', identifier
            print 'Catalogue toolchain: ', catalogue_toolchain
            print 'Runtime toolchain: ', runtime_toolchain
            print 'Module to run: ',
            jsonobjects.save(sys.stdout, operation.runmodule)
            print ''

        else:

            # Check toolchain
            if (runtime_toolchain != catalogue_toolchain):

                # Report the mismatch
                message = 'eumopps.catalogue.commandrun toolchain mismatch: catalogue toolchain = {0}, runtime toolchain = {1}'.format(str(catalogue_toolchain), str(runtime_toolchain))
                
                # If mismatches allowed it's just a warning, otherwise it's an error
                if allow_toolchain_mismatch:
                    print 'WARNING: ', message
                else:
                    raise ValueError(message)

            # Run it
            return operation.run()

    else:

        # If running in batch mode we allow out-of-range operations (at tail end of last batch)
        # but otherwise it's an error
        if use_batch_size:
            print 'eumopps.catalogue.commandrun operation out of range: {0}'.format(index)
        else:
            raise ValueError('eumopps.catalogue.commandrun operation index out of range')


def main(argvalues=None):
    """Entry point to run from commandline."""

    # store timestamp first
    timestamp = datetime.now().isoformat() + ' eumopps.catalogue.commandrun ' + ' '.join(sys.argv[1:])

    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('cataloguepathname', help='Catalogue input file (NetCDF)')
    parser.add_argument('operationname', help='Module to run')
    parser.add_argument('index', type=int, help='Operation index')
    parser.add_argument('--use_batch_size', type=int, default=None, help='When this is set, the index is interpreted as an offset within a batch')
    parser.add_argument('--batch', type=int, default=None, help='Batch number when use_batch_size is set')
    parser.add_argument('--justlist', '--list', action='store_true', default=False, help='List inputs and outputs only')
    parser.add_argument('--allow_unversioned_code', action='store_true', default=False, help='Disable use of code version')
    args = parser.parse_args(argvalues)

    # disable code version if requested
    set_allow_unversioned_code(args.allow_unversioned_code)

    # run it
    commandrun(timestamp=timestamp,
               allow_toolchain_mismatch=args.allow_unversioned_code,
               cataloguepathname=args.cataloguepathname,
               operationname=args.operationname,
               index=args.index,
               use_batch_size=args.use_batch_size,
               batch=args.batch,
               justlist=args.justlist)

if __name__ == '__main__':
    main()
