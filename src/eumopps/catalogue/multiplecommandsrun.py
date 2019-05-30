"""Runs multiple commands lines from catalogue."""

import argparse
import importlib
import os
import numpy
import sys
from commandrun import commandrun
from datetime import datetime
from eumopps.version.svn import set_allow_unversioned_code
from eumopps.catalogue.fileio.formatnetcdf import CatalogueReaderNetCDF
from eumopps.jsonobjects import jsonobjects

def multiple_commands_run(timestamp, cataloguepathname, outeroperationmethod, inneroperationname, batch_size, batch, allow_toolchain_mismatch, justlist = None):
    """Run the same type of operation over multiple days"""

    if batch_size == 0:
            raise ValueError('eumopps.catalogue.multiple_commandsrun: batch_size has to be at least 1')

    if justlist:
	for operationindex in range(batch_size):
	    daily_inputs = commandrun(timestamp, cataloguepathname, inneroperationname, operationindex, batch_size, batch, False, allow_toolchain_mismatch)

	    print 'Batch ',batch
	    print 'Index', operationindex
	    print 'Daily inputs \n', daily_inputs
    else:
	list_of_daily_inputs = []
	for operationindex in range(batch_size):
	    daily_inputs = commandrun(timestamp, cataloguepathname, inneroperationname, operationindex, batch_size, batch, False, allow_toolchain_mismatch)  
	    list_of_daily_inputs.append(daily_inputs)
		
	[ modulename, functionname ] = outeroperationmethod.rsplit('.', 1)
	run_module = importlib.import_module(modulename)
	run_function = getattr(run_module, functionname)
	run_function(list_of_daily_inputs)


def main(argvalues=None):
    """Entry point to run from commandline."""

    # store timestamp first
    timestamp = datetime.now().isoformat() + ' eumopps.catalogue.commandrun ' + ' '.join(sys.argv[1:])

    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('cataloguepathname', help='Catalogue input file (NetCDF)')
    parser.add_argument('outeroperationmethod', help='Outer module to run')
    parser.add_argument('inneroperationname', help='Inner operation to run')
    parser.add_argument('batch_size', type=int, default=None, help='Batch size')
    parser.add_argument('batch', type=int, default=None, help='Batch number')
    parser.add_argument('--justlist', '--list', action='store_true', default=False, help='List inputs and outputs only')
    parser.add_argument('--allow_unversioned_code', action='store_true', default=False, help='Disable use of code version')
    args = parser.parse_args(argvalues)

    # disable code version if requested
    set_allow_unversioned_code(args.allow_unversioned_code)

    # run it
    multiple_commands_run(timestamp=timestamp,
               allow_toolchain_mismatch=args.allow_unversioned_code,
               cataloguepathname=args.cataloguepathname,
               outeroperationmethod=args.outeroperationmethod,
	       inneroperationname=args.inneroperationname,
               batch_size=args.batch_size,
               batch=args.batch,
               justlist=args.justlist)

if __name__ == '__main__':
    main()
