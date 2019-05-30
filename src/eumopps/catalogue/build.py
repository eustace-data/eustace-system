"""Generate detailed catalogue based on descriptors."""

__version__ = "$Revision: 1019 $"
__author__ = "Joel R. Mitchelson"

import argparse
import sys
import os.path

from eumopps.catalogue.fileio.formatjson import CatalogueReaderJSON
from eumopps.catalogue.fileio.formatnetcdf import CatalogueReaderNetCDF
from eumopps.catalogue.fileio.formatnetcdf import CatalogueWriterNetCDF
from eumopps import toolchain
from eumopps.version.svn import set_allow_unversioned_code

def build_from_json_file(descriptionfilename, cataloguefilename, dataset_filter=None, update=False, overwrite=False, pathdefault=None, allow_unversioned_code=False, nochecksum=False):
    """Build catalogue file from JSON description file."""

    # option to allow unversioned code for debug purposes
    if allow_unversioned_code:
        set_allow_unversioned_code(True)

    # input must exist
    if not os.path.exists(descriptionfilename):
        raise ValueError('Specified JSON file \"{0}\" does not exist'.format(descriptionfilename))

    # output should not exist unless we are updating it
    if os.path.exists(cataloguefilename) and (not update) and (not overwrite):
        raise ValueError('Catalogue \"{0}\" already exists.  To update and overwrite an existing catalogue, use the --update option.'.format(cataloguefilename))

    # read input descripton
    inputcatalogue = CatalogueReaderJSON().load(descriptionfilename)

    # build as requested
    build(inputcatalogue, cataloguefilename, update, dataset_filter, pathdefault, nochecksum)

def build(inputcatalogue, cataloguefilename, update=False, dataset_filter=None, pathdefault=None, nochecksum=False):
    """Build catalogue file from in-memory catalogue."""

    # restrict dataset if requested
    if dataset_filter:

        # Get index of data set requested
        datasetindex = inputcatalogue.datasetindex(dataset_filter)

        # Must exist
        if datasetindex is None:
            raise ValueError('requested dataset \"{0}\" not found'.format(dataset_filter))

        # Restrict to requested
        inputcatalogue.datasets = [ inputcatalogue.dataset[datasetindex] ]

    # search for data and compute checksums
    for dataset in inputcatalogue.datasets:

        # Store version info used for search
        dataset.toolchain = toolchain.versionlist()

        # Search for files
        sys.stderr.write('eustace.catalogue.build dataset: \"' + dataset.name + '\" search\n')
        dataset.search()

        # Do checksums
        if not nochecksum:
            sys.stderr.write('eustace.catalogue.build dataset: \"' + dataset.name + '\" checksum\n')
            dataset.checksum()

    # read update if requested
    if update:

        # The original catalogue
        sys.stderr.write('eustace.catalogue.build loading: \"{0}\"\n'.format(cataloguefilename))
        outputcatalogue = CatalogueReaderNetCDF().load(cataloguefilename)

        # Append new one
        outputcatalogue.datasets.extend(inputcatalogue.datasets)
        outputcatalogue.operations.extend(inputcatalogue.operations)

    else:

        outputcatalogue = inputcatalogue

    # expand operations and add output datasets
    # note use of only those operations that came from input catalogue
    for operation in inputcatalogue.operations:

        sys.stderr.write('eustace.catalogue.build operation: \"{0}\"\n'.format(operation.name if operation.name is not None else operation.runmodule))
        operation.resolve_operation_references(outputcatalogue, pathdefault)

    # assign new identity
    outputcatalogue.newidentity()

    # save
    writer = CatalogueWriterNetCDF()
    writer.save(cataloguefilename, outputcatalogue)

def main(argvalues=None):
    """Entry point to run stand-alone."""

    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset_filter', help='restrict to one data set', default=None)
    parser.add_argument('--allow_unversioned_code', required=False, default=False, action='store_true', help='allow unversioned code for debug purposes')
    overwrite_policy = parser.add_mutually_exclusive_group()
    overwrite_policy.add_argument('--update', required=False, default=False, action='store_true', help='update existing catalogue')
    overwrite_policy.add_argument('--overwrite', required=False, default=False, action='store_true', help='update existing catalogue')
    parser.add_argument('--nochecksum', required=False, default=False, action='store_true', help='avoid computation of checksum on input data sets')
    parser.add_argument('--pathdefault', required=False, help='default path for output')
    parser.add_argument('descriptionfilename', help='pathname of input (JSON text file)')
    parser.add_argument('cataloguefilename', help='pathname of output (NetCDF catalogue file)')
    args = parser.parse_args(argvalues)

    # call build method using command line args
    build_from_json_file(**vars(args))

if __name__ == '__main__':
    main()
