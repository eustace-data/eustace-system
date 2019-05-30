"""Test for operation parameter classes that do final resolution to single operations."""

import unittest
from tempfile import NamedTemporaryFile
from datetime import datetime
from uuid import UUID
from eumopps.catalogue.fileio.formatnetcdf import CatalogueWriterNetCDF
from eumopps.catalogue.catalogue import Catalogue
from eumopps.catalogue.storage import DataStorageFiles
from eumopps.catalogue.dataset import CatalogueDataSet, CatalogueDataSubset, CatalogueFileEntry
from eumopps.catalogue.operation import Operation, RunModule
from eumopps.catalogue.operationparameters import OperationFileReference
from eumopps.catalogue.operationparameters import OperationFileListReference
from eumopps.catalogue.operationparameters import OperationCatalogueID
from eumopps.netcdfobjects import netcdfobjects

class TestOperationFileReference(unittest.TestCase):
    
    def test_resolve_single_operation(self):
        
        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
                    name='MyExample',
                    path='/some/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['allthesame'])),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'examplefile_%Y%m%d_%H.bin']),
                            matches=[
                                CatalogueFileEntry(name='2017/examplefile_20171110_12.bin', time=datetime(2017, 10, 13, 12)),
                                CatalogueFileEntry(name='2017/examplefile_20171114_12.bin', time=datetime(2017, 11, 14, 12)),
                                CatalogueFileEntry(name='2017/examplefile_20171115_06.bin', time=datetime(2017, 11, 15, 6)),
                                CatalogueFileEntry(name='2017/examplefile_20171117_12.bin', time=datetime(2017, 11, 17, 12)) ]) ]) ])

        # Store it to temp file
        tempfile = NamedTemporaryFile(prefix='eumopps.catalogue.test.test_operationparameters.TestOperationFileReference.', suffix='.nc')
        CatalogueWriterNetCDF().save(tempfile.name, catalogue)

        # Test parameters (as if loaded from file)
        operation_parameters=[
            [ ],
            [ 0, 1, 1 ],
            [ 0, 1, 2 ],
            [ ],
            [ 0, 1, 3 ],
            [ ] ]

        # Check resolution
        self.assertIsNone(OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 0))
        self.assertEqual('/some/path/2017/examplefile_20171114_12.bin', OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 1))
        self.assertEqual('/some/path/2017/examplefile_20171115_06.bin', OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 2))
        self.assertIsNone(OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 3))
        self.assertEqual('/some/path/2017/examplefile_20171117_12.bin', OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 4))
        self.assertIsNone(OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 5))

    def test_save_load(self):

        operation_parameters=[
            [ ],
            [ 0, 1, 1 ],
            [ 0, 1, 2 ],
            [ ],
            [ 0, 1, 3 ],
            [ ] ]

        tempfile = NamedTemporaryFile(prefix='eumopps.catalogue.test.test_operationparameters.TestOperationFileReference.', suffix='.nc')

        ncfile = netcdfobjects.nc_open_for_save(tempfile.name)
        netcdfobjects.save_object(ncfile, key='testsave', value=OperationFileReference(operation_parameters))
        ncfile.close()       

        ncfile = netcdfobjects.nc_open_for_load(tempfile.name)
        result = netcdfobjects.load_object(ncfile.groups['testsave'])
        self.assertEqual(operation_parameters, result.operation_parameters)

class TestOperationFileListReference(unittest.TestCase):
    
    def test_resolve_single_operation(self):
        
        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
                    name='MyExample',
                    path='/some/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['allthesame'])),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'examplefile_%Y%m%d_%H.bin']),
                            matches=[
                                CatalogueFileEntry(name='2017/examplefile_20171110_12.bin', time=datetime(2017, 10, 13, 12)),
                                CatalogueFileEntry(name='2017/examplefile_20171114_12.bin', time=datetime(2017, 11, 14, 12)),
                                CatalogueFileEntry(name='2017/examplefile_20171115_06.bin', time=datetime(2017, 11, 15, 6)),
                                CatalogueFileEntry(name='2017/examplefile_20171117_12.bin', time=datetime(2017, 11, 17, 12)) ]) ]) ])

        # Store it to temp file
        tempfile = NamedTemporaryFile(prefix='eumopps.catalogue.test.test_operationparameters.TestOperationFileReference.', suffix='.nc')
        CatalogueWriterNetCDF().save(tempfile.name, catalogue)

        # Test parameters (as if loaded from file)
        operation_parameters=[
            [ ],
            [ 0, 1, 1, 0, 1, 2 ],
            [ ],
            [ 0, 1, 3, 0, 1, 2, 0, 1, 1 ],
            [ ] ]

        # Check resolution
        self.assertIsNone(OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 0))
        self.assertEqual([
                '/some/path/2017/examplefile_20171114_12.bin',
                '/some/path/2017/examplefile_20171115_06.bin' ],
                OperationFileListReference(operation_parameters).resolve_single_operation(tempfile.name, 1))
        self.assertIsNone(OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 2))
        self.assertEqual([
                '/some/path/2017/examplefile_20171117_12.bin',
                '/some/path/2017/examplefile_20171115_06.bin',
                '/some/path/2017/examplefile_20171114_12.bin' ],
                OperationFileListReference(operation_parameters).resolve_single_operation(tempfile.name, 3))
        self.assertIsNone(OperationFileReference(operation_parameters).resolve_single_operation(tempfile.name, 4))

class TestOperationCatalogueID(unittest.TestCase):

    def test_resolve_single_operation(self):

        # Test catalogue with known ID
        catalogue = Catalogue(identifier=str(UUID('{12345678-1234-5678-1234-567812345678}')))
        tempfile = NamedTemporaryFile(prefix='eumopps.catalogue.test.test_operationparameters.TestOperationCatalogueID.', suffix='.nc')
        CatalogueWriterNetCDF().save(tempfile.name, catalogue)       

        # Check the resolved ID matches (and operation index should be ignored here)
        result = OperationCatalogueID().resolve_single_operation(tempfile.name, 527)
        self.assertEqual('12345678-1234-5678-1234-567812345678', result)