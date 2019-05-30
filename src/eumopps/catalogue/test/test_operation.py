# Test operation methods

import unittest
from tempfile import NamedTemporaryFile
from eumopps.catalogue.placeholder import InputFileList
from eumopps.catalogue.placeholder import InputFile
from eumopps.catalogue.placeholder import OutputFile
from eumopps.catalogue.operation import Operation
from eumopps.catalogue.operation import RunModule
from eumopps.catalogue.operationparameters import OperationFileReference
from eumopps.catalogue.operationparameters import OperationFileListReference
from eumopps.catalogue.operationparameters import OperationCatalogueID
from eumopps.catalogue.step import StepDaily
from eumopps.catalogue.step import StepOnce
from eumopps.catalogue.catalogue import Catalogue
from eumopps.catalogue.dataset import CatalogueDataSet, CatalogueDataSubset, CatalogueFileEntry
from eumopps.catalogue.storage import DataStorageFiles
from eumopps.catalogue.fileio.formatnetcdf import CatalogueWriterNetCDF
from datetime import datetime
import os.path
from eumopps.version.svn import set_allow_unversioned_code

class OperationExample(RunModule):

    def __init__(self):
        pass

class TestOperation(unittest.TestCase):

    def setUp(self):

        # Catalogue tests check toolchain information which results in an error
        # when SVN repository has modifications.
        # This is a problem during system development because it makes it impossible
        # to see that all tests pass before committing to repository.
        # Workaround is to allow unversioned code when checking toolchain during testing.
        set_allow_unversioned_code(True)

    def test_init(self):

        a = Operation(runmodule=OperationExample(), step=StepDaily(start='20171113000000', end='20171118000000'))
        self.assertIsInstance(a.runmodule, OperationExample)
        self.assertIsInstance(a.step, StepDaily)


    def test_resolve_operation_references_single_filelist(self):
        """Example resolving references for one input data set, one filelist parameter and one output."""
        
        # Make a catalogue with one non-empty data set (as in input test)
        catalogue = Catalogue([ 

            CatalogueDataSet(
                    name='MyExample',
                    path='/some/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['allthesame'])),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'examplefile_%Y%m%d_%H.bin']),
                            matches=[
                                CatalogueFileEntry(name='2017/examplefile_20171113_12.bin', time=datetime(2017, 11, 13, 12)),
                                CatalogueFileEntry(name='2017/examplefile_20171114_12.bin', time=datetime(2017, 11, 14, 12)),
                                CatalogueFileEntry(name='2017/examplefile_20171114_18.bin', time=datetime(2017, 11, 14, 18)),
                                CatalogueFileEntry(name='2017/examplefile_20171115_12.bin', time=datetime(2017, 11, 15, 12)),
                                CatalogueFileEntry(name='2017/examplefile_20171117_12.bin', time=datetime(2017, 11, 17, 12)) ]) ]),

            ])

        
        # New datasets to be build by the operation
        newdatasets = [

            CatalogueDataSet(

                    name='SomeOutput',
                    path='/new/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['notused'])),
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['alsonotused'])),
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['%Y', 'my_output_%Y%m%d.bin']))])
            ]
    
        # Stepper across 6 days
        step = StepDaily(start='20171113000000', end='20171118000000')

        # A class with nested member structures
        op = OperationExample()
        op.information = { 
            'apples': 3, 
            'something_nested': { 
                'not_a_bird': OutputFile('SomeOutput', 2) 
                }, 
            'my_inputs': [ 
                InputFileList('MyExample', 1, 'skip') 
                ] 
            }

        # Build operation run object
        oprun = Operation(runmodule=op, step=step, newdatasets=newdatasets)

        # Resolve the example class
        result = oprun.resolve_operation_references(catalogue)

        # The static info should be unchanged
        self.assertEqual(3, result.information['apples'])
        
        # Get input and output refs
        self.assertEqual(1, len(result.information['my_inputs']))
        inputref = result.information['my_inputs'][0]
        self.assertIsInstance(inputref, OperationFileListReference)
        self.assertEqual(4, len(inputref.operation_parameters))
        self.assertEqual([ 3, 6, 3, 3 ], [ len(operation_refs) for operation_refs in inputref.operation_parameters ])
        self.assertEqual([ 0, 1, 0 ], inputref.operation_parameters[0])
        self.assertEqual([ 0, 1, 1, 0, 1, 2 ], inputref.operation_parameters[1])
        self.assertEqual([ 0, 1, 3 ], inputref.operation_parameters[2])
        self.assertEqual([ 0, 1, 4 ], inputref.operation_parameters[3])

    def test_resolve_operation_references_list_of_files(self):
        """Example resolving references for one input data set, a list of individual files, and one output."""
        
        # Make a catalogue with two data sets and presence of data like:
        # Date      :  12  13  14  15  16  17  18
        # ExampleOne:   -   +   +   +   -   +   -
        # ExampleTwo:   +   -   +   +   -   +   -
        catalogue = Catalogue([ 

            CatalogueDataSet(
                    name='ExampleOne',
                    path='/some/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['allthesame'])),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'exampleone_%Y%m%d.bin']),
                            matches=[
                                CatalogueFileEntry(name='2017/exampleone_20171113.bin', time=datetime(2017, 11, 13)),
                                CatalogueFileEntry(name='2017/exampleone_20171114.bin', time=datetime(2017, 11, 14)),
                                CatalogueFileEntry(name='2017/exampleone_20171115.bin', time=datetime(2017, 11, 15)),
                                CatalogueFileEntry(name='2017/exampleone_20171117.bin', time=datetime(2017, 11, 17)) ]) ]),
            CatalogueDataSet(
                    name='ExampleTwo',
                    path='/another/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['allthesame'])),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'exampletwo_%Y%m%d.bin']),
                            matches=[
                                CatalogueFileEntry(name='2017/exampletwo_20171112.bin', time=datetime(2017, 11, 12)),
                                CatalogueFileEntry(name='2017/exampletwo_20171114.bin', time=datetime(2017, 11, 14)),
                                CatalogueFileEntry(name='2017/exampletwo_20171115.bin', time=datetime(2017, 11, 15)),
                                CatalogueFileEntry(name='2017/exampletwo_20171117.bin', time=datetime(2017, 11, 17)) ]) ])
            ])

        
        # New datasets to be build by the operation
        newdatasets = [

            CatalogueDataSet(

                    name='SomeOutput',
                    path='/new/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['%Y', 'my_output_%Y%m%d.bin']))])
            ]
    
        # Stepper across 7 days
        step = StepDaily(start='20171112000000', end='20171118000000')

        # A class with nested member structures
        op = OperationExample()
        op.information = { 
            'apples': 3,  
            'the_result': OutputFile('SomeOutput'),
            'my_inputs': [ 
                InputFile('ExampleOne', 1, 'allowed'),
                InputFile('ExampleTwo', 1, 'allowed')
                ]
            }

        # Build operation run object
        oprun = Operation(runmodule=op, step=step, newdatasets=newdatasets)

        # Resolve the example class
        result = oprun.resolve_operation_references(catalogue)

        # The static info should be unchanged
        self.assertEqual(3, result.information['apples'])
        
        # Get input refs
        self.assertEqual(2, len(result.information['my_inputs']))
        refone = result.information['my_inputs'][0]
        reftwo = result.information['my_inputs'][1]

        # Check types and sizes
        self.assertIsInstance(refone, OperationFileReference)
        self.assertIsInstance(reftwo, OperationFileReference)
        self.assertEqual(7, len(refone.operation_parameters))
        self.assertEqual(7, len(reftwo.operation_parameters))

        ## Pattern of first data
        #self.assertEqual([ 0, 3, 3, 3, 0, 3, 0 ], [ len(operation_refs) for operation_refs in refone.operation_parameters ])
        #self.assertEqual([ 0, 3, 3, 3, 0, 3, 0 ], [ len(operation_refs) for operation_refs in refone.operation_parameters ])
        #self.assertEqual([ 0, 1, 0 ], refone.operation_parameters[1])
        #self.assertEqual([ 0, 1, 1 ], refone.operation_parameters[2])
        #self.assertEqual([ 0, 1, 2 ], refone.operation_parameters[3])
        #self.assertEqual([ 0, 1, 3 ], refone.operation_parameters[5])

        ## Pattern of second data
        #self.assertEqual([ 3, 0, 3, 3, 0, 3, 0 ], [ len(operation_refs) for operation_refs in reftwo.operation_parameters ])
        
        #self.assertEqual([ 1, 1, 0 ], reftwo.operation_parameters[0])
        #self.assertEqual([ 1, 1, 1 ], reftwo.operation_parameters[2])
        #self.assertEqual([ 1, 1, 2 ], reftwo.operation_parameters[3])
        #self.assertEqual([ 1, 1, 3 ], reftwo.operation_parameters[5])
        
        # update for missing indicator (-1) when missing data is allowed
        
        # Pattern of first data
        self.assertEqual([ 3, 3, 3, 3, 3, 3, 3 ], [ len(operation_refs) for operation_refs in refone.operation_parameters ])
        
        self.assertEqual([ 0, 1, -1 ], refone.operation_parameters[0])
        self.assertEqual([ 0, 1, 0 ], refone.operation_parameters[1])
        self.assertEqual([ 0, 1, 1 ], refone.operation_parameters[2])
        self.assertEqual([ 0, 1, 2 ], refone.operation_parameters[3])
        self.assertEqual([ 0, 1, -1 ], refone.operation_parameters[4])
        self.assertEqual([ 0, 1, 3 ], refone.operation_parameters[5])
        
        ## Pattern of second data
        self.assertEqual([ 3, 3, 3, 3, 3, 3, 3 ], [ len(operation_refs) for operation_refs in reftwo.operation_parameters ])
        
        self.assertEqual([ 1, 1, 0 ], reftwo.operation_parameters[0])
        self.assertEqual([ 1, 1, -1 ], reftwo.operation_parameters[1])
        self.assertEqual([ 1, 1, 1 ], reftwo.operation_parameters[2])
        self.assertEqual([ 1, 1, 2 ], reftwo.operation_parameters[3])
        self.assertEqual([ 1, 1, -1 ], reftwo.operation_parameters[4])
        self.assertEqual([ 1, 1, 3 ], reftwo.operation_parameters[5])

    def test_resolve_catalogue_id(self):
        """Integration test for resolving catalogue ID"""
        
        # Make a catalogue with an ID
        catalogue = Catalogue(identifier='12345678-1234-5678-1234-567812345678')
        tempfile = NamedTemporaryFile(prefix='eumopps.catalogue.test.test_operation.TestOperation.', suffix='.nc')
        CatalogueWriterNetCDF().save(tempfile.name, catalogue)

        # A class with nested member structures including catalogue ID
        opmodule = OperationExample()
        opmodule.information = { 
            'apples': 3, 
            'some_details': { 
                'my_identity': OperationCatalogueID()
                }
            }

        # Operation defined on the example class
        op = Operation(runmodule=opmodule, step=StepDaily(start='20171112000000', end='20171118000000'))

        # Resolve to list of operations
        op.resolve_operation_references(catalogue)

        # The static info should be unchanged
        # And should still have operation ID at this point
        self.assertEqual(3, op.runmodule.information['apples'])    
        self.assertIsInstance(op.runmodule.information['some_details']['my_identity'], OperationCatalogueID)

        # Then resolve to individual operation
        op.resolve_single_operation(tempfile.name, 5, 'timenow')

        # Static info should still be unchanged but should have ID now
        self.assertEqual(3, op.runmodule.information['apples'])    
        self.assertEqual('12345678-1234-5678-1234-567812345678', op.runmodule.information['some_details']['my_identity'])
