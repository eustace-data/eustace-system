# Test placeholder methods

import unittest

from eumopps.catalogue.placeholder import InputFileList
from eumopps.catalogue.placeholder import InputFile
from eumopps.catalogue.placeholder import OutputFile
from eumopps.catalogue.placeholder import OperationException
from eumopps.catalogue.operation import Operation
from eumopps.catalogue.operationparameters import OperationParameter
from eumopps.catalogue.operationparameters import OperationFileListReference
from eumopps.catalogue.operationparameters import OperationFileReference
from eumopps.catalogue.step import StepDaily
from eumopps.catalogue.catalogue import Catalogue
from eumopps.catalogue.dataset import CatalogueDataSet, CatalogueDataSubset, CatalogueFileEntry
from eumopps.catalogue.storage import DataStorageFiles
from datetime import datetime
import os.path

class TestInputFileList(unittest.TestCase):

    def test_init(self):

        result = InputFileList('aname', 3, InputFileList.MISSING_DATA_ALLOWED)
        self.assertEqual('aname', result.datasetname)
        self.assertEqual(3, result.subsetindex)
        self.assertEqual(InputFileList.MISSING_DATA_ALLOWED, result.missing_data)

    def test_resolve_skip(self):

        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
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
                                CatalogueFileEntry(name='2017/examplefile_20171117_12.bin', time=datetime(2017, 11, 17, 12)) ]) ]) ])

        # Make an input stepper that works across 6 days (of which the catalogue has data for 4)
        step = StepDaily(start='20171113000000', end='20171118000000')
        self.assertEqual(6, step.count())

        # Build class
        example = InputFileList('MyExample', 1, InputFileList.MISSING_DATA_SKIP)
        
        # Attempt find
        request_skip = [ False, False, False, False, False, False ]
        dataref = example.operation_input_resolve(request_skip, catalogue, step)

        # Should be a parameters object (which can later be resolved to single operation)
        self.assertIsInstance(dataref, OperationParameter)
        self.assertIsInstance(dataref, OperationFileListReference)

        # Check skip requests worked ok
        self.assertEqual([ False, False, False, True, False, True ], request_skip)

        # Should have one on 13th, 15th, 17th and two on 14th
        self.assertEqual(6, len(dataref.operation_parameters))
        self.assertEqual([ 3, 6, 3, 0, 3, 0 ], [ len(operation_refs) for operation_refs in dataref.operation_parameters ])

        # Should refer to given dataset
        self.assertEqual([ 0, 1, 0 ], dataref.operation_parameters[0])
        self.assertEqual([ 0, 1, 1, 0, 1, 2 ], dataref.operation_parameters[1])
        self.assertEqual([ 0, 1, 3 ], dataref.operation_parameters[2])
        self.assertEqual([ 0, 1, 4 ], dataref.operation_parameters[4])

    def test_resolve_span_days(self):

        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
                    name='MyExample',
                    path='/some/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['allthesame'])),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'examplefile_%Y.bin']),
                            matches=[
                                CatalogueFileEntry(name='2000/examplefile_2000.bin', time=datetime(2000, 1, 1)),
                                CatalogueFileEntry(name='2001/examplefile_2001.bin', time=datetime(2001, 1, 1)),
                                CatalogueFileEntry(name='2002/examplefile_2002.bin', time=datetime(2002, 1, 1)),
                                CatalogueFileEntry(name='2003/examplefile_2003.bin', time=datetime(2003, 1, 1)),
                                 ]) ]) ])

        # Make an input stepper that works across 6 days (of which the catalogue has data for 4)
        step = StepDaily(start='20011231000000', end='20020102000000')
        self.assertEqual(3, step.count())

        # Build class
        example = InputFileList('MyExample', 1, InputFileList.MISSING_DATA_SKIP)
        
        # Attempt find
        request_skip = [ False, False, False ]
        dataref = example.operation_input_resolve(request_skip, catalogue, step)

        # Should be a parameters object (which can later be resolved to single operation)
        self.assertIsInstance(dataref, OperationParameter)
        self.assertIsInstance(dataref, OperationFileListReference)

        # Check skip requests worked ok
        self.assertEqual([ False, False, False ], request_skip)

        # Should have one on each step
        self.assertEqual(3, len(dataref.operation_parameters))
        self.assertEqual([ 3, 3, 3 ], [ len(operation_refs) for operation_refs in dataref.operation_parameters ])

        # Should refer to given dataset
        self.assertEqual([ 0, 1, 1 ], dataref.operation_parameters[0])
        self.assertEqual([ 0, 1, 2 ], dataref.operation_parameters[1])
        self.assertEqual([ 0, 1, 2 ], dataref.operation_parameters[2])

    def test_resolve_fail_missing(self):

        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
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
                                CatalogueFileEntry(name='2017/examplefile_20171117_12.bin', time=datetime(2017, 11, 17, 12)) ]) ]) ])

        # Make an input stepper that works across 6 days (of which the catalogue has data for 4)
        step = StepDaily(start='20171113000000', end='20171118000000')
        self.assertEqual(6, step.count())

        # Build class - default behaviour should be to disallow missing data
        example = InputFileList('MyExample', 1)
        
        # Attempt find - should raise exception
        with self.assertRaises(OperationException):
            request_skip = [ False, False, False, False, False, False ]
            example.operation_input_resolve(request_skip, catalogue, step)

    def test_resolve_allow_missing(self):

        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
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
                                CatalogueFileEntry(name='2017/examplefile_20171117_12.bin', time=datetime(2017, 11, 17, 12)) ]) ]) ])

        # Make an input stepper that works across 6 days (of which the catalogue has data for 4)
        step = StepDaily(start='20171113000000', end='20171118000000')
        self.assertEqual(6, step.count())

        # Build class - default behaviour should be to disallow missing data
        example = InputFileList('MyExample', 1, missing_data=InputFileList.MISSING_DATA_ALLOWED)
        
        # Attempt find
        request_skip = [ False, False, False, False, False, False ]
        dataref = example.operation_input_resolve(request_skip, catalogue, step)

        # Should be a parameters object (which can later be resolved to single operation)
        self.assertIsInstance(dataref, OperationParameter)
        self.assertIsInstance(dataref, OperationFileListReference)
        
        # Should not request any skip even though some items lack all data
        self.assertEqual([ False, False, False, False, False, False ], request_skip)

        # Should have one on 13th, 15th, 17th and two on 14th
        self.assertEqual(6, len(dataref.operation_parameters))
        #self.assertEqual([ 3, 6, 3, 0, 3, 0 ], [ len(operation_refs) for operation_refs in dataref.operation_parameters ])
        self.assertEqual([ 3, 6, 3, 3, 3, 3 ], [ len(operation_refs) for operation_refs in dataref.operation_parameters ]) # now valid operation refs (albeit using -1 to indicate a missing inputfile)
        
        # Should refer to given dataset
        self.assertEqual([ 0, 1, 0 ], dataref.operation_parameters[0])
        self.assertEqual([ 0, 1, 1, 0, 1, 2 ], dataref.operation_parameters[1])
        self.assertEqual([ 0, 1, 3 ], dataref.operation_parameters[2])
        self.assertEqual([ 0, 1, 4 ], dataref.operation_parameters[4])
        
        # We now use -1 to indicate missing inputs to prevent operations with no inputs from being excluded from the catalogue on NetCDF write
        self.assertEqual([ 0, 1, -1 ], dataref.operation_parameters[3])
        self.assertEqual([ 0, 1, -1 ], dataref.operation_parameters[5])

class TestInputFile(unittest.TestCase):

    def test_init(self):

        a = InputFile('bob', 3)
        self.assertEqual('bob', a.datasetname)
        self.assertEqual(3, a.subsetindex)
        self.assertEqual('not_allowed', a.missing_data)

        b = InputFileList('freda', 7, InputFile.MISSING_DATA_SKIP)
        self.assertEqual('freda', b.datasetname)
        self.assertEqual(7, b.subsetindex)
        self.assertEqual('skip', b.missing_data)

    def test_resolve_fail_multiple(self):
        """Check it fails if there are multiple inputs"""

        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
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
                                CatalogueFileEntry(name='2017/examplefile_20171117_12.bin', time=datetime(2017, 11, 17, 12)) ]) ]) ])

        # Make an input stepper that works across 6 days (of which the catalogue has data for 4)
        step = StepDaily(start='20171113000000', end='20171118000000')
        self.assertEqual(6, step.count())

        # Build class
        example = InputFile('MyExample', 1, InputFileList.MISSING_DATA_SKIP)
        
        # Attempt find
        request_skip = [ False, False, False, False, False, False ]
        with self.assertRaises(OperationException):
            example.operation_input_resolve(request_skip, catalogue, step)

    def test_resolve_skip(self):
        """Check it fails if there are multiple inputs"""

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

        # Make an input stepper that works across 6 days (of which the catalogue has data for 4)
        step = StepDaily(start='20171113000000', end='20171118000000')
        self.assertEqual(6, step.count())

        # Build class
        example = InputFile('MyExample', 1, InputFileList.MISSING_DATA_SKIP)

        # Attempt find
        request_skip = [ False, False, False, False, False, False ]
        dataref = example.operation_input_resolve(request_skip, catalogue, step)

        # Should be a parameters object (which can later be resolved to single operation)
        self.assertIsInstance(dataref, OperationParameter)
        self.assertIsInstance(dataref, OperationFileReference)

        # Check skip requested
        self.assertEqual([ True, False, False, True, False, True ], request_skip)

        # Should have one on 14th, 15th, 17th
        # There's also one on 10th but that isn't used
        self.assertEqual(6, len(dataref.operation_parameters))
        self.assertEqual([ 0, 3, 3, 0, 3, 0 ], [ len(operation_refs) for operation_refs in dataref.operation_parameters ])

        # Should refer to given dataset
        self.assertEqual([ 0, 1, 1 ], dataref.operation_parameters[1])
        self.assertEqual([ 0, 1, 2 ], dataref.operation_parameters[2])
        self.assertEqual([ 0, 1, 3 ], dataref.operation_parameters[4])

    def test_resolve_allow_missing(self):
        """Check it fails if there are multiple inputs"""

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

        # Make an input stepper that works across 6 days (of which the catalogue has data for 4)
        step = StepDaily(start='20171113000000', end='20171118000000')
        self.assertEqual(6, step.count())

        # Build class
        example = InputFile('MyExample', 1, InputFileList.MISSING_DATA_ALLOWED)
        
        # Attempt find
        request_skip = [ False, False, False, False, False, False ]
        dataref = example.operation_input_resolve(request_skip, catalogue, step)
        
        # Should be a parameters object (which can later be resolved to single operation)
        self.assertIsInstance(dataref, OperationParameter)
        self.assertIsInstance(dataref, OperationFileReference)
        
        # Check no skip requested
        self.assertEqual([ False, False, False, False, False, False ], request_skip)

        # Should have one on 14th, 15th, 17th
        # There's also one on 10th but that isn't used
        self.assertEqual(6, len(dataref.operation_parameters))
        #self.assertEqual([ 0, 3, 3, 0, 3, 0 ], [ len(operation_refs) for operation_refs in dataref.operation_parameters ])
        self.assertEqual([ 3, 3, 3, 3, 3, 3 ], [ len(operation_refs) for operation_refs in dataref.operation_parameters ]) # There should now always be at least 3 operation refs for MISSING_DATA_ALLOWED

        # Should refer to given dataset
        self.assertEqual([ 0, 1, 1 ], dataref.operation_parameters[1])
        self.assertEqual([ 0, 1, 2 ], dataref.operation_parameters[2])
        self.assertEqual([ 0, 1, 3 ], dataref.operation_parameters[4])

class TestOutputFile(unittest.TestCase):

    def test_init(self):

         a = OutputFile('fflp', 7)
         self.assertEqual('fflp', a.datasetname)
         self.assertEqual(7, a.subsetindex)

    def test_operation_output_resolve(self):

        # Make a catalogue with one non-empty data set (as in input test)
        # and one empty data set ready for outputs
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

            CatalogueDataSet(
                    name='SomeOutput',
                    path='/new/path',
                    subsets=[
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['notused'])),
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['alsonotused'])),
                        CatalogueDataSubset(layout=DataStorageFiles(patterns=['%Y', 'my_output_%Y%m%d.bin']))])
            ])
    
        # Stepper across 6 days
        step = StepDaily(start='20171113000000', end='20171118000000')

        # This should output to subset 2 of second dataset
        example_output = OutputFile('SomeOutput', 2)

        # Request indices 2, 3, 5
        result = example_output.operation_output_resolve(catalogue, step, [ 2, 3, 5 ])

        # Should populate the catalogue with matches
        self.assertEqual(3, len(catalogue.datasets[1].subsets[2].matches))

        # Check the details appended to catalogue
        self.assertEqual(('2017', 'my_output_20171115.bin'), os.path.split(catalogue.datasets[1].subsets[2].matches[0].name))
        self.assertEqual(('2017', 'my_output_20171116.bin'), os.path.split(catalogue.datasets[1].subsets[2].matches[1].name))
        self.assertEqual(('2017', 'my_output_20171118.bin'), os.path.split(catalogue.datasets[1].subsets[2].matches[2].name))
        self.assertEqual(datetime(2017, 11, 15), catalogue.datasets[1].subsets[2].matches[0].time)
        self.assertEqual(datetime(2017, 11, 16), catalogue.datasets[1].subsets[2].matches[1].time)
        self.assertEqual(datetime(2017, 11, 18), catalogue.datasets[1].subsets[2].matches[2].time)

        # Also check the reference list correctly describes these
        self.assertEqual(3, len(result.operation_parameters))
        self.assertEqual([ 3, 3, 3 ], [ len(operation_refs) for operation_refs in result.operation_parameters ])
        self.assertEqual([ 1, 2, 0 ], result.operation_parameters[0])
        self.assertEqual([ 1, 2, 1 ], result.operation_parameters[1])
        self.assertEqual([ 1, 2, 2 ], result.operation_parameters[2])
