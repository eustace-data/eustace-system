from eustace.analysis.advanced_standard.examples.placeholder import AnalysisStepIndex, AnnualBatchDays
from eustace.analysis.advanced_standard.examples.placeholder import AnnualBatchDaysOutput
from eumopps.catalogue.step import StepDaily
from eustace.timeutils.epoch import EPOCH
from datetime import timedelta

import unittest

from eumopps.catalogue.catalogue import Catalogue
from eumopps.catalogue.dataset import CatalogueDataSet, CatalogueDataSubset, CatalogueFileEntry
from eumopps.catalogue.storage import DataStorageFiles

import os.path

from eustace.analysis.advanced_standard.examples.step import StepAnnual

from datetime import datetime

class TestAnalysisStepIndex(unittest.TestCase):
    
    def test_step(self):
        """Check that AnalysisStepIndex generates operation parameter indices as days since the analysis EPOCH"""
        
        operation_arg_maker = AnalysisStepIndex()
        
        # Get the analysis reference date
        reference_date = EPOCH
        
        # Dummy variables that are passed to AnalysisStepIndex but not used
        unused_request_skip = None
        unused_catalogue = None
        
        # One day, 1000 days from EPOCH
        shift_1000_start = reference_date+timedelta(days = 1000)
        shift_1000_end = reference_date+timedelta(days = 1000)
        step_1000 = StepDaily(start=shift_1000_start, end=shift_1000_end)
        
        operation_parameters_1000 = operation_arg_maker.operation_input_resolve(unused_request_skip, unused_catalogue, step_1000)
        self.assertItemsEqual(operation_parameters_1000.operation_parameters, [1000,])

        # Three day sequence, 2000, 2001, 2002 days from EPOCH
        shift_2000_start = reference_date+timedelta(days = 2000)
        shift_2000_end = reference_date+timedelta(days = 2002)
        step_2000 = StepDaily(start=shift_2000_start, end=shift_2000_end)
        
        operation_parameters_2000 = operation_arg_maker.operation_input_resolve(unused_request_skip, unused_catalogue, step_2000)
        self.assertItemsEqual(operation_parameters_2000.operation_parameters, [2000, 2001, 2002])
        
        
class TestAnnualBatchDays(unittest.TestCase):
    
    def test_find_references(self):
        
        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet(
                    name='MyExampleDataset',
                    path='/some/path',
                    subsets=[
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['allthesame.bin']),
                             matches=[
                                CatalogueFileEntry(name='allthesame.bin', time=datetime(2000, 1, 1))],
                            ),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'examplefile_%Y.bin']),
                            matches=[
                                CatalogueFileEntry(name='2000/examplefile_2000.bin', time=datetime(2000, 1, 1)),
                                CatalogueFileEntry(name='2001/examplefile_2001.bin', time=datetime(2001, 1, 1)),
                                CatalogueFileEntry(name='2002/examplefile_2002.bin', time=datetime(2002, 1, 1)),
                                CatalogueFileEntry(name='2003/examplefile_2003.bin', time=datetime(2003, 1, 1)),
                                 ]) ]) ])
                                 
        catalogue = Catalogue([ CatalogueDataSet(
                    name='MyExampleDataset',
                    path='/some/path',
                    subsets=[
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['allthesame.bin']),
                             matches=[
                                CatalogueFileEntry(name='allthesame.bin', time=None)],
                                ),
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['%Y', 'examplefile_%Y%m%d.bin']),
                            matches=[
                                CatalogueFileEntry(name='2001/examplefile_20010101.bin', time=datetime(2001, 1, 1)),
                                CatalogueFileEntry(name='2001/examplefile_20010101.bin', time=datetime(2001, 1, 2)),
                                CatalogueFileEntry(name='2001/examplefile_20010101.bin', time=datetime(2001, 1, 3)),
                                CatalogueFileEntry(name='2001/examplefile_20010101.bin', time=datetime(2002, 1, 2)),
                                 ]) ]) ])
        
        
        step = StepAnnual(start = datetime(2001, 01, 01), end = datetime(2002, 02, 01))
        print step.count()
        inputhandler = AnnualBatchDays('MyExampleDataset', subsetindex=1, missing_data=AnnualBatchDays.MISSING_DATA_ALLOWED)
        request_skip = [False] * step.count()
        print inputhandler.find_references(request_skip, catalogue, step)
        print inputhandler
        
        
        # now need to test day number computation
        
        
class TestAnnualBatchDaysOutput(unittest.TestCase):
    
    def test_operation_output_resolve(self):
        
        # Make a catalogue with one data set and two subsets
        catalogue = Catalogue([ CatalogueDataSet( name='MyExampleDataset',
                                                  path='/some/path',
                                                  subsets=[ CatalogueDataSubset(layout=DataStorageFiles(patterns=['%Y', 'examplefile_%Y%m%d.bin']) ) ],
                                                 )
                               ])
    
        step = StepAnnual(start = datetime(2001, 12, 30), end = datetime(2002, 01, 02))
        
        output_placeholder = AnnualBatchDaysOutput('MyExampleDataset', 0)
        print output_placeholder.subsetindex
        operation_indices = [0,1]
        output_entries = output_placeholder.operation_output_resolve(catalogue, step, operation_indices)
        
        print output_entries
        print zip([(entry.name, entry.time) for entry in catalogue.datasets[0].subsets[0].matches])
        
        # Check the details appended to catalogue
        self.assertEqual(('2001', 'examplefile_20011230.bin'), os.path.split(catalogue.datasets[0].subsets[0].matches[0].name))
        self.assertEqual(('2001', 'examplefile_20011231.bin'), os.path.split(catalogue.datasets[0].subsets[0].matches[1].name))
        self.assertEqual(('2002', 'examplefile_20020101.bin'), os.path.split(catalogue.datasets[0].subsets[0].matches[2].name))
        self.assertEqual(('2002', 'examplefile_20020102.bin'), os.path.split(catalogue.datasets[0].subsets[0].matches[3].name))
        self.assertEqual(datetime(2001, 12, 30), catalogue.datasets[0].subsets[0].matches[0].time)
        self.assertEqual(datetime(2001, 12, 31), catalogue.datasets[0].subsets[0].matches[1].time)
        self.assertEqual(datetime(2002, 01, 01), catalogue.datasets[0].subsets[0].matches[2].time)
        self.assertEqual(datetime(2002, 01, 02), catalogue.datasets[0].subsets[0].matches[3].time)
    