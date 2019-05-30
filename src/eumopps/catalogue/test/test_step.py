"""Test the step classes."""

import unittest
from datetime import datetime
from eumopps.catalogue.step import StepOnce
from eumopps.catalogue.step import StepDaily

class TestStepOnce(unittest.TestCase):

    def test_init(self):

        s = StepOnce()

    def test_count(self):

        self.assertEqual(1, StepOnce().count())

    def test_index_at_time(self):

        self.assertEqual(0, StepOnce().index_at_time(datetime(2017, 3, 22)))
        self.assertEqual(0, StepOnce().index_at_time(datetime(1789, 8, 4)))

    def create_output_entry(self):

       entry = StepOnce().create_output_entry([ 'parentfolder', 'someinfo.file'], 6789)
       self.assertIsNone(entry.time)
       self.assertEqual('parentfolder/someinfo.file', entry.name)


class TestStepDaily(unittest.TestCase):

    def test_count(self):

        s = StepDaily(datetime(1893, 4, 2), datetime(1955, 9, 12))
        self.assertEqual(22808, s.count())

    def test_index_at_time(self):

        s = StepDaily(datetime(1893, 4, 2), datetime(1955, 9, 12))
        self.assertEqual(11629, s.index_at_time(datetime(1925, 2, 3)))


    def test_create_output_entry(self):

        s = StepDaily(datetime(1893, 4, 2), datetime(1955, 9, 12))
        entry = s.create_output_entry([ 'somefolder/%Y/', 'someinfo_%Y%m%d.nonsense' ], 11629)
        self.assertEqual('somefolder/1925/someinfo_19250203.nonsense', entry.name)
        self.assertEqual(datetime(1925, 2, 3), entry.time)

    def test_count_from_strings(self):

        s = StepDaily('18930402000000', '19550912000000')
        self.assertEqual(22808, s.count())

    def test_index_at_time_from_strings(self):

        s = StepDaily('18930402000000', '19550912000000')
        self.assertEqual(11629, s.index_at_time(datetime(1925, 2, 3)))


    def test_create_output_entry_from_strings(self):

        s = StepDaily('18930402000000', '19550912000000')
        entry = s.create_output_entry([ 'somefolder/%Y/', 'someinfo_%Y%m%d.nonsense' ], 11629)
        self.assertEqual('somefolder/1925/someinfo_19250203.nonsense', entry.name)
        self.assertEqual(datetime(1925, 2, 3), entry.time)

