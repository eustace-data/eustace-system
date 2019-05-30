""" Unit tests for loading data set descriptors and searching for matching files."""

# pylint: disable=missing-docstring,invalid-name

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

import unittest
from ..namepatterns import NameCandidate
from ..storage import DataStorageFiles
from ..dataset import CatalogueDataSubset
from datetime import datetime


class TestFileBase(unittest.TestCase):

    @staticmethod
    def simulated_file_search_simple(addpath):
        filename_list = ['A19990305.1123.somefile',
                         'Junk!!!',
                         'A18501211.0643.somefile',
                         'A20030219.1856.somefile']
        return [NameCandidate(addpath + '/' + filename) for filename in filename_list]

    def test_find_files_simple(self):
        f = DataStorageFiles(['Freda/A%Y%m%d.%H%M.somefile'])
        candidates = TestFileBase.simulated_file_search_simple('Freda')
        result = CatalogueDataSubset()
        f.match_files(result, datasetpath='Freda', candidates=candidates)
        self.assertEqual(3, len(result.matches))
        self.assertEqual(0, len(result.value_errors))
        self.assertEqual(0, len(result.archive_unused))
        self.assertEqual(datetime(1999, 3, 5, 11, 23), result.matches[0].time)
        self.assertEqual(datetime(1850, 12, 11, 6, 43), result.matches[1].time)
        self.assertEqual(datetime(2003, 2, 19, 18, 56), result.matches[2].time)
        self.assertEqual([], result.matches[0].tags)
        self.assertEqual([], result.matches[1].tags)
        self.assertEqual([], result.matches[2].tags)
        self.assertEqual(True, candidates[0].is_matched())
        self.assertEqual(False, candidates[1].is_matched())
        self.assertEqual(True, candidates[2].is_matched())
        self.assertEqual(True, candidates[3].is_matched())

    @staticmethod
    def simulated_file_search_nested(addpath):
        return [
            NameCandidate(addpath + '/MID_254/midwhatisthis254.b1.20150923.170235.cdf'),
            NameCandidate(addpath + '/MID_254/midwhatisthis254.b1.20150923.170345.cdf'),
            NameCandidate(addpath + '/MID_254/midwhatisthis254.b1.20140108.085426.cdf'),
            NameCandidate(addpath + '/MID_254/bob'),
            NameCandidate(addpath + '/notinpathatall/weird'),
            NameCandidate(addpath + '/MID_254/midanotherone254.b1.20150923.171825.cdf')
        ]

    def test_find_files_nested(self):

        candidates = TestFileBase.simulated_file_search_nested('somepath')
        f = DataStorageFiles(patterns=['somepath/%A_%B/', '%a%C%B.b1.%Y%m%d.%H%M%S.cdf'])
        result = CatalogueDataSubset()
        f.match_files(result, datasetpath='somepath', candidates=candidates)

        self.assertEqual(True, candidates[0].is_matched())
        self.assertEqual(True, candidates[1].is_matched())
        self.assertEqual(True, candidates[2].is_matched())
        self.assertEqual(False, candidates[3].is_matched())
        self.assertEqual(False, candidates[4].is_matched())
        self.assertEqual(True, candidates[5].is_matched())
        self.assertEqual(4, len(result.matches))
        self.assertEqual(0, len(result.value_errors))
        self.assertEqual(0, len(result.archive_unused))
        self.assertEqual(['MID', '254', 'whatisthis'], result.matches[0].tags)
        self.assertEqual(['MID', '254', 'whatisthis'], result.matches[1].tags)
        self.assertEqual(['MID', '254', 'whatisthis'], result.matches[2].tags)
        self.assertEqual(['MID', '254', 'anotherone'], result.matches[3].tags)
        self.assertEqual(datetime(2015, 9, 23, 17, 2, 35), result.matches[0].time)
        self.assertEqual(datetime(2015, 9, 23, 17, 3, 45), result.matches[1].time)
        self.assertEqual(datetime(2014, 1, 8, 8, 54, 26), result.matches[2].time)
        self.assertEqual(datetime(2015, 9, 23, 17, 18, 25), result.matches[3].time)

        # self.assertEqual(2, len(tags_group.groups))
        # self.assertEqual(['MID', '254', 'anotherone'], tags_group.groups[0].tags)
        # self.assertEqual(1, len(tags_group.groups[0].times))
        # self.assertEqual('20150923171825', tags_group.groups[0].times[0])
        # self.assertEqual(['MID', '254', 'whatisthis'], tags_group.groups[1].tags)
        # self.assertEqual(3, len(tags_group.groups[1].times))
        # self.assertEqual('20150923170235', tags_group.groups[1].times[0])
        # self.assertEqual('20150923170345', tags_group.groups[1].times[1])
        # self.assertEqual('20140108085426', tags_group.groups[1].times[2])
