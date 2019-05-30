"""Unit tests for file name/path pattern parsing and searching."""

# pylint: disable=missing-docstring,invalid-name

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

import unittest
from ..namepatterns import nested_search_patterns
from ..namepatterns import NameCandidate


class TestNamePatterns(unittest.TestCase):

    def test_nested_search_patterns_MODIS(self): # noqa - allow uppercase characters in function name

        candidates = [NameCandidate('MODIS/2007/05/03/GT_MYD_2P/GT_SSD-L2-MYD11_LST_2-20070503_235500-CUOL-0.01X0.01-V1.0.nc')]
        self.assertEqual(False, candidates[0].is_matched())
        result = nested_search_patterns(
            ['MODIS/%Y/%m/%d/',
             'GT_MYD_2P/GT_SSD-L2-MYD11_LST_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V1.0.nc'],
            candidates)
        self.assertEqual(True, candidates[0].is_matched())
        self.assertEqual(1, len(result))
        self.assertEqual(
            'MODIS/2007/05/03/GT_MYD_2P/GT_SSD-L2-MYD11_LST_2-20070503_235500-CUOL-0.01X0.01-V1.0.nc', result[0].name)
        self.assertEqual({'year': '2007', 'month': '05', 'day': '03', 'hour': '23', 'minute': '55'}, result[0].fields)

        # self.assertEqual(len(result.tags), 0)
        # self.assertEqual(result.time, datetime(2007, 5, 3, 23, 55))

    def test_nested_search_patterns_ARM(self): # noqa - allow uppercase characters in function name

        candidates = [NameCandidate('ARM__/ENA_C1/enametC1.b1.20141123.000000.cdf')]
        result = nested_search_patterns(['ARM__/%A_%B/', '%a%C%B.b1.%Y%m%d.000000.cdf'], candidates)
        self.assertEqual(True, candidates[0].is_matched())
        self.assertEqual(1, len(result))
        self.assertEqual('ARM__/ENA_C1/enametC1.b1.20141123.000000.cdf', result[0].name)
        self.assertEqual({'A': 'ENA', 'B': 'C1', 'C': 'met', 'year': '2014', 'month': '11', 'day': '23'}, result[0].fields)

    def test_nested_search_patterns_AASTA(self): # noqa - allow uppercase characters in function name

        candidates = [NameCandidate('2000/02/03/20000203234800-DMI_METNO-L2P_GHRSST-STskin-GAC_polar_SST_IST-noaa14_01357_12162-v02.0-fv01.0.nc')]
        patterns = ['%Y/%m/%d/', '%Y%m%d%H%M%S-DMI_METNO-L2P_GHRSST-STskin-GAC_polar_SST_IST-%A_?????_?????-v02.0-fv01.0.nc']
        result = nested_search_patterns(patterns, candidates)
        self.assertEqual(1, len(result))
        self.assertEqual('2000/02/03/20000203234800-DMI_METNO-L2P_GHRSST-STskin-GAC_polar_SST_IST-noaa14_01357_12162-v02.0-fv01.0.nc', result[0].name)
        self.assertEqual({'A': 'noaa14', 'year': '2000', 'month': '02', 'day': '03', 'hour': '23', 'minute': '48', 'second': '00'}, result[0].fields)
