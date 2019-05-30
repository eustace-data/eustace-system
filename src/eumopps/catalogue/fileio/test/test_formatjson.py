"""Unit tests for file name/path pattern parsing and searching."""

# pylint: disable=missing-docstring,invalid-name

__version__ = "$Revision: 1335 $"
__author__ = "Joel R. Mitchelson"

import unittest
from ..formatjson import CatalogueReaderJSON
from ...storage import DataStorageFiles
from ...storage import DataStorageTarGZip
import StringIO

class TestCatalogueReaderJSON(unittest.TestCase):

    def test_files(self):
        inputstring = \
            '{ "python_class": "eumopps.catalogue.catalogue.Catalogue", ' \
            '"datasets": [ ' \
            '  { "python_class": "eumopps.catalogue.dataset.CatalogueDataSet", ' \
            '    "path": "/gws/nopw/j04/eustace/data/incoming/MODIS", "name": "Aqua-MODIS", ' \
            '    "subsets" : ' \
            '    [ { "python_class": "eumopps.catalogue.dataset.CatalogueDataSubset", ' \
            '        "layout": '\
            '        { ' \
            '          "python_class": "eumopps.catalogue.storage.DataStorageFiles", ' \
            '          "patterns": [ "%Y/%m/%d/", "GT_MYD_2P/GT_SSD-L2-MYD11_LST_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V1.0.nc" ] ' \
            '        }' \
            '     } ] ' \
            '  } ] ' \
            '}'

        dataset_list = CatalogueReaderJSON().loadstream(StringIO.StringIO(inputstring)).datasets
        self.assertEqual(1, len(dataset_list))
        d = dataset_list[0]
        self.assertEqual('Aqua-MODIS', d.name)
        self.assertEqual('/gws/nopw/j04/eustace/data/incoming/MODIS', d.path)
        self.assertTrue(d.subsets)
        self.assertEqual(1, len(d.subsets))
        self.assertIsInstance(d.subsets[0].layout, DataStorageFiles)
        self.assertEqual(['%Y/%m/%d/', 'GT_MYD_2P/GT_SSD-L2-MYD11_LST_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V1.0.nc'], d.subsets[0].layout.patterns)

    def test_tar_gzip(self):

        inputstring = \
            '{ "python_class": "eumopps.catalogue.catalogue.Catalogue", ' \
            '  "datasets": [ ' \
            '  { "python_class": "eumopps.catalogue.dataset.CatalogueDataSet", ' \
            '    "path": "/Some/Path/Here", "name": "Bob", ' \
            '    "subsets" : ' \
            '    [ { "python_class": "eumopps.catalogue.dataset.CatalogueDataSubset", ' \
            '        "layout": '\
            '        { ' \
            '          "python_class": "eumopps.catalogue.storage.DataStorageTarGZip", ' \
            '          "archivepathname": "SomeSubFolder/OtherFolder/%Y%m%d_myname.tar.gz", ' \
            '          "internalpathname": "a/b/c.%Y%m%d.%H%M.nc" ' \
            '        }' \
            '     } ] ' \
            '  } ] ' \
            '}'

        dataset_list = CatalogueReaderJSON().loadstream(StringIO.StringIO(inputstring)).datasets
        self.assertEqual(1, len(dataset_list))
        d = dataset_list[0]
        self.assertEqual('Bob', d.name)
        self.assertEqual('/Some/Path/Here', d.path)
        self.assertTrue(d.subsets)
        self.assertEqual(1, len(d.subsets))
        self.assertIsInstance(d.subsets[0].layout, DataStorageTarGZip)
        self.assertEqual('SomeSubFolder/OtherFolder/%Y%m%d_myname.tar.gz', d.subsets[0].layout.pathname)
        self.assertEqual('a/b/c.%Y%m%d.%H%M.nc', d.subsets[0].layout.internalpathname)
