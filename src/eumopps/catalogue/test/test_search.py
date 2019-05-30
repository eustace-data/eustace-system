""" Unit tests for searching for files."""

# pylint: disable=missing-docstring,invalid-name

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

import unittest
import os
from ..search import find_all_files


class TestSearch(unittest.TestCase):

    def test_find_all_files(self):
        codepath = os.path.dirname(os.path.abspath(__file__))
        basepath = os.path.join(codepath, 'data/search')
        result = find_all_files(basepath)
        filtered_result = [
            filename for filename in result if '.svn' not in filename]
        self.assertEqual(3, len(filtered_result))
        self.assertTrue('a.somefile' in result)
        self.assertTrue('directory/b.somefile' in result)
        self.assertTrue('directory/c.somefile' in result)
