"""Test SVN communications."""

# pylint: disable=missing-docstring

import unittest
from .. import svn
from .. import codepath


class TestSVN(unittest.TestCase):

    def test_repository_state_exception(self):
        try:
            raise svn.RepositoryStateException()
        except svn.RepositoryStateException as result:
            self.assertEqual(
                'Local copy of SVN has modifications or is in a mixed state.', str(result))

    def test_get_revision_string(self):
        result = svn.get_revision_string(codepath.get_module_path(svn))
        self.assertTrue(isinstance(result, str))

    def test_parse_revision_string(self):
        self.assertRaises(svn.RepositoryStateException, svn.parse_revision_string, "224:232")
        self.assertRaises(svn.RepositoryStateException, svn.parse_revision_string, "224:232M")
        self.assertRaises(svn.RepositoryStateException, svn.parse_revision_string, "")
        self.assertRaises(svn.RepositoryStateException, svn.parse_revision_string, None)
        self.assertEqual(100789, svn.parse_revision_string("100789"))

    def test_get_revision_id_for_number(self):
        self.assertRaises(ValueError, svn.get_revision_id_for_number, "224:232")
        self.assertEqual("R000123", svn.get_revision_id_for_number(123))
