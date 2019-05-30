"""
Test path retrieval.
"""

# pylint: disable=missing-docstring

import unittest
import os
import eumopps.version.codepath
import eumopps.version.test


class TestCodePath(unittest.TestCase):

    def test_get_module_path(self):
        """Check the path retrieval method as best we can - can only do relative paths."""

        # Path at base of this package
        basepath = os.path.normpath(os.path.join(eumopps.version.codepath.get_module_path(eumopps), '..'))

        # Our path should be like [basepath]/eumopps/version/test
        testpath = eumopps.version.codepath.get_module_path(eumopps.version.test)

        # Check the same
        self.assertEqual(os.path.join(basepath, 'eumopps', 'version', 'test'), testpath)
