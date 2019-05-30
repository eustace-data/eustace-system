"""Tests of code which generates example formats."""

__version__ = "$Revision: 824 $"
__author__ = "Joel R. Mitchelson"

import unittest
from eustace.outputformats.examples import example_global_field
from eustace.outputformats.examples import example_station_data
from eustace.outputformats.verify.verify_globalfield import VerifyGlobalField
from eustace.outputformats.verify.verify_stations import VerifyStationData
import tempfile
import StringIO
import os
import shutil
import sys

class VerifierArguments(object):
    """Helper class to provide namespace for verification arguments."""

    def __init__(self):
        pass


class TestExampleGlobalField(unittest.TestCase):
    """Test global field example format."""

    def test_write_read_verify(self):
        """Write an example format and verify against specification document."""

        # make temporary directory to hold results
        directory = tempfile.mkdtemp(prefix='eustace_test_examples_')

        # use try-finally to remove directory even if test fails
        try:

            example_global_field.write_example_global_field(
                'some_source_file.py',
                'R090807',
                directory,
                institution='Test Institution Name')

            # pylint: disable=attribute-defined-outside-init
            args = VerifierArguments()
            args.pathname = os.path.join(directory, '2015/tas_example_eustace_0_20151105.nc')
            outputstream = StringIO.StringIO()
            # outputstream = sys.stderr # use stderr if diagnostics needed
            self.assertTrue(VerifyGlobalField().run(args, outputstream))

        finally:

            # remove temporary directory
            shutil.rmtree(directory)


class TestExampleStationData(unittest.TestCase):
    """Test station data example format."""

    def test_write_read_verify(self):
        """Write an example format and verify against specification document."""

        # make temporary directory to hold results
        directory = tempfile.mkdtemp(prefix='eustace_test_examples_')

        # use try-finally to remove directory even if test fails
        try:

            example_station_data.write_example_station_data(
                'some_source_file.py',
                'R090807',
                directory,
                institution='Test Institution Name')

            # pylint: disable=attribute-defined-outside-init
            args = VerifierArguments()
            args.pathnameprefix = os.path.join(directory, 'eustace_stations_example_R090807')
            outputstream = StringIO.StringIO()
            # outputstream = sys.stderr # use stderr if diagnostics needed
            self.assertTrue(VerifyStationData().run(args, outputstream))

        finally:

            # remove temporary directory
            shutil.rmtree(directory)
