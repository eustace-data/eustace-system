"""Integration tests for entire system."""

import unittest
import os
import subprocess
from StringIO import StringIO
import eumopps.catalogue.build
import eumopps.catalogue.commandcount
import eumopps.catalogue.commandrun
from eumopps.catalogue.fileio.test.test_formatnetcdf import NamedTemporaryDirectory

class TestSystem(unittest.TestCase):

    def test_daily_and_fixed(self):
        """Run end-to-end test.
           Demonstrates batch operation using example files in data/daily and data/onefile as inputs.
           Runs examplemodule on each daily input using onefile as fixed input."""

        # Change to test directory as file locations are relative to this
        testdirectory = os.path.dirname(__file__)
        if testdirectory is not '':            
            os.chdir(testdirectory)
        

        # Output dir
        outputdirectory = NamedTemporaryDirectory()

        # Temporary catalogue file to use in output directory
        cataloguefilename = os.path.join(outputdirectory.pathname, 'test_catalogue.nc')

        # Generate catalogue
        eumopps.catalogue.build.main(['--allow_unversioned_code', 'test_system_inputs.json', cataloguefilename ])

        # Update with operations
        eumopps.catalogue.build.main(['--allow_unversioned_code', '--pathdefault', outputdirectory.pathname, '--update', 'test_system_operation.json', cataloguefilename])

        # TEMP: inspect output
        # print '\n' + subprocess.Popen(['ncdump', cataloguefilename], stdout=subprocess.PIPE).communicate()[0]

        # Get operation count
        commandcount_string = StringIO()
        eumopps.catalogue.commandcount.main(argvalues=[cataloguefilename, 'eumopps.catalogue.test.examplemodule.run'], outputstream=commandcount_string)
        commandcount = int(commandcount_string.getvalue())

        # Should be 5 of them (one per day)
        self.assertEqual(5, commandcount)

        # Execute commands
        for commandindex in range(commandcount):
            eumopps.catalogue.commandrun.main(['--allow_unversioned_code', cataloguefilename, 'eumopps.catalogue.test.examplemodule.run', str(commandindex)])

        # Read result file contents
        result0 = open(os.path.join(outputdirectory.pathname, 'example_20170418_result.txt'), 'r').read()
        result1 = open(os.path.join(outputdirectory.pathname, 'example_20170419_result.txt'), 'r').read()
        result2 = open(os.path.join(outputdirectory.pathname, 'example_20170420_result.txt'), 'r').read()
        result3 = open(os.path.join(outputdirectory.pathname, 'example_20170421_result.txt'), 'r').read()
        result4 = open(os.path.join(outputdirectory.pathname, 'example_20170422_result.txt'), 'r').read()

        # Check values
        self.assertEqual('My Results\n' +
                         '20170418\n' +
                         'one\n' +
                         'Move along now, nothing to see here.\n',
                         result0)
        self.assertEqual('My Results\n' +
                         '20170419\n' +
                         'two\n' +
                         'Move along now, nothing to see here.\n',
                         result1)
        self.assertEqual('My Results\n' +
                         '20170420\n' +
                         'three\n' +
                         'Move along now, nothing to see here.\n',
                         result2)
        self.assertEqual('My Results\n' +
                         '20170421\n' +
                         'four\n' +
                         'Move along now, nothing to see here.\n',
                         result3)
        self.assertEqual('My Results\n' +
                         '20170422\n' +
                         'five\n' +
                         'Move along now, nothing to see here.\n',
                         result4)


    def test_daily_input_list(self):
        """Run end-to-end test.
           Demonstrates batch operation using example files in data/twosources and data/onefile as inputs.
           Runs examplemodule_inputlist using list of two inputs."""

        # Change to test directory as file locations are relative to this
        testdirectory = os.path.dirname(__file__)
        os.chdir(testdirectory)

        # Output dir
        outputdirectory = NamedTemporaryDirectory()

        # Temporary catalogue file to use in output directory
        cataloguefilename = os.path.join(outputdirectory.pathname, 'test_catalogue.nc')

        # Generate catalogue
        eumopps.catalogue.build.main(['--allow_unversioned_code', '--pathdefault', outputdirectory.pathname, 'test_system_input_list.json', cataloguefilename ])

        # TEMP: inspect output
        # print '\n' + subprocess.Popen(['ncdump', cataloguefilename], stdout=subprocess.PIPE).communicate()[0]

        # Get operation count
        commandcount_string = StringIO()
        eumopps.catalogue.commandcount.main(argvalues=[cataloguefilename, 'eumopps.catalogue.test.examplemodule_input_list.run'], outputstream=commandcount_string)
        commandcount = int(commandcount_string.getvalue())

        # Should be 6 of them (one per day)
        self.assertEqual(6, commandcount)

        # Execute commands
        for commandindex in range(commandcount):
            eumopps.catalogue.commandrun.main(['--allow_unversioned_code', cataloguefilename, 'eumopps.catalogue.test.examplemodule_input_list.run', str(commandindex)])

        # Read result file contents
        result0 = open(os.path.join(outputdirectory.pathname, 'example_20171210_result.txt'), 'r').read()
        result1 = open(os.path.join(outputdirectory.pathname, 'example_20171211_result.txt'), 'r').read()
        result2 = open(os.path.join(outputdirectory.pathname, 'example_20171212_result.txt'), 'r').read()
        result3 = open(os.path.join(outputdirectory.pathname, 'example_20171213_result.txt'), 'r').read()
        result4 = open(os.path.join(outputdirectory.pathname, 'example_20171214_result.txt'), 'r').read()
        result5 = open(os.path.join(outputdirectory.pathname, 'example_20171215_result.txt'), 'r').read()

        # Check values
        self.assertEqual('My Results\n' +
                         '20171210\n' +
                         '--\n' +
                         'sourcetwo: ten\n' +
                         'Move along now, nothing to see here.\n',
                         result0)
        self.assertEqual('My Results\n' +
                         '20171211\n' +
                         'sourceone: eleven\n' +
                         '--\n' +
                         'Move along now, nothing to see here.\n',
                         result1)
        self.assertEqual('My Results\n' +
                         '20171212\n' +
                         'sourceone: twelve\n' +
                         '--\n' +
                         'Move along now, nothing to see here.\n',
                         result2)
        self.assertEqual('My Results\n' +
                         '20171213\n' +
                         'sourceone: thirteen\n' +
                         'sourcetwo: thirteen\n' +
                         'Move along now, nothing to see here.\n',
                         result3)
        self.assertEqual('My Results\n' +
                         '20171214\n' +
                         '--\n' +
                         '--\n' +
                         'Move along now, nothing to see here.\n',
                         result4)
        self.assertEqual('My Results\n' +
                         '20171215\n' +
                         '--\n' +
                         'sourcetwo: fifteen\n' +
                         'Move along now, nothing to see here.\n',
                         result5)

