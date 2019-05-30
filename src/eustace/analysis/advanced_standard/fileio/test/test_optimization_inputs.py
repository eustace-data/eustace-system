"""Integration tests for entire system."""

import unittest
import os
import subprocess
from StringIO import StringIO
import eumopps.catalogue.build
import eumopps.catalogue.commandcount
import eumopps.catalogue.commandrun
from eumopps.catalogue.fileio.test.test_formatnetcdf import NamedTemporaryDirectory
import dateutil.parser

class TestSystem(unittest.TestCase):

    def test_summaries_create(self):
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
        eumopps.catalogue.build.main(['--allow_unversioned_code', '--pathdefault', outputdirectory.pathname, '--update', 'test_optimization_input_list.json', cataloguefilename])

        # Get operation count
        commandcount_string = StringIO()
        eumopps.catalogue.commandcount.main(argvalues=[cataloguefilename, 'input_summaries_create'], outputstream=commandcount_string)
        commandcount = int(commandcount_string.getvalue())
        
        # Should be 5 of them (one per input time)
        self.assertEqual(5, commandcount)

        # Execute commands
        for commandindex in range(commandcount):
            eumopps.catalogue.commandrun.main(['--allow_unversioned_code', cataloguefilename, 'input_summaries_create', str(commandindex)])


        from  eumopps.catalogue.fileio import formatjson
        result0 = formatjson.CatalogueReaderJSON().load(os.path.join(outputdirectory.pathname, 'example_20170418_inputs.json'))
        result1 = formatjson.CatalogueReaderJSON().load(os.path.join(outputdirectory.pathname, 'example_20170419_inputs.json'))
        result2 = formatjson.CatalogueReaderJSON().load(os.path.join(outputdirectory.pathname, 'example_20170420_inputs.json'))
        result3 = formatjson.CatalogueReaderJSON().load(os.path.join(outputdirectory.pathname, 'example_20170421_inputs.json'))
        result4 = formatjson.CatalogueReaderJSON().load(os.path.join(outputdirectory.pathname, 'example_20170422_inputs.json'))
        
        self.assertEqual( "20170418", dateutil.parser.parse(result0[0]).strftime("%Y%m%d") )
        self.assertEqual( 'data/daily/exampledaily_20170418.txt', result0[1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', result0[1][0]['inputfixed'] )
        
        self.assertEqual( "20170419", dateutil.parser.parse(result1[0]).strftime("%Y%m%d") )
        self.assertEqual( 'data/daily/exampledaily_20170419.txt', result1[1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', result1[1][0]['inputfixed'] )
        
        self.assertEqual( "20170420", dateutil.parser.parse(result2[0]).strftime("%Y%m%d") )
        self.assertEqual( 'data/daily/exampledaily_20170420.txt', result2[1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', result2[1][0]['inputfixed'] )
        
        self.assertEqual( "20170421", dateutil.parser.parse(result3[0]).strftime("%Y%m%d") )
        self.assertEqual( 'data/daily/exampledaily_20170421.txt', result3[1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', result3[1][0]['inputfixed'] )
        
        self.assertEqual( "20170422", dateutil.parser.parse(result4[0]).strftime("%Y%m%d") )
        self.assertEqual( 'data/daily/exampledaily_20170422.txt', result4[1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', result4[1][0]['inputfixed'] )
        
        # need to test reading the full data set

    def test_input_list(self):
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
        eumopps.catalogue.build.main(['--allow_unversioned_code', '--pathdefault', outputdirectory.pathname, '--update', 'test_optimization_input_list.json', cataloguefilename])

        # Get operation count
        commandcount_string = StringIO()
        eumopps.catalogue.commandcount.main(argvalues=[cataloguefilename, 'input_summaries_create'], outputstream=commandcount_string)
        commandcount = int(commandcount_string.getvalue())
        
        # Should be 5 of them (one per input time)
        self.assertEqual(5, commandcount)

        # Execute commands
        for commandindex in range(commandcount):
            eumopps.catalogue.commandrun.main(['--allow_unversioned_code', cataloguefilename, 'input_summaries_create', str(commandindex)])

        # Get operation count
        commandcount_string = StringIO()
        eumopps.catalogue.commandcount.main(argvalues=[cataloguefilename, 'input_summary_merge'], outputstream=commandcount_string)
        commandcount = int(commandcount_string.getvalue())
        # Should be 1 of them (one per output file)
        self.assertEqual(1, commandcount)
        eumopps.catalogue.commandrun.main([cataloguefilename, 'input_summary_merge', 0, '--allow_unversioned_code']) 
        
        from  eumopps.catalogue.fileio import formatjson    
        input_list = formatjson.CatalogueReaderJSON().load(os.path.join(outputdirectory.pathname, 'merged_input_summary.json'))
        
        self.assertEqual( "20170418", dateutil.parser.parse(input_list[0][0]).strftime("%Y%m%d") )
        self.assertEqual( 'data/daily/exampledaily_20170418.txt', input_list[0][1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', input_list[0][1][0]['inputfixed'] )
        
        self.assertEqual( "20170419", dateutil.parser.parse(input_list[1][0]).strftime("%Y%m%d" ) )
        self.assertEqual( 'data/daily/exampledaily_20170419.txt', input_list[1][1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', input_list[1][1][0]['inputfixed'] )
        
        self.assertEqual( "20170420", dateutil.parser.parse(input_list[2][0]).strftime("%Y%m%d" ) )
        self.assertEqual( 'data/daily/exampledaily_20170420.txt', input_list[2][1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', input_list[2][1][0]['inputfixed'] )
        
        self.assertEqual( "20170421", dateutil.parser.parse(input_list[3][0]).strftime("%Y%m%d" ) )
        self.assertEqual( 'data/daily/exampledaily_20170421.txt', input_list[3][1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', input_list[3][1][0]['inputfixed'] )
        
        self.assertEqual( "20170422", dateutil.parser.parse(input_list[4][0]).strftime("%Y%m%d" ) )
        self.assertEqual( 'data/daily/exampledaily_20170422.txt', input_list[4][1][0]['inputdaily'] )
        self.assertEqual( 'data/onefile/thisfile.txt', input_list[4][1][0]['inputfixed'] )
        