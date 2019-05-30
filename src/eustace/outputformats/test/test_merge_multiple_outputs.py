"""Testing the implementation of the generalized merging"""

import importlib
import os
from StringIO import StringIO
import sys
import unittest
from ..merge_multiple_outputs import return_daily_outputs
from ..merge_multiple_outputs import merge_single_day_outputs
from netCDF4 import Dataset
from datetime import datetime
from eumopps.timeutils.timebase import TimeBaseDays
from eustace.surfaceairmodel.correlation_ranges import CorrelationRanges
from eustace.surfaceairmodel.fileio.ocean import OCEAN_CORRELATED_UNCERTAINTIES
from eustace.surfaceairmodel.fileio.ocean import OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS
from eustace.surfaceairmodel.fileio.ocean import OCEAN_CORRELATED_UNCERTAINTIES_TIMES
from eustace.surfaceairmodel.fileio.ocean import UNITS
from shutil import copyfile

class TestMergeMultipleOutputs(unittest.TestCase):

    def test_return_daily_outputs(self):

	input_dictionary = {'catalogue_id': 'NO-CATALOGUE', 'institution': 'MetOffice', 
			    'output_main':'baba_main', 'output_ancillary':'baba_ancillary', 
                            'correlation_indexes':[.1,.2], 'comment_prefix':'Gargantua is the rule', 'definitions_module': 'eustace.outputformats.merge_multiple_outputs',
                            'primary_fields_list':['I','you'], 'uncertainty_fields_list':['notsureI','notsureyou'], 'ancillary_fields_list':['maybehime?'],
                            'correlation_ranges':'a python object',
			    'input_main':'a file', 'input_ancillary':'another file'}

	self.assertDictEqual(input_dictionary,return_daily_outputs(**input_dictionary))

    def test_merge_single_day_outputs(self):
	""" Testing correct formats of merged output files: using the ocean case (the only one with necessary merging)"""

	# preparing input parameters for the ocean case

	source_labels=['AATSR','ATSR2']
	inputs = {'catalogue_id' : 'NO_CATALOGUE',
		  'institution' : 'Met Office',
                  'output_main': ['some_main_file.nc'],
		  'output_ancillary':['some_ancillary_file.nc'],
                  'correlation_indexes':'[[0.],[1.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]]',
                  'comment_prefix':'Merging from ',
                  'definitions_module':'eustace.outputformats.definitions',
		  'primary_fields_list':['TAS'],
                  'uncertainty_fields_list':['TASUNCERTAINTY'],
                  'ancillary_fields_list': ['SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_RANDOM',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED2',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC2',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER0',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER1',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER2',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER3',
					    'SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER4'],
		 'correlation_ranges' : CorrelationRanges(OCEAN_CORRELATED_UNCERTAINTIES, UNITS, OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS, OCEAN_CORRELATED_UNCERTAINTIES_TIMES),
		 'input_main' : [None,None],
		 'input_ancillary' : [None,None]}

	# build output variable objects for testing purposes

	module = importlib.import_module(inputs['definitions_module'])
	primary_fields_list = [getattr(module, name) for name in inputs['primary_fields_list']]
	primary_uncertainties_list = [getattr(module, name) for name in inputs['uncertainty_fields_list']]
	ancillary_fields_list = [getattr(module, name) for name in inputs['ancillary_fields_list']]

	# 1st test case: no input to merge: raise an exception

	old_stdout = sys.stdout
	result = StringIO()
	sys.stdout = result
	merge_single_day_outputs(**inputs)  
	#self.assertEquals(sys.stdout.getvalue().strip(),'eustace.outputformats.merge_single_day_outputs WARNING:No input files to merge')
        self.assertEquals(sys.stdout.getvalue().splitlines()[3],'eustace.outputformats.merge_single_day_outputs WARNING:No input files to merge')


	# Redirect again the std output to screen
 
	sys.stdout = old_stdout

	# 2n case: we have some output to merge

	basepath = '/gws/nopw/j04/eustace/data/internal/satgrid_lst/test_0.25_0.25/'
	outputs_main=['tas_ocean_eustace_0_AATSR_20020714.nc','tas_ocean_eustace_0_ATSR2_20020714.nc']
	outputs_ancillary=['ancillary_ocean_eustace_0_AATSR_20020714.nc','ancillary_ocean_eustace_0_ATSR2_20020714.nc']
	source_labels=['AATSR','ATSR2']

	suffix='copy.'

	outputs_main_copy = []
	outputs_ancillary_copy = []
	
	# outputs from different sources gets removed after the merging: we need to copy them to preserve original test data

	for name_main, name_ancillary in zip(outputs_main,outputs_ancillary):
	  oldname_main = basepath+name_main
	  oldname_ancillary = basepath+name_ancillary
	  newname_main = basepath+suffix+name_main
	  newname_ancillary = basepath+suffix+name_ancillary
	  copyfile(oldname_main, newname_main)
	  copyfile(oldname_ancillary, newname_ancillary)
	  outputs_main_copy.append(newname_main)
	  outputs_ancillary_copy.append(newname_ancillary)

	copied_main = [basepath+suffix+outputs_main[0],basepath+suffix+outputs_main[1]]
	copied_ancillary = [basepath+suffix+outputs_ancillary[0],basepath+suffix+outputs_ancillary[1]]
	final_main = [basepath+'copy.tas_ocean_eustace_0_20020714.nc']
	final_ancillary = [basepath+'copy.ancillary_ocean_eustace_0_20020714.nc']

	inputs['output_main'] = final_main
	inputs['output_ancillary'] = final_ancillary
	inputs['input_main'] = copied_main
	inputs['input_ancillary'] = copied_ancillary

	merge_single_day_outputs(**inputs)
	
	# Test the result of merging: time of creation, details, variable names, correlation time/scales values and units
	dataset = Dataset(final_main[0],'r')

	date = TimeBaseDays(datetime(1850,01,01)).datetime_to_number(datetime(2002,07,14))
	self.assertEqual(dataset.variables['time'][...][0],date)
	self.assertEqual(dataset.institution,'Met Office')
	self.assertEqual(dataset.comment,'Merging from '+source_labels[0]+', '+source_labels[1]+' retrievals.')
	self.assertEqual(dataset.source,'EUSTACE Catalogue NO_CATALOGUE')
	for variable in primary_fields_list+primary_uncertainties_list:
	  self.assertTrue(dataset.variables.has_key(variable.name))
	dataset.close()

	dataset = Dataset(final_ancillary[0],'r')

	self.assertEqual(dataset.variables['time'][...][0],date)
	self.assertEqual(dataset.institution,'Met Office')
	self.assertEqual(dataset.comment,'Merging from '+source_labels[0]+', '+source_labels[1]+' retrievals.')
	self.assertEqual(dataset.source,'EUSTACE Catalogue NO_CATALOGUE')
	for variable in primary_fields_list+primary_uncertainties_list:
	    for source in source_labels:
	      self.assertTrue(dataset.variables.has_key(variable.name+'_'+source))
	for variable in ancillary_fields_list:
	    self.assertTrue(dataset.variables.has_key(variable.name))  
	    for source in source_labels:
	      self.assertTrue(dataset.variables.has_key(variable.name+'_'+source))
	    for key, length_scale, time_scale in zip(OCEAN_CORRELATED_UNCERTAINTIES, OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS, OCEAN_CORRELATED_UNCERTAINTIES_TIMES):
	      self.assertEqual(dataset.variables[key].length_scale,length_scale,"Testing correlation length value for "+key)
	      self.assertEqual(dataset.variables[key].length_scale_units,UNITS[0],"Testing correlation length units value for "+key)
	      self.assertEqual(dataset.variables[key].time_scale,time_scale,"Testing correlation time value for "+key)
	      self.assertEqual(dataset.variables[key].time_scale_units,UNITS[1],"Testing correlation time units value for "+key)
              for source in source_labels:
		tempkey = key+'_'+source
		self.assertEqual(dataset.variables[tempkey].length_scale,length_scale,"Testing correlation length value for "+tempkey)
		self.assertEqual(dataset.variables[tempkey].length_scale_units,UNITS[0],"Testing correlation length units value for "+tempkey)
		self.assertEqual(dataset.variables[tempkey].time_scale,time_scale,"Testing correlation time value for "+tempkey)
		self.assertEqual(dataset.variables[tempkey].time_scale_units,UNITS[1],"Testing correlation time units value for "+tempkey)

	dataset.close()
	for name in final_main+final_ancillary:
	  os.remove(name)
	  
	# 3d case: we have just one output to merge -> we need just renaming and update of correlation lengths infos

	basepath = '/gws/nopw/j04/eustace/data/internal/satgrid_lst/test_0.25_0.25/'
	outputs_main=['tas_ocean_eustace_0_AATSR_20030623.nc']
	outputs_ancillary=['ancillary_ocean_eustace_0_AATSR_20030623.nc']
	source_labels=['AATSR']

	suffix='copy.'

	outputs_main_copy = []
	outputs_ancillary_copy = []
	
	# outputs from different sources gets removed after the merging: we need to copy them to preserve original test data

	for name_main, name_ancillary in zip(outputs_main,outputs_ancillary):
	  oldname_main = basepath+name_main
	  oldname_ancillary = basepath+name_ancillary
	  newname_main = basepath+suffix+name_main
	  newname_ancillary = basepath+suffix+name_ancillary
	  copyfile(oldname_main, newname_main)
	  copyfile(oldname_ancillary, newname_ancillary)
	  outputs_main_copy.append(newname_main)
	  outputs_ancillary_copy.append(newname_ancillary)

	copied_main = [basepath+suffix+outputs_main[0],None]
	copied_ancillary = [basepath+suffix+outputs_ancillary[0],None]
	final_main = [basepath+'copy.tas_ocean_eustace_0_20030623.nc']
	final_ancillary = [basepath+'copy.ancillary_ocean_eustace_0_20030623.nc']

	inputs['output_main'] = final_main
	inputs['output_ancillary'] = final_ancillary
	inputs['input_main'] = copied_main
	inputs['input_ancillary'] = copied_ancillary

	merge_single_day_outputs(**inputs)
	
	# Test the result of merging: time of creation, details, variable names, correlation time/scales values and units
	dataset = Dataset(final_main[0],'r')

	date = TimeBaseDays(datetime(1850,01,01)).datetime_to_number(datetime(2003,06,23))
	self.assertEqual(dataset.variables['time'][...][0],date)
	self.assertEqual(dataset.institution,'Met Office')
	self.assertEqual(dataset.comment,'Merging from '+source_labels[0]+' retrievals.')
	for variable in primary_fields_list+primary_uncertainties_list:
	  self.assertTrue(dataset.variables.has_key(variable.name))
	dataset.close()

	dataset = Dataset(final_ancillary[0],'r')

	self.assertEqual(dataset.variables['time'][...][0],date)
	self.assertEqual(dataset.institution,'Met Office')
	self.assertEqual(dataset.comment,'Merging from '+source_labels[0]+' retrievals.')
	for variable in ancillary_fields_list:
	    self.assertTrue(dataset.variables.has_key(variable.name))
	    for key, length_scale, time_scale in zip(OCEAN_CORRELATED_UNCERTAINTIES, OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS, OCEAN_CORRELATED_UNCERTAINTIES_TIMES):
	      self.assertEqual(dataset.variables[key].length_scale,length_scale,"Testing correlation length value for "+key)
	      self.assertEqual(dataset.variables[key].length_scale_units,UNITS[0],"Testing correlation length units value for "+key)
	      self.assertEqual(dataset.variables[key].time_scale,time_scale,"Testing correlation time value for "+key)
	      self.assertEqual(dataset.variables[key].time_scale_units,UNITS[1],"Testing correlation time units value for "+key)

	dataset.close()

	# Cleaning
	
	for name in final_main+final_ancillary:
	  os.remove(name)
