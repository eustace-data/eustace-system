"""Test of run method."""

import unittest
import tempfile
import os
from netCDF4 import Dataset
from ..run_ocean import run_day
from ..run_ocean import return_daily_inputs
from ..run_ocean import run_multiple_days
from eustace.surfaceairmodel.fileio.ocean import UNITS
from eustace.surfaceairmodel.fileio.ocean import OCEAN_CORRELATED_UNCERTAINTIES
from eustace.surfaceairmodel.fileio.ocean import OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS
from eustace.surfaceairmodel.fileio.ocean import OCEAN_CORRELATED_UNCERTAINTIES_TIMES
import eustaceconfig
from comparison.comparefields_ocean import comparefields_ocean

class TestRunOcean(unittest.TestCase):

    def test_run_day(self):
	"""Test single day processing"""

        # details
        catalogue_id = 'NOCATALOGUE'
        institution = 'TEST'
	
	# outputs
        outputs = [ tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.test.test_run_ocean', suffix='.nc') for index in range(2) ]
        output_main = outputs[0].name
        output_ancillary = outputs[1].name

	# model parameters
        input_model_parameters = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_v.5.0.dat')
        input_model_uncertainty = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_uncertainty_v.5.0.dat')
        input_model_stdev = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_STDEV_clim_parameters_v.5.0.dat')

        input_satellite = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/UoR_sst/l3c/AATSR/2005/03/200503211200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc')

        run_day(catalogue_id, institution, output_main, output_ancillary, input_model_parameters, input_model_uncertainty, input_model_stdev, input_satellite)

        # do full comparison
        issues = comparefields_ocean(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/MS8_version2/Ocean/data/v.5.0/2005/eustace_AATSR_satellite_ocean_20050321.nc'),
                                      outputs[0].name, outputs[1].name)
        # should have no issues
        self.assertEqual(0, len(issues))

        # testing correct value of correlation ranges
	dataset = Dataset(outputs[1].name, 'r')
	for j,k in enumerate(OCEAN_CORRELATED_UNCERTAINTIES):
	  self.assertEqual(OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS[j],dataset.variables[k].length_scale)
	  self.assertEqual(OCEAN_CORRELATED_UNCERTAINTIES_TIMES[j],dataset.variables[k].time_scale)
	  self.assertEqual(UNITS[0],dataset.variables[k].length_scale_units)
	  self.assertEqual(UNITS[1],dataset.variables[k].time_scale_units)
	dataset.close()

    def test_return_daily_inputs(self):
	
	# details
        catalogue_id = 'NOCATALOGUE'
        institution = 'TEST'

	# outputs
        output_main = 'Brombobo'
        output_ancillary = 'Brombobo_ancillary'


	# model parameters
        input_model_parameters = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_v.5.0.dat')
        input_model_uncertainty = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_uncertainty_v.5.0.dat')
        input_model_stdev = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_STDEV_clim_parameters_v.5.0.dat')

        input_satellite_ocean = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/UoR_sst/l3c/AATSR/2005/03/200503211200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc')

	dictionary = return_daily_inputs(catalogue_id, institution, output_main, output_ancillary, input_model_parameters, input_model_uncertainty, input_model_stdev, input_satellite_ocean)

	# check details
	self.assertEqual(dictionary['details']['catalogue_id'], catalogue_id)
	self.assertEqual(dictionary['details']['institution'], institution)

	# check input parameters
	self.assertItemsEqual(dictionary['input_parameters']['input_model_parameters'],input_model_parameters)
	self.assertEqual(dictionary['input_parameters']['input_model_uncertainty'], input_model_uncertainty)
	self.assertEqual(dictionary['input_parameters']['input_model_stdev'], input_model_stdev)

	# check input files
	self.assertEqual(dictionary['input_files']['input_satellite_ocean'], input_satellite_ocean)
	
	# check output files
	self.assertEqual(dictionary['output_files']['output_main'], output_main)
	self.assertEqual(dictionary['output_files']['output_ancillary'], output_ancillary)

    def test_run_multiple_days(self):

        # details
        catalogue_id = 'NOCATALOGUE'
        institution = 'TEST'

        # ouputs
        outputs = [ tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.test.test_run_ocean', suffix='.nc') for index in range(4) ]
        outputs_main = [ outputs[0].name, outputs[2].name ]
        outputs_ancillary = [ outputs[1].name, outputs[3].name ]

        # configuration
        input_model_parameters = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_v.5.0.dat')
        input_model_uncertainty = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_uncertainty_v.5.0.dat')
        input_model_stdev = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ocean/hires_AST_STDEV_clim_parameters_v.5.0.dat')
        inputs_satellite_ocean = [
            os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/UoR_sst/l3c/AATSR/2005/03/200503211200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc'),
            os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/UoR_sst/l3c/AATSR/2007/01/200701011200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc') ]

	# creating list  of daily inputs
	list_of_daily_inputs = []
	for index in range(len(outputs_main)):
	    daily_inputs = return_daily_inputs(catalogue_id, institution, outputs_main[index], outputs_ancillary[index], input_model_parameters, input_model_uncertainty, input_model_stdev,inputs_satellite_ocean[index])
	    list_of_daily_inputs.append(daily_inputs)


        # run
        run_multiple_days(list_of_daily_inputs)

        # do full comparison
        issuesA = comparefields_ocean(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/MS8_version2/Ocean/data/v.5.0/2005/eustace_AATSR_satellite_ocean_20050321.nc'),
                                      outputs[0].name, outputs[1].name)
        issuesB = comparefields_ocean(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/MS8_version2/Ocean/data/v.5.0/2007/eustace_AATSR_satellite_ocean_20070101.nc'),
                                      outputs[2].name, outputs[3].name)

        # should have no issues
        self.assertEqual(0, len(issuesA))
        self.assertEqual(0, len(issuesB))

        
        # testing correct value of correlation ranges
	for i in [outputs[1].name,outputs[3].name]:
	  dataset = Dataset(i, 'r')
	  for j,k in enumerate(OCEAN_CORRELATED_UNCERTAINTIES):
	    self.assertEqual(OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS[j],dataset.variables[k].length_scale)
	    self.assertEqual(OCEAN_CORRELATED_UNCERTAINTIES_TIMES[j],dataset.variables[k].time_scale)
	    self.assertEqual(UNITS[0],dataset.variables[k].length_scale_units)
	    self.assertEqual(UNITS[1],dataset.variables[k].time_scale_units)
          dataset.close()
