"""Test of run method."""

import unittest
import tempfile
import os
import netCDF4
from ..run_land import run_day
from ..run_land import run_multiple_days
from ..run_land import return_daily_inputs
from eustace.surfaceairmodel.fileio.land import UNITS
from eustace.surfaceairmodel.fileio.land import LAND_CORRELATED_UNCERTAINTIES
from eustace.surfaceairmodel.fileio.land import LAND_CORRELATED_UNCERTAINTIES_LENGTHS
from eustace.surfaceairmodel.fileio.land import LAND_CORRELATED_UNCERTAINTIES_TIMES
import eustaceconfig
from comparison.comparefields_land import comparefields_land

class TestRunLand(unittest.TestCase):

    def test_run_day(self):

        # details
        catalogue_id = 'NOCATALOGUE'
        institution = 'TEST'

	# temporary output
        outputs = [ tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.test.test_run_land', suffix='.nc') for index in range(2) ]
        output_main = outputs[0].name
        output_ancillary = outputs[1].name

	# model parameters
        input_model_1 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model1_coefs.txt')
	input_model_2 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model2_coefs.txt')
        input_model_3 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model3_coefs.txt')
	
	# input from satgrid_iris
	input_satellite_land_day = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.day.20120101.nc')
	input_satellite_land_night = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.night.20120101.nc')

	# snow
        input_satellite_land_snow ='/gws/nopw/j04/eustace/users/ejn2/Output/L3/snow/global/res_0.25_0.25/snow_global_0.25_0.25_20120101.nc'

	# DMI mask
        input_satellite_land_mask_north = '/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_nh.nc'
        input_satellite_land_mask_south = '/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_fitosisaf_sh.nc'

	# DEM
	input_satellite_land_DEM = '/gws/nopw/j04/eustace/users/ejn2/Output/L3/DEM/with_unc/global/DEM_global_0.25_0.25.nc'

	# FVC
        input_satellite_land_fvc = '/gws/nopw/j04/eustace/users/ejn2/Output/L3/FVC/with_unc/global/res_0.25_0.25/'

	# Fixed parameters
        frac_sat_obs_threshold = .2
        sampling_threshold = 3.
        
	run_day(catalogue_id, institution, 
		output_main, output_ancillary, 
		input_model_1, input_model_2, input_model_3,
		input_satellite_land_day, input_satellite_land_night, 
		input_satellite_land_snow,
		input_satellite_land_mask_north, input_satellite_land_mask_south,
		input_satellite_land_DEM,
		input_satellite_land_fvc,
		frac_sat_obs_threshold, 
		sampling_threshold)

        issues = comparefields_land(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/test_0.25_0.25/eustace_satellite_4.100_20120101.nc'),
				    outputs[0].name, outputs[1].name)
        # should have no issues
        self.assertEqual(0, len(issues))

        # testing correct value of correlation ranges
	dataset = netCDF4.Dataset(outputs[1].name, 'r')
	for j,k in enumerate(LAND_CORRELATED_UNCERTAINTIES):
	    self.assertEqual(LAND_CORRELATED_UNCERTAINTIES_LENGTHS[j],dataset.variables[k].length_scale)
	    self.assertEqual(LAND_CORRELATED_UNCERTAINTIES_TIMES[j],dataset.variables[k].time_scale)
	    self.assertEqual(UNITS[0],dataset.variables[k].length_scale_units)
	    self.assertEqual(UNITS[1],dataset.variables[k].time_scale_units)
	dataset.close()

    def test_return_daily_inputs(self):
	
	# details
        catalogue_id = 'NOCATALOGUE'
        institution = 'TEST'

        output_main = 'Brombobo'
        output_ancillary = 'Brombobo_ancillary'

	# model parameters
        input_model_1 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model1_coefs.txt')
	input_model_2 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model2_coefs.txt')
        input_model_3 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model3_coefs.txt')
	
	# input from satgrid_iris
	input_satellite_land_day = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.day.20120101.nc')
	input_satellite_land_night = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.night.20120101.nc')

	# snow
        input_satellite_land_snow ='/gws/nopw/j04/eustace/users/ejn2/Output/L3/snow/global/res_0.25_0.25/snow_global_0.25_0.25_20120101.nc'

	# DMI mask
        input_satellite_land_mask_north = '/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_nh.nc'
        input_satellite_land_mask_south = '/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_fitosisaf_sh.nc'

	# DEM
	input_satellite_land_DEM = '/gws/nopw/j04/eustace/users/ejn2/Output/L3/DEM/with_unc/global/DEM_global_0.25_0.25.nc'

	# FVC
        input_satellite_land_fvc = '/gws/nopw/j04/eustace/users/ejn2/Output/L3/FVC/with_unc/global/res_0.25_0.25/'

	# Fixed parameters
        frac_sat_obs_threshold = .2
        sampling_threshold = 3.

	dictionary = return_daily_inputs(catalogue_id,institution,
			output_main,output_ancillary,
			input_model_1,input_model_2,input_model_3,
			input_satellite_land_day,input_satellite_land_night,
			input_satellite_land_snow,
			input_satellite_land_mask_north,input_satellite_land_mask_south,
			input_satellite_land_DEM,
			input_satellite_land_fvc,
			frac_sat_obs_threshold, 
			sampling_threshold)

	# check details
	self.assertEqual(dictionary['details']['catalogue_id'], catalogue_id)
	self.assertEqual(dictionary['details']['institution'], institution)

	# check input parameters
	self.assertItemsEqual(dictionary['input_parameters']['input_model_parameters'],[input_model_1,input_model_2,input_model_3])
	self.assertEqual(dictionary['input_parameters']['frac_sat_obs_threshold'], frac_sat_obs_threshold)
	self.assertEqual(dictionary['input_parameters']['sampling_threshold'], sampling_threshold)

	# check input files
	self.assertEqual(dictionary['input_files']['input_fvc'], input_satellite_land_fvc)
	self.assertEqual(dictionary['input_files']['input_DEM'], input_satellite_land_DEM)
	self.assertEqual(dictionary['input_files']['input_snow'], input_satellite_land_snow)
	self.assertItemsEqual(dictionary['input_files']['input_mask'], [input_satellite_land_mask_north,input_satellite_land_mask_south])
	self.assertItemsEqual(dictionary['input_files']['inputs_day'], [input_satellite_land_day,input_satellite_land_night])

	# check output files
	self.assertEqual(dictionary['output_files']['output_main'], output_main)
	self.assertEqual(dictionary['output_files']['output_ancillary'], output_ancillary)


    def test_run_multiple_days(self):

        # details
        catalogue_id = 'NOCATALOGUE'
        institution = 'TEST'

	# outputs
	outputs = [ tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.test.test_run_land', suffix='.nc') for index in range(6) ]
        outputs_main = [ outputs[0].name, outputs[2].name, outputs[4].name]
        outputs_ancillary = [ outputs[1].name, outputs[3].name, outputs[5].name]

	# model parameters
        input_model_1 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model1_coefs.txt')
        input_model_2 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model2_coefs.txt')
	input_model_3 = os.path.join(eustaceconfig.WORKSPACE_PATH, 'users/ejn2/Output/sat_lsat/May2017/MYG_model3_coefs.txt')
        frac_sat_obs_threshold = .2
        sampling_threshold = 3.
	
	# FVC
        input_fvc = '/gws/nopw/j04/eustace/users/ejn2/Output/L3/FVC/with_unc/global/res_0.25_0.25/'

	# DEM
        input_DEM = '/gws/nopw/j04/eustace/users/ejn2/Output/L3/DEM/with_unc/global/DEM_global_0.25_0.25.nc'

	# DMI mask
        input_mask = ['/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_nh.nc',
                      '/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_fitosisaf_sh.nc']
        
        inputs_snow=['/gws/nopw/j04/eustace/users/ejn2/Output/L3/snow/global/res_0.25_0.25/snow_global_0.25_0.25_20120101.nc',
                     '/gws/nopw/j04/eustace/users/ejn2/Output/L3/snow/global/res_0.25_0.25/snow_global_0.25_0.25_20070101.nc',
                     '/gws/nopw/j04/eustace/users/ejn2/Output/L3/snow/global/res_0.25_0.25/snow_global_0.25_0.25_20130101.nc']
        
	inputs_day = [[os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.day.20120101.nc'),
		       os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.night.20120101.nc')],
		      [os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.day.20070101.nc'),
                       os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.night.20070101.nc')],
                      [os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.day.20130101.nc'),
                       os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/res_0.25/satgrid.night.20130101.nc')]]        

	# creating list  of daily inputs
	list_of_daily_inputs = []
	for index in range(len(outputs_main)):
	    daily_inputs = return_daily_inputs(catalogue_id,institution,
						outputs_main[index],outputs_ancillary[index],
						input_model_1,input_model_2,input_model_3,
						inputs_day[index][0],inputs_day[index][1],
						inputs_snow[index],
						input_mask[0],input_mask[1],
						input_DEM,
						input_fvc,
						frac_sat_obs_threshold, 
						sampling_threshold)
	    list_of_daily_inputs.append(daily_inputs)

        # run
        run_multiple_days(list_of_daily_inputs)
        # do full comparison
        issuesA = comparefields_land(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/test_0.25_0.25/eustace_satellite_4.100_20120101.nc'),
                                     outputs[0].name, outputs[1].name)
        issuesB = comparefields_land(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/test_0.25_0.25/eustace_satellite_4.100_20070101.nc'),
                                      outputs[2].name, outputs[3].name)
        issuesC = comparefields_land(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/satgrid_lst_update/test_0.25_0.25/eustace_satellite_4.100_20130101.nc'),
                                      outputs[4].name, outputs[5].name)

        # should have no issues
        self.assertEqual(0, len(issuesA))
        self.assertEqual(0, len(issuesB))
        self.assertEqual(0, len(issuesC))

        # testing correct value of correlation ranges
	for i in [outputs[1].name,outputs[3].name,outputs[5].name]:
	  dataset = netCDF4.Dataset(i, 'r')
	  for j,k in enumerate(LAND_CORRELATED_UNCERTAINTIES):

	    self.assertEqual(LAND_CORRELATED_UNCERTAINTIES_LENGTHS[j],dataset.variables[k].length_scale)
	    self.assertEqual(LAND_CORRELATED_UNCERTAINTIES_TIMES[j],dataset.variables[k].time_scale)
	    self.assertEqual(UNITS[0],dataset.variables[k].length_scale_units)
	    self.assertEqual(UNITS[1],dataset.variables[k].time_scale_units)
	  dataset.close()