from ..model_land import ModelLand
from eumopps.timeutils.timebase import TimeBaseDays
from datetime import datetime
import json
import netCDF4
import os
import tempfile
import unittest

EPOCH=datetime(1850,1,1)

# Path to files need to be updated in the future
FILENAMES_INPUT=[['/gws/nopw/j04/eustace/users/ejn2/Output/L3/FVC/with_unc/global/res_0.25_0.25/',
                          '/gws/nopw/j04/eustace/users/ejn2/Output/L3/DEM/with_unc/global/DEM_global_0.25_0.25.nc',
                          ['/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_nh.nc',
                           '/gws/nopw/j04/eustace/data/incoming/DMI_Masks/mask_land_water_icecap_fitosisaf_sh.nc'],
                          '/gws/nopw/j04/eustace/users/ejn2/Output/L3/snow/global/res_0.25_0.25/snow_global_0.25_0.25_20110104.nc',
                          ['/gws/nopw/j04/eustace/data/internal/satgrid_lst_update/res_0.25/satgrid.day.20110104.nc',
                           '/gws/nopw/j04/eustace/data/internal/satgrid_lst_update/res_0.25/satgrid.night.20110104.nc']]]
FILENAMES_OUTPUT=[['./','eustace.surfaceairmodel.test.TestModelLand_20110104']]

class TestModelLand(unittest.TestCase):

    def test_read_input_file_time(self):
        
        result = ModelLand.read_input_file_time('/gws/nopw/j04/eustace/data/internal/satgrid_lst/res_0.25/satgrid.day.20041101.nc')
        self.assertEqual(datetime(2004, 11, 01), result)

    def test_build_idl_commandline(self):

        result = ModelLand.build_idl_commandline('bob.txt')
        self.assertEqual(5, len(result))
        self.assertEqual('idl', result[0])
        self.assertEqual('-e', result[1])
        self.assertTrue(result[2].endswith('model_land_idl.pro'))
        self.assertEqual('-args', result[3])
        self.assertEqual('bob.txt', result[4])

    def test_build_config_dictionary(self):

        model = ModelLand(['a','b','c'], 1.43, .46)

        details = model.build_config_dictionary(FILENAMES_OUTPUT, FILENAMES_INPUT)
        self.assertEqual(['a','b','c'], details['model']['filename_model_clim_parameters'])
        self.assertEqual(1.43, details['model']['frac_sat_obs_threshold'])
        self.assertEqual(.46, details['model']['sampling_threshold'])
        self.assertEqual(1, len(details['data']))
        self.assertEqual(2011, details['data'][0]['date']['year'])
        self.assertEqual(   1, details['data'][0]['date']['month'])
        self.assertEqual(  4, details['data'][0]['date']['day'])
        self.assertEqual(FILENAMES_INPUT[0][0], details['data'][0]['dir_input_fvc'])
        self.assertEqual(FILENAMES_INPUT[0][1], details['data'][0]['filename_input_dem'])
        self.assertEqual(FILENAMES_INPUT[0][2][0], details['data'][0]['filename_input_mask'][0])
        self.assertEqual(FILENAMES_INPUT[0][2][1], details['data'][0]['filename_input_mask'][1])
        self.assertEqual(FILENAMES_INPUT[0][3], details['data'][0]['filename_input_snow'])
        self.assertEqual(FILENAMES_INPUT[0][4][0], details['data'][0]['filename_input_lst'][0])        
        self.assertEqual(FILENAMES_INPUT[0][4][1], details['data'][0]['filename_input_lst'][1])   
        self.assertEqual(FILENAMES_OUTPUT[0][1], details['data'][0]['filename_output_id'])
        self.assertEqual(FILENAMES_OUTPUT[0][0], details['data'][0]['dir_output'])

    def test_write_config_file(self):

        # Temporary file to write to
        example = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.test.test_model_land', suffix='.json')

        # Make model and write

        model = ModelLand(['a','b','c'], 1.43, .46)
        
        model.write_config_file(example.name, FILENAMES_OUTPUT, FILENAMES_INPUT)

        # Read JSON and check results
        with open(example.name, 'r') as examplefile:
            details = json.load(examplefile)
	    self.assertEqual(['a','b','c'], details['model']['filename_model_clim_parameters'])
	    self.assertEqual(1.43, details['model']['frac_sat_obs_threshold'])
	    self.assertEqual(.46, details['model']['sampling_threshold'])
	    self.assertEqual(1, len(details['data']))
	    self.assertEqual(2011, details['data'][0]['date']['year'])
	    self.assertEqual(   1, details['data'][0]['date']['month'])
	    self.assertEqual(  4, details['data'][0]['date']['day'])
	    self.assertEqual(FILENAMES_INPUT[0][0], details['data'][0]['dir_input_fvc'])
	    self.assertEqual(FILENAMES_INPUT[0][1], details['data'][0]['filename_input_dem'])
            self.assertEqual(FILENAMES_INPUT[0][2][0], details['data'][0]['filename_input_mask'][0])
            self.assertEqual(FILENAMES_INPUT[0][2][1], details['data'][0]['filename_input_mask'][1])
	    self.assertEqual(FILENAMES_INPUT[0][3], details['data'][0]['filename_input_snow'])
	    self.assertEqual(FILENAMES_INPUT[0][4][0], details['data'][0]['filename_input_lst'][0])        
	    self.assertEqual(FILENAMES_INPUT[0][4][1], details['data'][0]['filename_input_lst'][1])   
	    self.assertEqual(FILENAMES_OUTPUT[0][1], details['data'][0]['filename_output_id'])
	    self.assertEqual(FILENAMES_OUTPUT[0][0], details['data'][0]['dir_output'])

    def test_run_idl(self):

        outputfile = tempfile.NamedTemporaryFile(prefix=FILENAMES_OUTPUT[0][1]+'_', suffix='.nc')
        
        # Should be updated to final directory for parameters
        model = ModelLand(['/gws/nopw/j04/eustace/users/ejn2/Output/sat_lsat/May2017/MYG_model1_coefs.txt',
                            '/gws/nopw/j04/eustace/users/ejn2/Output/sat_lsat/May2017/MYG_model2_coefs.txt',
                            '/gws/nopw/j04/eustace/users/ejn2/Output/sat_lsat/May2017/MYG_model3_coefs.txt'],
                           .2,3.)
        
        model.run_idl( [[tempfile.gettempdir()+'/',os.path.basename(outputfile.name)]], FILENAMES_INPUT )
        
        #check netcdf output produced something
        dataset = netCDF4.Dataset(outputfile.name)

        lat=dataset.variables['lat']
        self.assertEqual('latitude', lat.standard_name)
        lat=dataset.variables['lon']
        self.assertEqual('longitude', lat.standard_name)
        tasmin=dataset.variables['tasmin']
        self.assertEqual('air_temperature', tasmin.standard_name)
        self.assertEqual('Minimum daily surface air temperature (K)', tasmin.long_name)
        tasmax=dataset.variables['tasmax']
        self.assertEqual('air_temperature', tasmax.standard_name)
        self.assertEqual('Maximum daily surface air temperature (K)', tasmax.long_name)
        unc_rand_tasmin=dataset.variables['unc_rand_tasmin']
        self.assertEqual('Random uncertainty on minimum daily surface air temperature (K)', unc_rand_tasmin.long_name)
        unc_rand_tasmax=dataset.variables['unc_rand_tasmax']
        self.assertEqual('Random uncertainty on maximum daily surface air temperature (K)', unc_rand_tasmax.long_name)
        unc_corr_atm_tasmin=dataset.variables['unc_corr_atm_tasmin']
        self.assertEqual('Locally correlated uncertainty (atmospheric scale lengths) on minimum daily surface air temperature (K)', unc_corr_atm_tasmin.long_name)
	unc_corr_sfc_tasmin=dataset.variables['unc_corr_sfc_tasmin']
        self.assertEqual('Locally correlated uncertainty (surface scale lengths) on minimum daily surface air temperature (K)', unc_corr_sfc_tasmin.long_name)
        unc_corr_atm_tasmax=dataset.variables['unc_corr_atm_tasmax']
        self.assertEqual('Locally correlated uncertainty (atmospheric scale lengths) on maximum daily surface air temperature (K)', unc_corr_atm_tasmax.long_name)
	unc_corr_sfc_tasmax=dataset.variables['unc_corr_sfc_tasmax']
        self.assertEqual('Locally correlated uncertainty (surface scale lengths) on maximum daily surface air temperature (K)', unc_corr_sfc_tasmax.long_name)
        unc_sys_tasmin=dataset.variables['unc_sys_tasmin']
        self.assertEqual('Systematic uncertainty on minimum daily surface air temperature (K)', unc_sys_tasmin.long_name)
        unc_sys_tasmax=dataset.variables['unc_sys_tasmax']
        self.assertEqual('Systematic uncertainty on maximum daily surface air temperature (K)', unc_sys_tasmax.long_name)
        tasmin_model_number=dataset.variables['tasmin_model_number']
        self.assertEqual('Model number used for estimating Tmin from satellite data', tasmin_model_number.long_name)
        tasmax_model_number=dataset.variables['tasmax_model_number']
        self.assertEqual('Model number used for estimating Tmax from satellite data', tasmax_model_number.long_name)
        
        file_time=dataset.variables['time']
        self.assertEqual('time', file_time.standard_name)
        self.assertEqual('Time (days)', file_time.long_name)
        t=TimeBaseDays(EPOCH)
        file_date=t.number_to_datetime(int(file_time[0]))
        self.assertEqual(datetime.strptime('20110104','%Y%m%d'),file_date)
