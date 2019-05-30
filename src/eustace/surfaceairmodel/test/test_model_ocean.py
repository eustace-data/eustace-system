from ..model_ocean import ModelOcean
import unittest
import tempfile
from datetime import datetime
import netCDF4
import json

class TestModelOcean(unittest.TestCase):

    def test_read_netcdf_time(self):
        
        result = ModelOcean.read_netcdf_time('/gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/AATSR/2007/02/200702231200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc')
        self.assertEqual(datetime(2007, 2, 23, 11, 0), result)

    def test_build_idl_commandline(self):

        result = ModelOcean.build_idl_commandline('bob.txt')
        self.assertEqual(5, len(result))
        self.assertEqual('idl', result[0])
        self.assertEqual('-e', result[1])
        self.assertTrue(result[2].endswith('model_ocean_idl.pro'))
        self.assertEqual('-args', result[3])
        self.assertEqual('bob.txt', result[4])

    def test_build_config_dictionary(self):

        model = ModelOcean('a', 'b', 'c')
        details = model.build_config_dictionary([ 'd' ], [ '/gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/AATSR/2007/02/200702231200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc' ])
        self.assertEqual('a', details['model']['filename_parameters'])
        self.assertEqual('b', details['model']['filename_uncertainty'])
        self.assertEqual('c', details['model']['filename_stdev'])
        self.assertEqual(1, len(details['data']))
        self.assertEqual(2007, details['data'][0]['date']['year'])
        self.assertEqual(   2, details['data'][0]['date']['month'])
        self.assertEqual(  23, details['data'][0]['date']['day'])
        self.assertEqual('/gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/AATSR/2007/02/200702231200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc', details['data'][0]['filename_input'])
        self.assertEqual('d', details['data'][0]['filename_output'])

    def test_write_config_file(self):

        # Temporary file to write to
        example = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.test.test_model_ocean', suffix='.json')

        # Make model and write
        model = ModelOcean('a', 'b', 'c')
        model.write_config_file(example.name, [ 'd' ], [ '/gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/AATSR/2007/02/200702231200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc' ])

        # Read JSON and check results
        with open(example.name, 'r') as examplefile:
            details = json.load(examplefile)
        self.assertEqual('a', details['model']['filename_parameters'])
        self.assertEqual('b', details['model']['filename_uncertainty'])
        self.assertEqual('c', details['model']['filename_stdev'])
        self.assertEqual(1, len(details['data']))
        self.assertEqual(2007, details['data'][0]['date']['year'])
        self.assertEqual(   2, details['data'][0]['date']['month'])
        self.assertEqual(  23, details['data'][0]['date']['day'])
        self.assertEqual('/gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/AATSR/2007/02/200702231200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc', details['data'][0]['filename_input'])
        self.assertEqual('d', details['data'][0]['filename_output'])

    def test_run_idl(self):

        outputfile = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.test.TestModelOcean_', suffix='.nc')

        model = ModelOcean('/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_v.5.0.dat',
                           '/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ocean/hires_AST_clim_parameters_uncertainty_v.5.0.dat',
                           '/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ocean/hires_AST_STDEV_clim_parameters_v.5.0.dat')        

        model.run_idl( [ outputfile.name ], [ '/gws/nopw/j04/eustace/data/incoming/UoR_sst/l3c/AATSR/2007/01/200701011200-ESACCI-L3C_GHRSST-AATSR-EXP1.2-v02.0-fv1.0.nc' ] )        

        # check netcdf output produced something
        dataset = netCDF4.Dataset(outputfile.name)
        tas = dataset.variables['tas']
        self.assertEqual('air_temperature', tas.standard_name)
