"""Tests for SatelliteFileNamesMODIS
   and also a quick check that MODIS files can be loaded to Iris cubes."""

import unittest
from ..modis import SatelliteFieldNamesMODIS
from eustaceconfig import WORKSPACE_PATH
import os
import numpy
from iris import load_cubes
from iris import load_cube
from iris import load_raw
from iris import Constraint
from iris.coord_systems import GeogCS
from iris.coords import AuxCoord

class TestSatelliteFieldNamesMODIS(unittest.TestCase):

    def test_init(self):

        fields = SatelliteFieldNamesMODIS()
        self.assertEqual(fields.primary, 'LST')
        self.assertEqual(fields.qc, 'QC')
        self.assertEqual(3, len(fields.uncertainty_fields))
        self.assertEqual(1, len(fields.uncertainty_scalars))
        self.assertTrue('LST_unc_ran' in fields.uncertainty_fields)
        self.assertTrue('LST_unc_loc_atm' in fields.uncertainty_fields)
        self.assertTrue('LST_unc_loc_sfc' in fields.uncertainty_fields)
        self.assertTrue('LST_unc_sys' in fields.uncertainty_scalars)

    def test_MODIS_format_iris_compatibility(self):

        # Filenames
        filename_lst=os.path.join(WORKSPACE_PATH, 'data/incoming/MODIS/2010/01/01/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-20100101_120000-CUOL-0.01X0.01-V2.0.nc')
        filename_aux=os.path.join(WORKSPACE_PATH, 'data/incoming/MODIS/2010/01/01/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-20100101_120000-CUOL-0.01X0.01-V2.0.nc')

        # Load them
        cube_LST = load_cube(filename_lst, 'surface_temperature')
        cube_LST_unc_ran = load_cube(filename_aux, constraint=Constraint(cube_func=lambda cube: cube.var_name == 'LST_unc_ran'))
        cube_LST_unc_loc_atm = load_cube(filename_aux, constraint=Constraint(cube_func=lambda cube: cube.var_name == 'LST_unc_loc_atm'))
        cube_LST_unc_loc_sfc = load_cube(filename_aux, constraint=Constraint(cube_func=lambda cube: cube.var_name == 'LST_unc_loc_sfc'))
        cube_LST.add_aux_coord(AuxCoord(points=cube_LST_unc_ran.data, var_name='LST_unc_ran'), (0,1))
        cube_LST.add_aux_coord(AuxCoord(points=cube_LST_unc_loc_atm.data, var_name='LST_unc_loc_atm'), (0,1))
        cube_LST.add_aux_coord(AuxCoord(points=cube_LST_unc_loc_sfc.data, var_name='LST_unc_loc_sfc'), (0,1))

        # Generate summary text
        summarytext = str(cube_LST)

        # This is how it should begin
        expected_summary = 'surface_temperature / (K)           (-- : 2030; -- : 1354)\n' + \
                           '     Auxiliary coordinates:\n' + \
                           '          LST_unc_loc_atm               x          x\n' + \
                           '          LST_unc_loc_sfc               x          x\n' + \
                           '          LST_unc_ran                   x          x\n' + \
                           '          latitude                      x          x\n' + \
                           '          longitude                     x          x\n'

        # Check it matches
        self.assertTrue(summarytext.startswith(expected_summary))

