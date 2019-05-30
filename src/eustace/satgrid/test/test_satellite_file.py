"""Tests for satellite_file module."""
# pylint: disable=missing-docstring,bad-whitespace

__version__ = "$Revision: 472 $"
__author__ = "Joel R. Mitchelson"

import unittest
import os
import numpy
import eustaceconfig
from ..satellite_file import SatelliteFilename
from ..satellite_file import SatelliteFile


class TestSatelliteFile(unittest.TestCase):

    def test_load(self):
        filename_lst = os.path.join(
            eustaceconfig.WORKSPACE_PATH, 'data/incoming/MODIS/2010/03/04/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-20100304_120500-CUOL-0.01X0.01-V2.0.nc')
        filename_aux = os.path.join(
            eustaceconfig.WORKSPACE_PATH, 'data/incoming/MODIS/2010/03/04/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-20100304_120500-CUOL-0.01X0.01-V2.0.nc')
        satfile = SatelliteFile(SatelliteFilename(filename_lst, filename_aux))
        self.assertTrue(satfile.handle_lst)
        self.assertTrue(satfile.handle_aux)

    def test_get_variables(self):

        # load it
        filename_lst = os.path.join(
            eustaceconfig.WORKSPACE_PATH, 'data/incoming/MODIS/2010/03/04/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-20100304_060500-CUOL-0.01X0.01-V2.0.nc')
        filename_aux = os.path.join(
            eustaceconfig.WORKSPACE_PATH, 'data/incoming/MODIS/2010/03/04/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-20100304_060500-CUOL-0.01X0.01-V2.0.nc')
        satfile = SatelliteFile(SatelliteFilename(filename_lst, filename_aux))

        # check LST
        lst = satfile.get_variable('LST')
        self.assertTrue(isinstance(lst, numpy.ndarray))
        self.assertEqual(lst.dtype, numpy.float32)
        # expected shape based on ncdump
        self.assertEqual((2030, 1354), lst.shape)

        # check QC
        qc = satfile.get_variable('QC') # pylint: disable=invalid-name
        self.assertTrue(isinstance(qc, numpy.ndarray))
        self.assertEqual(qc.dtype, numpy.int32)
        # expected shape based on ncdump
        self.assertEqual((2030, 1354), qc.shape)

        # check lat
        lat = satfile.get_variable('lat')
        self.assertTrue(isinstance(lat, numpy.ndarray))
        self.assertEqual(lat.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lat.shape)

        # check lon
        lon = satfile.get_variable('lon')
        self.assertTrue(isinstance(lon, numpy.ndarray))
        self.assertEqual(lon.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lon.shape)

        # check LST_unc_ran (random uncertainty)
        lst_unc_ran = satfile.get_variable('LST_unc_ran')
        self.assertTrue(isinstance(lst_unc_ran, numpy.ndarray))
        self.assertEqual(lst_unc_ran.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lst_unc_ran.shape)

        # check LST_unc_loc (locally correlated uncertainty, atmospheric component)
        lst_unc_loc_atm = satfile.get_variable('LST_unc_loc_atm')
        self.assertTrue(isinstance(lst_unc_loc_atm, numpy.ndarray))
        self.assertEqual(lst_unc_loc_atm.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lst_unc_loc_atm.shape)

        # check LST_unc_loc (locally correlated uncertainty, surface component)
        lst_unc_loc_sfc = satfile.get_variable('LST_unc_loc_sfc')
        self.assertTrue(isinstance(lst_unc_loc_sfc, numpy.ndarray))
        self.assertEqual(lst_unc_loc_sfc.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lst_unc_loc_sfc.shape)

        # check LST_unc_sys (systematic uncertainty: just one number per file)
        lst_unc_sys = satfile.get_variable('LST_unc_sys')
        self.assertTrue(isinstance(lst_unc_sys, numpy.ndarray))
        self.assertEqual(lst_unc_sys.dtype, numpy.float32)
        self.assertEqual((1, ), lst_unc_sys.shape)

    def test_find_qc_flag_indices(self):
        testdata = numpy.array([[0, 1, 3, 7, 2],
                                [1, 0, 8, 3, 3],
                                [7, 2, 6, 3, 1]], numpy.int32)
        file_a = SatelliteFile.find_qc_flag_indices(testdata, 0xFFFF, 0x0008)
        self.assertEqual([7], file_a.tolist())
        file_b = SatelliteFile.find_qc_flag_indices(testdata, 0x0004, 0x0004)
        self.assertEqual([3, 10, 12], file_b.tolist())
        file_c = SatelliteFile.find_qc_flag_indices(testdata, 0x0003, 0x0002)
        self.assertEqual([4, 11, 12], file_c.tolist())
