"""Test post-processing module"""
from eustace.outputformats.definitions import TASMIN, TASMAX, TASMINUNCERTAINTY, TASMAXUNCERTAINTY
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMIN
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMAX
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_TMINMODELNUMBER
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_TMAXMODELNUMBER
from ..post_processing import apply_land_post_processing
import numpy
import numpy.ma
import numpy.testing 
import netCDF4
import os
import tempfile
import unittest

TEST_MASKS = [numpy.array([[[False, True,  False],[True, False, False],[False, False, False],[False, False, True]]]),
              numpy.array([[[False, False, False],[False, False, False],[False, False, False],[True, True, False]]])]

TEST_FUNDAMENTAL_FIELDS = [numpy.array([[[3. , -111.,  0.],[-111.,  0., 10.],[4., 18., 4.],[   7.,    7., -111.]]]),
          	           numpy.array([[[23.,    0., 10.],[  77., 30.,  2.],[1., 22., 3.],[-111., -111.,    1.]]])]

TEST_FUNDAMENTAL_UNCERTAINTIES = [numpy.array([[[.3,  -111.,  .1],[ -111., .9,  .1],[.4,   .18,   .4],[  .7,     .7, -111.]]]),
                                   numpy.array([[[.3,     .2,  .2],[    .7, .3, .22],[.11,  .22,  .66],[-111., -111.,    1.]]])]

TEST_ANCILLARY_1 = [numpy.array([[[.1,  -111.,  .1],[ -111., .9,  .1],[.4,   .18,   .4],[  .1,     .7, -111.]]]),
                    numpy.array([[[.2,     .4,  .4],[    .7, .1, .22],[.11,  .22,   .6],[-111., -111.,    1.]]])]

TEST_ANCILLARY_2 = [numpy.array([[[0,  -111.,  0],[ -111., 0,  0],[0,   2,  0],[1,        1, -111.]]]),
                    numpy.array([[[0,     1,  .0],[    0,  2,  2],[1,   2,  0],[-111., -111.,    0]]])]


EXPECTED_MASKS = [numpy.array([[[False, True,  False],[True, False, True],[True, False, True],[False, False, True]]]),
                  numpy.array([[[False, False, False],[False, False, True],[True, False, True],[True, True, False]]])]

EXPECTED_FUNDAMENTAL_FIELDS = [numpy.array([[[3. , -111.,  0.],[-111.,  0., -111.],[-111., 18., -111.],[   7.,    7., -111.]]]),
          	               numpy.array([[[23.,    0., 10.],[  77., 30., -111.],[-111., 22., -111.],[-111., -111.,    1.]]])]

EXPECTED_FUNDAMENTAL_UNCERTAINTIES = [numpy.array([[[.3,  -111.,  .1],[ -111., .9, -111.],[-111., .18, -111.],[  .7,     .7, -111.]]]),
                                       numpy.array([[[.3,     .2,  .2],[    .7, .3, -111.],[-111., .22, -111.],[-111., -111.,    1.]]])]

EXPECTED_ANCILLARY_1 = [numpy.array([[[.1,  -111.,  .1],[ -111., .9,  -111.],[-111.,  .18,   -111.],[  .1,     .7, -111.]]]),
                        numpy.array([[[.2,     .4,  .4],[    .7, .1,  -111.],[-111.,  .22,   -111.],[-111., -111.,    1.]]])]

EXPECTED_ANCILLARY_2 = [numpy.array([[[0,  -111.,  0],[ -111., 0,  -111.],[-111.,   2,  -111.],[1,        1, -111.]]]),
                        numpy.array([[[0,     1,  .0],[    0,  2,  -111.],[-111.,   2,  -111.],[-111., -111.,    0]]])]

class TestPostProcessing(unittest.TestCase):
  def test_apply_land_post_processing(self):
      """Test the correctness of post-processing of satstace land output"""
    
      results = {TASMIN.name: numpy.ma.MaskedArray(data = TEST_FUNDAMENTAL_FIELDS[0], mask = TEST_MASKS[0], fill_value = -111.),
		 TASMAX.name: numpy.ma.MaskedArray(data = TEST_FUNDAMENTAL_FIELDS[1], mask = TEST_MASKS[1], fill_value = -111.),
		 TASMINUNCERTAINTY.name: numpy.ma.MaskedArray(data = TEST_FUNDAMENTAL_UNCERTAINTIES[0], mask = TEST_MASKS[0], fill_value = -111.),
		 TASMAXUNCERTAINTY.name: numpy.ma.MaskedArray(data = TEST_FUNDAMENTAL_UNCERTAINTIES[1], mask = TEST_MASKS[1], fill_value = -111.),
                 SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMIN.name: numpy.ma.MaskedArray(data = TEST_ANCILLARY_1[0], mask = TEST_MASKS[0], fill_value = -111.),
                 SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMAX.name: numpy.ma.MaskedArray(data = TEST_ANCILLARY_1[1], mask = TEST_MASKS[1], fill_value = -111.),
		 SURFACEAIRMODEL_ANCILLARY_LAND_TMINMODELNUMBER.name: numpy.ma.MaskedArray(data = TEST_ANCILLARY_2[0], mask = TEST_MASKS[0], fill_value = -111.),
		 SURFACEAIRMODEL_ANCILLARY_LAND_TMAXMODELNUMBER.name: numpy.ma.MaskedArray(data = TEST_ANCILLARY_2[1], mask = TEST_MASKS[1], fill_value = -111.)}

      apply_land_post_processing(results, [SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMIN, SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMAX, SURFACEAIRMODEL_ANCILLARY_LAND_TMINMODELNUMBER, SURFACEAIRMODEL_ANCILLARY_LAND_TMAXMODELNUMBER])
    
      for index, variables in enumerate(zip([TASMIN, TASMAX], [TASMINUNCERTAINTY, TASMAXUNCERTAINTY], [SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMIN, SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMAX], [SURFACEAIRMODEL_ANCILLARY_LAND_TMINMODELNUMBER, SURFACEAIRMODEL_ANCILLARY_LAND_TMAXMODELNUMBER])):
	numpy.testing.assert_array_equal(results[variables[0].name][:].mask, EXPECTED_MASKS[index])
	numpy.testing.assert_array_equal(results[variables[0].name][:].data, EXPECTED_FUNDAMENTAL_FIELDS[index])
	numpy.testing.assert_array_equal(results[variables[1].name][:].mask, EXPECTED_MASKS[index])
	numpy.testing.assert_array_equal(results[variables[1].name][:].data, EXPECTED_FUNDAMENTAL_UNCERTAINTIES[index])
	numpy.testing.assert_array_equal(results[variables[2].name][:].mask, EXPECTED_MASKS[index])
	numpy.testing.assert_array_equal(results[variables[2].name][:].data, EXPECTED_ANCILLARY_1[index])
	numpy.testing.assert_array_equal(results[variables[3].name][:].mask, EXPECTED_MASKS[index])
	numpy.testing.assert_array_equal(results[variables[3].name][:].data, EXPECTED_ANCILLARY_2[index])