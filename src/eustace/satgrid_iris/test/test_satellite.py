"""Tests for satellite_file module."""
# pylint: disable=missing-docstring,bad-whitespace

__version__ = "$Revision: 978 $"
__author__ = "Joel R. Mitchelson"

import unittest
import os
import numpy
import eustaceconfig
from datetime import datetime
from iris.cube import Cube
from iris.cube import CubeList
from iris.coords import Coord
from iris.coords import DimCoord
from iris.coords import AuxCoord
from ..satellite import SatelliteInputList
from ..satellite import SatelliteCollection
from ..satellite import SatelliteFieldNames
from ..satellite import ERRORVALUE
from ..modis import SatelliteFieldNamesMODIS

class TestSatelliteInputList(unittest.TestCase):

    class PathPatternDate(object):
        """Data class for file name patterns."""

        def __init__(self, path, pattern_lst, pattern_aux, date):
            self.path = path
            self.pattern_lst = pattern_lst
            self.pattern_aux = pattern_aux
            self.date = date

    def test_build_one_filename(self):
        result = SatelliteInputList(TestSatelliteInputList.PathPatternDate('bob', 'afile.%Y-%m-%d.%H_%M.nc', 'auxstuff.%Y%m%d.%H%M.nc', '193311031233'))
        self.assertEqual('single', result.time_mode)
        self.assertEqual(1, len(result.filenames))
        self.assertEqual('bob/afile.1933-11-03.12_33.nc', result.filenames[0][0])
        self.assertEqual('bob/auxstuff.19331103.1233.nc', result.filenames[0][1])

    def test_build_wholeday_filenames(self):

        result = SatelliteInputList(TestSatelliteInputList.PathPatternDate('who', 'me.%Y%m%d.%H%M.nc', 'yes.%Y%m%d.%H%M.nc', '20010723'))
        self.assertEqual('day', result.time_mode)
        self.assertEqual(288, len(result.filenames))
        self.assertEqual('who/me.20010723.0000.nc', result.filenames[0][0])
        self.assertEqual('who/yes.20010723.0000.nc', result.filenames[0][1])
        self.assertEqual('who/me.20010723.0005.nc', result.filenames[1][0])
        self.assertEqual('who/yes.20010723.0005.nc', result.filenames[1][1])
        self.assertEqual('who/me.20010723.0110.nc', result.filenames[14][0])
        self.assertEqual('who/yes.20010723.0110.nc', result.filenames[14][1])
        self.assertEqual('who/me.20010723.2350.nc', result.filenames[286][0])
        self.assertEqual('who/yes.20010723.2350.nc', result.filenames[286][1])
        self.assertEqual('who/me.20010723.2355.nc', result.filenames[287][0])
        self.assertEqual('who/yes.20010723.2355.nc', result.filenames[287][1])

        # causes problems if relying on exception thrown by
        # datetime.strptime(args.date, '%Y%m%d%H%M')
        # to detect single file / whole day mode
        result = SatelliteInputList(TestSatelliteInputList.PathPatternDate('a', 'b.%Y%m%d.%H%M.nc', 'c.%Y%m%d.%H%M.nc', '20101102'))
        self.assertEqual('day', result.time_mode)
        self.assertEqual(288, len(result.filenames))
        self.assertEqual('a/b.20101102.0000.nc', result.filenames[0][0])
        self.assertEqual('a/c.20101102.0000.nc', result.filenames[0][1])
        self.assertEqual('a/b.20101102.2355.nc', result.filenames[287][0])
        self.assertEqual('a/c.20101102.2355.nc', result.filenames[287][1])

class TestSatelliteCollection(unittest.TestCase):

    # Field names used in example to test QC filter
    class FieldNamesForTest(SatelliteFieldNames):

        def __init__(self):

            super(TestSatelliteCollection.FieldNamesForTest, self).__init__(
                primary='EXAMPLE',
                qc='QC',
                coordinate_fields=['latitude','longitude'],
                uncertainty_fields=[ ],
                uncertainty_scalars=[ ],
                output_primary='exampleoutput',
                output_uncertainty=[ ],
                output_uncertainty_binned_scalars=[ ])

    def test_load_modis_example(self):

        # example filenames for testing
        filename_lst = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/MODIS/2010/03/04/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-20100304_060500-CUOL-0.01X0.01-V2.0.nc')
        filename_aux = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/MODIS/2010/03/04/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-20100304_060500-CUOL-0.01X0.01-V2.0.nc')

        # load
        result = SatelliteCollection.load(SatelliteFieldNamesMODIS(), [filename_lst, filename_aux])

        # check properties are populated appropriately
        self.assertEqual(result.fields.__dict__, SatelliteFieldNamesMODIS().__dict__)
        self.assertTrue(isinstance(result.cube, Cube))
	self.assertFalse(result.safety_flag)
	self.assertTrue(isinstance(result.safety_mask,numpy.ndarray))
	self.assertEqual(result.safety_mask.shape,(1,))
	self.assertFalse(result.safety_mask[0])

        # get cube to check individual variables
        cube = result.cube

        # check main variable
        lst = cube.data
        self.assertTrue(isinstance(lst, numpy.ndarray))
        self.assertEqual(lst.dtype, numpy.float32)
        # expected shape based on ncdump
        self.assertEqual((2030, 1354), lst.shape)

        # check QC
        qc = cube.coords('QC')[0].points  # pylint: disable=invalid-name
        self.assertTrue(isinstance(qc, numpy.ndarray))
        self.assertEqual(qc.dtype, numpy.int32)
        # expected shape based on ncdump
        self.assertEqual((2030, 1354), qc.shape)

        # check lat
        lat = cube.coords('latitude')[0].points
        self.assertTrue(isinstance(lat, numpy.ndarray))
        self.assertEqual(lat.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lat.shape)

        # check lon
        lon = cube.coords('longitude')[0].points
        self.assertTrue(isinstance(lon, numpy.ndarray))
        self.assertEqual(lon.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lon.shape)

        # check LST_unc_ran (random uncertainty)
        lst_unc_ran = cube.coords('LST_unc_ran')[0].points
        self.assertTrue(isinstance(lst_unc_ran, numpy.ndarray))
        self.assertEqual(lst_unc_ran.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lst_unc_ran.shape)

        # check LST_unc_loc_atm (locally correlated uncertainty, atm component)
        lst_unc_loc_atm = cube.coords('LST_unc_loc_atm')[0].points
        self.assertTrue(isinstance(lst_unc_loc_atm, numpy.ndarray))
        self.assertEqual(lst_unc_loc_atm.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lst_unc_loc_atm.shape)

        # check LST_unc_loc_atm (locally correlated uncertainty, atm component)
        lst_unc_loc_sfc = cube.coords('LST_unc_loc_sfc')[0].points
        self.assertTrue(isinstance(lst_unc_loc_sfc, numpy.ndarray))
        self.assertEqual(lst_unc_loc_sfc.dtype, numpy.float32)
        self.assertEqual((2030, 1354), lst_unc_loc_sfc.shape)

        # check LST_unc_sys (systematic uncertainty: just one number per file)
        lst_unc_sys = cube.coords('LST_unc_sys')[0].points
        self.assertTrue(isinstance(lst_unc_sys, numpy.ndarray))
        self.assertEqual(lst_unc_sys.dtype, numpy.float32)
        self.assertEqual((1, ), lst_unc_sys.shape)

    def test_apply_filter(self):

        # create cube with just one field and some QC flags
        latitude = DimCoord.from_regular(zeroth=-45.0, step=5.0, count=3, var_name='latitude', circular=False)
        longitude = DimCoord.from_regular(zeroth=-10.0, step=5.0, count=5, var_name='longitude', circular=False)
        axes_coordinates = [ (latitude, 0) , (longitude, 1) ]
        fielddimensions = (len(latitude.points), len(longitude.points))

        # a set of data (all flagged as valid)
        zerodata = numpy.ma.masked_array(
            data=numpy.zeros(fielddimensions, dtype=numpy.float32), 
            mask=numpy.zeros(fielddimensions, dtype=numpy.bool))

        # example of QC flags
        qc = numpy.array([[0, 1, 3, 7, 2],
                          [1, 0, 8, 3, 3],
                          [7, 2, 6, 3, 1]], numpy.int32)

        # build into cube
        cube = Cube(zerodata, var_name='EXAMPLE', dim_coords_and_dims=axes_coordinates)
        cube.add_aux_coord(AuxCoord(qc, var_name='QC'), data_dims=(0,1))

        # copies of cube for tests
        collection_a = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
        collection_b = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
        collection_c = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
    
        # run with different filters

        collection_a.apply_filter(collection_a.get_filter_from_qc_flags(qc_mask=0xFFFF, qc_filter=0x0008))
        numpy.testing.assert_equal(
            numpy.array([ [ True, True, True, True, True ],
                          [ True, True, False, True, True ],
                          [ True, True, True, True, True ] ]),
            collection_a.cube.data.mask)

        collection_b.apply_filter(collection_b.get_filter_from_qc_flags(qc_mask=0x0004, qc_filter=0x0004))
        numpy.testing.assert_equal(
            numpy.array([ [ True, True, True, False, True ],
                          [ True, True, True, True, True ],
                          [ False, True, False, True, True ] ]),
            collection_b.cube.data.mask)

        collection_c.apply_filter(collection_c.get_filter_from_qc_flags(qc_mask=0x0003, qc_filter=0x0002))
        numpy.testing.assert_equal(
            numpy.array([ [ True, True, True, True, False ],
                          [ True, True, True, True, True ],
                          [ True, False, False, True, True ] ]),
            collection_c.cube.data.mask)

    def test_generate_safety_mask(self):

        # create cube with just one field and some QC flags
        latitude = DimCoord.from_regular(zeroth=-45.0, step=5.0, count=3, var_name='right_latitude', circular=False)
        longitude = DimCoord.from_regular(zeroth=-10.0, step=5.0, count=5, var_name='right_longitude', circular=False)
        axes_coordinates = [ (latitude, 0) , (longitude, 1) ]
        fielddimensions = (len(latitude.points), len(longitude.points))

        # a set of data (all flagged as valid)
        zerodata = numpy.ma.masked_array(
            data=numpy.zeros(fielddimensions, dtype=numpy.float32), 
            mask=numpy.zeros(fielddimensions, dtype=numpy.bool))

	# faulty coordinates with non overlapping mask
	wrong_latitude = numpy.array([[0, ERRORVALUE, 3, 7, 2],
                                      [1, 0, ERRORVALUE, ERRORVALUE, 3],
                                      [7, 2, 6, 3, ERRORVALUE]], numpy.float32)

	wrong_longitude = numpy.array([[ERRORVALUE, 1, 3, 7, ERRORVALUE],
                                       [1, 0, 8, 3, 3],
                                       [7, 2, 6, ERRORVALUE, 1]], numpy.float32)

        # build into cube
        cube = Cube(zerodata, var_name='EXAMPLE', dim_coords_and_dims=axes_coordinates)
	cube.add_aux_coord(AuxCoord(wrong_latitude, var_name='latitude'), data_dims=(0,1))
	cube.add_aux_coord(AuxCoord(wrong_longitude, var_name='longitude'), data_dims=(0,1))

        # copies of cube for tests
        collection = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
    
        # input coordinates with non overlapping masks

	collection.generate_safety_mask()
	numpy.testing.assert_equal(
            numpy.array([ [ True, True, False, False, True ],
                          [ False, False, True, True, False ],
                          [ False, False, False, True, True ] ]),
            collection.safety_mask)

	# faulty coordinates with overlapping mask

	wrong_latitude = numpy.array([[ERRORVALUE, 1, 3, 7, ERRORVALUE],
                                      [1, 0, 2, 4, 3],
                                      [7, 2, 6, ERRORVALUE, 1]], numpy.float32)

	wrong_longitude = numpy.array([[ERRORVALUE, 1, 3, 7, ERRORVALUE],
                                       [1, 0, 8, 3, 3],
                                       [7, 2, 6, ERRORVALUE, 1]], numpy.float32)

        # build into cube
        cube = Cube(zerodata, var_name='EXAMPLE', dim_coords_and_dims=axes_coordinates)
	cube.add_aux_coord(AuxCoord(wrong_latitude, var_name='latitude'), data_dims=(0,1))
	cube.add_aux_coord(AuxCoord(wrong_longitude, var_name='longitude'), data_dims=(0,1))

        # copies of cube for tests
        collection = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
    
        # input coordinates with non overlapping masks

	collection.generate_safety_mask()
	self.assertTrue(hasattr(collection,'safety_mask'))
	numpy.testing.assert_equal(
            numpy.array([ [ True, False, False, False, True ],
                          [ False, False, False, False, False ],
                          [ False, False, False, True, False ] ]),
            collection.safety_mask)




    def test_apply_filter_with_safety_mask_enabled(self):
        # create cube with just one field and some QC flags
        latitude = DimCoord.from_regular(zeroth=-45.0, step=5.0, count=3, var_name='right_latitude', circular=False)
        longitude = DimCoord.from_regular(zeroth=-10.0, step=5.0, count=5, var_name='right_longitude', circular=False)
        axes_coordinates = [ (latitude, 0) , (longitude, 1) ]
        fielddimensions = (len(latitude.points), len(longitude.points))

        # a set of data (all flagged as valid)
        zerodata = numpy.ma.masked_array(
            data=numpy.zeros(fielddimensions, dtype=numpy.float32), 
            mask=numpy.zeros(fielddimensions, dtype=numpy.bool))

	# faulty coordinates with non overlapping mask
	wrong_latitude = numpy.array([[0, ERRORVALUE, 3, 7, 2],
                                      [1, 0, ERRORVALUE, ERRORVALUE, 3],
                                      [7, 2, 6, 3, ERRORVALUE]], numpy.float32)

	wrong_longitude = numpy.array([[ERRORVALUE, 1, 3, 7, ERRORVALUE],
                                       [1, 0, 8, 3, 3],
                                       [7, 2, 6, ERRORVALUE, 1]], numpy.float32)

	# example of QC flags
        qc = numpy.array([[0, 1, 3, 7, 2],
                          [1, 0, 8, 3, 3],
                          [7, 2, 6, 3, 1]], numpy.int32)

        # build into cube
        cube = Cube(zerodata, var_name='EXAMPLE', dim_coords_and_dims=axes_coordinates)
        cube.add_aux_coord(AuxCoord(qc, var_name='QC'), data_dims=(0,1))
	cube.add_aux_coord(AuxCoord(wrong_latitude, var_name='latitude'), data_dims=(0,1))
	cube.add_aux_coord(AuxCoord(wrong_longitude, var_name='longitude'), data_dims=(0,1))

        # copies of cube for tests
        collection_a = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
        collection_b = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
        collection_c = SatelliteCollection(TestSatelliteCollection.FieldNamesForTest(), cube.copy())
    
        # run with different filters
	collection_a.generate_safety_mask()
	collection_b.generate_safety_mask()
	collection_c.generate_safety_mask()

        collection_a.apply_filter(collection_a.get_filter_from_qc_flags(qc_mask=0xFFFF, qc_filter=0x0008))
        numpy.testing.assert_equal(
            numpy.array([ [ True, True, True, True, True ],
                          [ True, True, True, True, True ],
                          [ True, True, True, True, True ] ]),
            collection_a.cube.data.mask)
        collection_b.apply_filter(collection_b.get_filter_from_qc_flags(qc_mask=0x0004, qc_filter=0x0004))
        numpy.testing.assert_equal(
            numpy.array([ [ True, True, True, False, True ],
                          [ True, True, True, True, True ],
                          [ False, True, False, True, True ] ]),
            collection_b.cube.data.mask)

        collection_c.apply_filter(collection_c.get_filter_from_qc_flags(qc_mask=0x0003, qc_filter=0x0002))
        numpy.testing.assert_equal(
            numpy.array([ [ True, True, True, True, True ],
                          [ True, True, True, True, True ],
                          [ True, False, False, True, True ] ]),
            collection_c.cube.data.mask)