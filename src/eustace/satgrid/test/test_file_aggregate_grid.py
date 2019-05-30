"""Tests for aggregator module."""
# pylint: disable=duplicate-code,bad-whitespace

__version__ = "$Revision: 472 $"
__author__ = "Joel R. Mitchelson"

import unittest
import tempfile
import math
import numpy
from netCDF4 import Dataset
from ..file_aggregate_grid import SatelliteFileAggregator
from ..satellite_file import SatelliteFile
from ..satellite_file import SatelliteFilename
from ..grid import GridAxis


class TestSatelliteFileAggregator(unittest.TestCase):
    """Test SatelliteFileAggregator class."""

    # Use attribute assignment outside of class as shortcut for creating
    # a simulated command argument namespace, and allow short names as required
    # pylint: disable=attribute-defined-outside-init, invalid-name, too-many-instance-attributes, too-many-locals, too-many-statements

    def test_init(self):
        """Test constructor."""

        # Short name helps write tests concisely
        # pylint: disable=invalid-name

        ag = SatelliteFileAggregator(axis_lat=GridAxis(-180.0, 0.25, 1440), axis_lon=GridAxis(-90.0, 0.1, 720), qc_mask_obs=15,
                                     qc_filter_obs=3, qc_mask_valid=23, qc_filter_valid=2, observation_name='LST', max_bin_obs=4)
        self.assertEqual(-180.0, ag.grid.axis_lat.start)
        self.assertEqual(0.25, ag.grid.axis_lat.increment)
        self.assertEqual(1440, ag.grid.axis_lat.length)
        self.assertEqual(-90.0, ag.grid.axis_lon.start)
        self.assertEqual(0.1, ag.grid.axis_lon.increment)
        self.assertEqual(720, ag.grid.axis_lon.length)
        self.assertEqual(15, ag.qc_mask_obs)
        self.assertEqual(3, ag.qc_filter_obs)
        self.assertEqual(23, ag.qc_mask_valid)
        self.assertEqual(2, ag.qc_filter_valid)
        self.assertEqual('LST', ag.observation_name)
        # self.assertEqual(4, ag.outputgrid.max_bin_obs)

    class TestArguments(object):
        """Shell namespace to hold arguments for test of building from command line."""

        def __init(self):
            pass

    def test_build_from_command_arguments(self):
        """Test the interpretation of command line arguments namespace."""

        args = TestSatelliteFileAggregator.TestArguments()
        args.x0 = -33.33
        args.xn = 79
        args.xs = 0.56
        args.y0 = 5.0
        args.yn = 10000
        args.ys = 100.2
        args.wrapcoords = True
        args.qc_mask_obs = 78911
        args.qc_filter_obs = 87
        args.qc_mask_valid = 2233
        args.qc_filter_valid = 25
        args.obsname = 'bob'
        args.maxbinobs = 7

        ag = SatelliteFileAggregator.from_command_arguments(args)
        self.assertEqual(-33.33, ag.grid.axis_lon.start)
        self.assertEqual(79, ag.grid.axis_lon.length)
        self.assertEqual(0.56, ag.grid.axis_lon.increment)
        self.assertEqual(True, ag.grid.axis_lon.wrap)
        self.assertEqual(5.0, ag.grid.axis_lat.start)
        self.assertEqual(10000, ag.grid.axis_lat.length)
        self.assertEqual(100.2, ag.grid.axis_lat.increment)
        self.assertEqual(True, ag.grid.axis_lat.wrap)
        self.assertEqual(78911, ag.qc_mask_obs)
        self.assertEqual(87, ag.qc_filter_obs)
        self.assertEqual(2233, ag.qc_mask_valid)
        self.assertEqual(25, ag.qc_filter_valid)
        self.assertEqual('bob', ag.observation_name)
        # self.assertEqual(7, ag.outputgrid.max_bin_obs)

    def test_aggregate_file(self):
        """Example aggregating files."""

        # temp file (usually in /tmp folder): will be destroyed when
        # this variable goes out of scope (even if tests fail)
        generated_file_lst = tempfile.NamedTemporaryFile(
            prefix='satgrid.lst', suffix='.nc')
        generated_file_aux = tempfile.NamedTemporaryFile(
            prefix='satgrid.aux', suffix='.nc')
        # print 'Output file: ', generated_file.name

        # build a netCDF-style file
        r = Dataset(generated_file_lst.name, 'w', 'NETCDF4')
        r.createDimension('nj', 2)
        r.createDimension('ni', 5)
        lat = r.createVariable('lat', 'f4', ('nj', 'ni'))
        lat[:] = numpy.array([[-0.1, -0.1, 0.1, 0.1, 1.1],
                              [10.1, 10.2, 10.6, 15.1, 16.1]], numpy.float32)
        lon = r.createVariable('lon', 'f4', ('nj', 'ni'))
        lon[:] = numpy.array([[1.0, 2.0, 3.0, 3.9, 5.0],
                              [1.0, 3.0, 2.0, 4.0, 6.0]], numpy.float32)
        lst = r.createVariable('LST', 'f4', ('nj', 'ni',))
        lst[:] = numpy.array([[18.0, 5.0, 33.0, 44.0, 55.0],
                              [-1.0, 23.0, -89.0, 22.0, 65.0]], numpy.float32)
        qc = r.createVariable('QC', 'i4', ('nj', 'ni',))
        qc[:] = numpy.array([[0, 3, 3, 3, 0],
                             [3, 1, 0, 3, 0]], numpy.int32)
        r.close()

        # and accompanying aux file with uncertainty information
        aux = Dataset(generated_file_aux.name, 'w', 'NETCDF4')
        aux.createDimension('ns', 1)
        aux.createDimension('nj', 2)
        aux.createDimension('ni', 5)
        lst_unc_sys = aux.createVariable('LST_unc_sys', 'f4', ('ns',))
        lst_unc_sys[:] = numpy.float32(0.2)
        lst_unc_ran = aux.createVariable('LST_unc_ran', 'f4', ('nj', 'ni',))
        lst_unc_ran[:] = numpy.array([[0.1, 0.2, 0.3, 0.4, 0.5],
                                      [0.9, 0.8, 0.7, 0.6, 0.5]], numpy.float32)
        lst_unc_loc_atm = aux.createVariable('LST_unc_loc_atm', 'f4', ('nj', 'ni',))
        lst_unc_loc_atm[:] = numpy.array([[0.1, 0.2, 0.3, 0.4, 0.5],
                                          [0.9, 0.8, 0.7, 0.6, 0.5]], numpy.float32)
        lst_unc_loc_sfc = aux.createVariable('LST_unc_loc_sfc', 'f4', ('nj', 'ni',))
        lst_unc_loc_sfc[:] = numpy.array([[0.2, 0.4, 0.6, 0.8, 1.0],
                                          [1.8, 1.6, 1.5, 1.2, 1.0]], numpy.float32)
        aux.close()

        # make aggregator
        ag = SatelliteFileAggregator(
            axis_lat=GridAxis(-5.0, 5.0, 5),
            axis_lon=GridAxis(0.0, 1.0, 6),
            qc_mask_obs=1, qc_filter_obs=1, qc_mask_valid=3, qc_filter_valid=3,
            observation_name='LST',
            max_bin_obs=4)

        # put this file onto it
        ag.aggregate_file(SatelliteFile(SatelliteFilename(generated_file_lst.name, generated_file_aux.name)))

        # the five (lat, lon) coords and values with qc=3 and grid mappings are:
        # ( -0.1, 2.0):  5  ->  ( 0, 2 )
        # (  0.1, 3.0): 33      ( 1, 3 )
        # (  0.1, 3.9): 44      ( 1, 3 )
        # ( 10.1, 1.0): -1      ( 3, 1 )
        # ( 15.1, 4.0): 22      ( 4, 4 )

        # In addition we have one co-ord with qc=1 which is observable but not valid
        # (10.2, 3.0)      ->   ( 3, 3 )

        # count of valid obs
        numpy.testing.assert_equal(
            ag.outputgrid.count,
            numpy.array([[0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 2, 0, 0],
                         [0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0]], numpy.int32))

        # count of total obs
        # note difference in box (3, 3) due to co-ord observable but not valid
        numpy.testing.assert_equal(
            ag.obsgrid.count,
            numpy.array([[0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 2, 0, 0],
                         [0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 1, 0, 0],
                         [0, 0, 0, 0, 1, 0]], numpy.int32))

        numpy.testing.assert_almost_equal(
            ag.outputgrid.values_sum,
            numpy.array([[0.0, 0.0, 5.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 77.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 22.0, 0.0]], numpy.float32))

        # put this file onto it *again* : should double everything!
        ag.aggregate_file(SatelliteFile(SatelliteFilename(
            generated_file_lst.name, generated_file_aux.name)))

        numpy.testing.assert_equal(
            ag.outputgrid.count,
            numpy.array([[0, 0, 2, 0, 0, 0],
                         [0, 0, 0, 4, 0, 0],
                         [0, 0, 0, 0, 0, 0],
                         [0, 2, 0, 0, 0, 0],
                         [0, 0, 0, 0, 2, 0]], numpy.int32))

        numpy.testing.assert_almost_equal(
            ag.outputgrid.values_sum,
            numpy.array([[0.0, 0.0, 10.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 154.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, -2.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 44.0, 0.0]], numpy.float32))

        numpy.testing.assert_equal(
            ag.obsgrid.count,
            numpy.array([[0, 0, 2, 0, 0, 0],
                         [0, 0, 0, 4, 0, 0],
                         [0, 0, 0, 0, 0, 0],
                         [0, 2, 0, 2, 0, 0],
                         [0, 0, 0, 0, 2, 0]], numpy.int32))

    def test_aggregate_uncertainties_all_valid(self):
        """Compare with results from spreadsheet, assuming all obs are valid."""

        observations = numpy.empty((5, 5), numpy.float32)
        observations.fill(288.5)

        uncertainty_random = numpy.array(
            [[0.144, 0.143, 0.147, 0.147, 0.143],
             [0.143, 0.143, 0.143, 0.151, 0.142],
             [0.143, 0.143, 0.141, 0.146, 0.141],
             [0.142, 0.141, 0.146, 0.146, 0.141],
             [0.144, 0.147, 0.145, 0.146, 0.151]], numpy.float32)

        uncertainty_local_atm = numpy.array(
            [[0.032, 0.041, 0.029, 0.055, 0.051],
             [0.000, 0.001, 0.020, 0.114, 0.014],
             [0.000, 0.002, 0.004, 0.048, 0.104],
             [0.002, 0.025, 0.044, 0.031, 0.111],
             [0.028, 0.057, 0.046, 0.090, 0.047]], numpy.float32)

        # Note these are not from excel spreadsheet
        # as it had space for only one uncertainty component
        # They are just invented
        #   sum of squares = 0.29
        #   rms = sqrt(0.29 / 25) = 0.107703296143
        uncertainty_local_sfc = numpy.array(
            [[0.000, 0.000, 0.000, 0.000, 0.000],
             [0.000, 0.200, 0.000, 0.000, 0.000],
             [0.000, 0.000, 0.000, 0.000, 0.000],
             [0.000, 0.000, 0.000, 0.500, 0.000],
             [0.000, 0.000, 0.000, 0.000, 0.000]], numpy.float32)

        uncertainty_systematic = numpy.array([0.3], numpy.float32)

        coord_lat = numpy.array(
            [[50.70, 50.70, 50.70, 50.70, 50.70],
             [50.71, 50.71, 50.71, 50.71, 50.71],
             [50.72, 50.72, 50.72, 50.72, 50.72],
             [50.73, 50.73, 50.73, 50.73, 50.73],
             [50.74, 50.74, 50.74, 50.74, 50.74]], numpy.float32)

        coord_lon = numpy.array(
            [[-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45]], numpy.float32)

        # temp file (usually in /tmp folder): will be destroyed when
        # this variable goes out of scope (even if tests fail)
        generated_file_lst = tempfile.NamedTemporaryFile(
            prefix='satgrid.lst', suffix='.nc')
        generated_file_aux = tempfile.NamedTemporaryFile(
            prefix='satgrid.aux', suffix='.nc')
        # print 'Output file: ', generated_file.name

        # build a netCDF-style file
        r = Dataset(generated_file_lst.name, 'w', 'NETCDF4')
        r.createDimension('nj', 5)
        r.createDimension('ni', 5)
        lat = r.createVariable('lat', 'f4', ('nj', 'ni'))
        lat[:] = coord_lat
        lon = r.createVariable('lon', 'f4', ('nj', 'ni'))
        lon[:] = coord_lon
        lst = r.createVariable('LST', 'f4', ('nj', 'ni',))
        lst[:] = observations
        qc = r.createVariable('QC', 'i4', ('nj', 'ni',))
        qc[:] = numpy.zeros((5, 5), numpy.int32)
        r.close()

        # and accompanying aux file with uncertainty information
        aux = Dataset(generated_file_aux.name, 'w', 'NETCDF4')
        aux.createDimension('ns', 1)
        aux.createDimension('nj', 5)
        aux.createDimension('ni', 5)
        lst_unc_sys = aux.createVariable('LST_unc_sys', 'f4', ('ns',))
        lst_unc_sys[:] = uncertainty_systematic
        lst_unc_ran = aux.createVariable('LST_unc_ran', 'f4', ('nj', 'ni',))
        lst_unc_ran[:] = uncertainty_random
        lst_unc_loc_atm = aux.createVariable('LST_unc_loc_atm', 'f4', ('nj', 'ni',))
        lst_unc_loc_atm[:] = uncertainty_local_atm
        lst_unc_loc_sfc = aux.createVariable('LST_unc_loc_sfc', 'f4', ('nj', 'ni',))
        lst_unc_loc_sfc[:] = uncertainty_local_sfc
        aux.close()

        # make aggregator with:
        #
        #   latitude grid box boundaries  50.00, 50.25, 50.50, 50.75, 60.00
        #   longitude grid box boundaries -4.00, -3.75, -3.50, -3.25, -3.00
        #
        # so that all of our data falls into box (2, 2)
        #
        # assume all obs valid
        #
        ag = SatelliteFileAggregator(
            axis_lat=GridAxis(50.0, 0.25, 4),
            axis_lon=GridAxis(-4.0, 0.25, 4),
            qc_mask_obs=8192, qc_filter_obs=0, qc_mask_valid=8193, qc_filter_valid=0,
            observation_name='LST',
            max_bin_obs=25)

        # put this file onto it
        satfile0 = SatelliteFile(SatelliteFilename(generated_file_lst.name, generated_file_aux.name))
        ag.aggregate_file(satfile0)

        # check counts
        numpy.testing.assert_equal(
            numpy.array(
                [[0,  0,  0,  0],
                 [0,  0,  0,  0],
                 [0,  0, 25,  0],
                 [0,  0,  0,  0]], numpy.int32),
            ag.outputgrid.count)

        # check random noise (taken from excel spreadsheet)
        self.assertAlmostEqual(0.029, (math.sqrt(ag.outputgrid.uncertainty_random_sum_sq[2][2] / 25.0) / math.sqrt(25.0)), places=3)

        # check local (taken from excel spreadsheet)
        self.assertAlmostEqual(0.052, math.sqrt(ag.outputgrid.uncertainty_local_atm_sum_sq[2][2] / 25.0), places=3)

        # check local sfc (from hand computation)
        self.assertAlmostEqual(0.29, ag.outputgrid.uncertainty_local_sfc_sum_sq[2][2], places=3)

        # check systematic error
        numpy.testing.assert_almost_equal(
            numpy.array(
                [[0,  0,  0,     0],
                 [0,  0,  0,     0],
                 [0,  0,  2.25,  0],
                 [0,  0,  0,     0]], numpy.float32),
            ag.outputgrid.uncertainty_systematic_sum_sq,
            decimal=6)

        # check count of all obs is also total
        numpy.testing.assert_array_equal(
            numpy.array(
                [[0,  0,  0,  0],
                 [0,  0,  0,  0],
                 [0,  0, 25,  0],
                 [0,  0,  0,  0]], numpy.int32),
            ag.obsgrid.count)

        # check others zero
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[0][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[0][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[0][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[0][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[1][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[1][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[1][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[1][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[2][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[2][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[2][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[3][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[3][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[3][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_random_sum_sq[3][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[0][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[0][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[0][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[0][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[1][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[1][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[1][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[1][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[2][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[2][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[2][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[3][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[3][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[3][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_atm_sum_sq[3][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[0][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[0][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[0][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[0][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[1][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[1][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[1][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[1][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[2][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[2][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[2][3])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[3][0])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[3][1])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[3][2])
        self.assertAlmostEqual(0.0, ag.outputgrid.uncertainty_local_sfc_sum_sq[3][3])

        # same checks on output fields
        fields0 = ag.compute_fields_first_order()

        # go through file again to compute devsq
        ag.aggregate_file_second_order(fields0, satfile0)

        # do second-order stuff
        fields1 = ag.compute_fields_second_order()

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.029,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            fields0['tsmean_unc_ran'],
            decimal=3)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.052,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            fields0['tsmean_unc_loc_atm'],
            decimal=3)

        #   rms = sqrt(0.29 / 25) = 0.107703296143
        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.107703,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            fields0['tsmean_unc_loc_sfc'],
            decimal=6)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.3,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            fields0['tsmean_unc_sys'],
            decimal=6)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.0,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            fields1['tsmean_unc_spl'],
            decimal=6)

    def test_aggregate_sampling_uncertainty_with_partial_obs_valid(self):
        """Compare with results from spreadsheet, with invalid points."""

        # Obs as in spreadsheet but with addition row,
        # to check behaviour of grid box with one obs
        X = 250.00 # noqa
        observations = numpy.array(
            [[X,      X,      X,      280.05, X     ],            # noqa
             [X,      X,      280.21, 279.43, 280.14],
             [279.73, 279.73, 280.27, 279.81, 279.39],
             [280.00, 280.13, X,      280.08, 279.85],
             [280.07, X,      X,      X,      279.97],
             [X,      X,      290.00, X,      X]], numpy.float32) # noqa

        qc_cloud = numpy.array(
            [[1, 1, 1, 0, 1],
             [1, 1, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [0, 1, 1, 1, 0],
             [1, 1, 0, 1, 1]], numpy.int32)

        uncertainty_random = numpy.zeros((6, 5), numpy.float32)

        uncertainty_local_atm = numpy.zeros((6, 5), numpy.float32)

        uncertainty_local_sfc = numpy.zeros((6, 5), numpy.float32)

        uncertainty_systematic = numpy.array([0.0], numpy.float32)

        coord_lat = numpy.array(
            [[50.70, 50.70, 50.70, 50.70, 50.70],
             [50.71, 50.71, 50.71, 50.71, 50.71],
             [50.72, 50.72, 50.72, 50.72, 50.72],
             [50.73, 50.73, 50.73, 50.73, 50.73],
             [50.74, 50.74, 50.74, 50.74, 50.74],
             [50.01, 50.02, 50.03, 50.04, 50.05]], numpy.float32)

        coord_lon = numpy.array(
            [[-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.49, -3.48, -3.47, -3.46, -3.45],
             [-3.74, -3.73, -3.72, -3.71, -3.70]], numpy.float32)

        # temp file (usually in /tmp folder): will be destroyed when
        # this variable goes out of scope (even if tests fail)
        generated_file_lst = tempfile.NamedTemporaryFile(
            prefix='satgrid.lst', suffix='.nc')
        generated_file_aux = tempfile.NamedTemporaryFile(
            prefix='satgrid.aux', suffix='.nc')
        # print 'Output file: ', generated_file.name

        # build a netCDF-style file
        r = Dataset(generated_file_lst.name, 'w', 'NETCDF4')
        r.createDimension('nj', 5)
        r.createDimension('ni', 6)
        lat = r.createVariable('lat', 'f4', ('nj', 'ni'))
        lat[:] = coord_lat
        lon = r.createVariable('lon', 'f4', ('nj', 'ni'))
        lon[:] = coord_lon
        lst = r.createVariable('LST', 'f4', ('nj', 'ni',))
        lst[:] = observations
        qc = r.createVariable('QC', 'i4', ('nj', 'ni',))
        qc[:] = qc_cloud
        r.close()

        # and accompanying aux file with uncertainty information
        aux = Dataset(generated_file_aux.name, 'w', 'NETCDF4')
        aux.createDimension('ns', 1)
        aux.createDimension('nj', 5)
        aux.createDimension('ni', 6)
        lst_unc_sys = aux.createVariable('LST_unc_sys', 'f4', ('ns',))
        lst_unc_sys[:] = uncertainty_systematic
        lst_unc_ran = aux.createVariable('LST_unc_ran', 'f4', ('nj', 'ni',))
        lst_unc_ran[:] = uncertainty_random
        lst_unc_loc_atm = aux.createVariable('LST_unc_loc_atm', 'f4', ('nj', 'ni',))
        lst_unc_loc_atm[:] = uncertainty_local_atm
        lst_unc_loc_sfc = aux.createVariable('LST_unc_loc_sfc', 'f4', ('nj', 'ni',))
        lst_unc_loc_sfc[:] = uncertainty_local_sfc
        aux.close()

        # make aggregator with:
        #
        #   latitude grid box boundaries  50.00, 50.25, 50.50, 50.75, 60.00
        #   longitude grid box boundaries -4.00, -3.75, -3.50, -3.25, -3.00
        #
        # so that all of our data falls into box (2, 2)
        # except for single obs in (0, 1)
        #
        ag = SatelliteFileAggregator(
            axis_lat=GridAxis(50.0, 0.25, 4),
            axis_lon=GridAxis(-4.0, 0.25, 4),
            qc_mask_obs=512, qc_filter_obs=0, qc_mask_valid=1, qc_filter_valid=0,
            observation_name='LST',
            max_bin_obs=25)

        # put this file onto it
        satfile0 = SatelliteFile(SatelliteFilename(generated_file_lst.name, generated_file_aux.name))
        ag.aggregate_file(satfile0)

        # check counts

        # 15 valid obs pixels in box (2, 2)
        # and one in box (0, 1)
        numpy.testing.assert_equal(
            numpy.array(
                [[0,  1,  0,  0],
                 [0,  0,  0,  0],
                 [0,  0, 15,  0],
                 [0,  0,  0,  0]], numpy.int32),
            ag.outputgrid.count)

        # Total obs is still 25 in (2, 2) and 5 possible in (0, 1)
        numpy.testing.assert_equal(
            numpy.array(
                [[0,  5,  0,  0],
                 [0,  0,  0,  0],
                 [0,  0, 25,  0],
                 [0,  0,  0,  0]], numpy.int32),
            ag.obsgrid.count)

        # systematic, random, locally correlated errors should all be zero
        numpy.testing.assert_almost_equal(
            numpy.zeros((4, 4), numpy.float32),
            ag.outputgrid.uncertainty_random_sum_sq,
            decimal=6)
        numpy.testing.assert_almost_equal(
            numpy.zeros((4, 4), numpy.float32),
            ag.outputgrid.uncertainty_local_atm_sum_sq,
            decimal=6)
        numpy.testing.assert_almost_equal(
            numpy.zeros((4, 4), numpy.float32),
            ag.outputgrid.uncertainty_local_sfc_sum_sq,
            decimal=6)
        numpy.testing.assert_almost_equal(
            numpy.zeros((4, 4), numpy.float32),
            ag.outputgrid.uncertainty_systematic_sum_sq,
            decimal=6)

        # same checks on output fields
        fields0 = ag.compute_fields_first_order()

        # go through file again to compute devsq
        ag.aggregate_file_second_order(fields0, satfile0)

        # do second-order stuff
        fields1 = ag.compute_fields_second_order()

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, 0.0, -1,    -1],
                 [-1,  -1, -1,    -1],
                 [-1,  -1,  0.0,  -1],
                 [-1,  -1, -1,    -1]], numpy.float32),
            fields0['tsmean_unc_ran'],
            decimal=3)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, 0.0, -1,    -1],
                 [-1, -1,  -1,    -1],
                 [-1, -1,   0.0,  -1],
                 [-1, -1,  -1,    -1]], numpy.float32),
            fields0['tsmean_unc_loc_atm'],
            decimal=3)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, 0.0, -1,    -1],
                 [-1, -1,  -1,    -1],
                 [-1, -1,   0.0,  -1],
                 [-1, -1,  -1,    -1]], numpy.float32),
            fields0['tsmean_unc_loc_sfc'],
            decimal=3)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, 0.0, -1,    -1],
                 [-1,  -1, -1,    -1],
                 [-1,  -1,  0.0,  -1],
                 [-1,  -1, -1,    -1]], numpy.float32),
            fields0['tsmean_unc_sys'],
            decimal=6)

        # variance should be computed,
        # but and is invalid for (0, 1) as it has only one valid obs
        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,         -1],
                 [-1, -1, -1,         -1],
                 [-1, -1,  0.070211,  -1],
                 [-1, -1, -1,         -1]], numpy.float32),
            fields1['tsvariance'],
            decimal=6)

        # sampling uncertainty taken from spreadsheet for (2, 2)
        # but invalid for (0, 1) as it has only one valid obs
        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,         -1],
                 [-1, -1, -1,         -1],
                 [-1, -1,  0.110406,  -1],
                 [-1, -1, -1,         -1]], numpy.float32),
            fields1['tsmean_unc_spl'],
            decimal=6)

    def test_run(self):
        """Compare with results from spreadsheet, with cloud cover."""

        # As in spreadsheet but with additional column to check single item
        # behaviour
        X = 250.00 # noqa
        observations = numpy.array(
            [[X,      X,      X,      280.05, X,      X],                 # noqa
             [X,      X,      280.21, 279.43, 280.14, X],
             [279.73, 279.73, 280.27, 279.81, 279.39, 260.0],
             [280.00, 280.13, X,      280.08, 279.85, X],
             [280.07, X,      X,      X,      279.97, X]], numpy.float32) # noqa

        qc_cloud = numpy.array(
            [[1, 1, 1, 0, 1, 1],
             [1, 1, 0, 0, 0, 1],
             [0, 0, 0, 0, 0, 0],
             [0, 0, 1, 0, 0, 1],
             [0, 1, 1, 1, 0, 1]], numpy.int32)

        uncertainty_random = numpy.zeros((5, 6), numpy.float32)

        uncertainty_local_atm = numpy.zeros((5, 6), numpy.float32)

        uncertainty_local_sfc = numpy.zeros((5, 6), numpy.float32)

        uncertainty_systematic = numpy.array([0.0], numpy.float32)

        coord_lat = numpy.array(
            [[50.70, 50.70, 50.70, 50.70, 50.70, 50.26],
             [50.71, 50.71, 50.71, 50.71, 50.71, 50.27],
             [50.72, 50.72, 50.72, 50.72, 50.72, 50.28],
             [50.73, 50.73, 50.73, 50.73, 50.73, 50.29],
             [50.74, 50.74, 50.74, 50.74, 50.74, 50.30]], numpy.float32)

        coord_lon = numpy.array(
            [[-3.49, -3.48, -3.47, -3.46, -3.45, -3.74],
             [-3.49, -3.48, -3.47, -3.46, -3.45, -3.73],
             [-3.49, -3.48, -3.47, -3.46, -3.45, -3.72],
             [-3.49, -3.48, -3.47, -3.46, -3.45, -3.71],
             [-3.49, -3.48, -3.47, -3.46, -3.45, -3.70]], numpy.float32)

        # temp file (usually in /tmp folder): will be destroyed when
        # this variable goes out of scope (even if tests fail)
        generated_file_lst = tempfile.NamedTemporaryFile(
            prefix='satgrid.lst', suffix='.nc')
        generated_file_aux = tempfile.NamedTemporaryFile(
            prefix='satgrid.aux', suffix='.nc')
        # print 'Output file: ', generated_file.name

        # build a netCDF-style file
        r = Dataset(generated_file_lst.name, 'w', 'NETCDF4')
        r.createDimension('nj', 6)
        r.createDimension('ni', 5)
        lat = r.createVariable('lat', 'f4', ('nj', 'ni'))
        lat[:] = coord_lat
        lon = r.createVariable('lon', 'f4', ('nj', 'ni'))
        lon[:] = coord_lon
        lst = r.createVariable('LST', 'f4', ('nj', 'ni',))
        lst[:] = observations
        qc = r.createVariable('QC', 'i4', ('nj', 'ni',))
        qc[:] = qc_cloud
        r.close()

        # and accompanying aux file with uncertainty information
        aux = Dataset(generated_file_aux.name, 'w', 'NETCDF4')
        aux.createDimension('ns', 1)
        aux.createDimension('nj', 6)
        aux.createDimension('ni', 5)
        lst_unc_sys = aux.createVariable('LST_unc_sys', 'f4', ('ns',))
        lst_unc_sys[:] = uncertainty_systematic
        lst_unc_ran = aux.createVariable('LST_unc_ran', 'f4', ('nj', 'ni',))
        lst_unc_ran[:] = uncertainty_random
        lst_unc_loc_atm = aux.createVariable('LST_unc_loc_atm', 'f4', ('nj', 'ni',))
        lst_unc_loc_atm[:] = uncertainty_local_atm
        lst_unc_loc_sfc = aux.createVariable('LST_unc_loc_sfc', 'f4', ('nj', 'ni',))
        lst_unc_loc_sfc[:] = uncertainty_local_sfc
        aux.close()

        # make aggregator with:
        #
        #   latitude grid box boundaries  50.00, 50.25, 50.50, 50.75, 60.00
        #   longitude grid box boundaries -4.00, -3.75, -3.50, -3.25, -3.00
        #
        # so that all of our data falls into box (2, 2)
        # except single item in (1, 1)
        #
        ag = SatelliteFileAggregator(
            axis_lat=GridAxis(50.0, 0.25, 4),
            axis_lon=GridAxis(-4.0, 0.25, 4),
            qc_mask_obs=2, qc_filter_obs=0, qc_mask_valid=1, qc_filter_valid=0,
            observation_name='LST',
            max_bin_obs=25)

        # this file name
        filename0 = SatelliteFilename(
            generated_file_lst.name, generated_file_aux.name)

        # list of files
        satfiles = [filename0]

        # run method
        result = ag.run(satfiles)

        # check results
        self.assertEqual(1, len(result.load))
        self.assertEqual('ok', result.load[0].status)

        # total of 25 items in (2,2) and 5 in (1,1)
        numpy.testing.assert_equal(
            numpy.array(
                [[0, 0, 0, 0],
                 [0, 5, 0, 0],
                 [0, 0, 25, 0],
                 [0, 0, 0, 0]], numpy.int32),
            result.fields['total_number_of_observations'])

        # only 15 in (2,2) are valid, and 1 in (1,1)
        numpy.testing.assert_equal(
            numpy.array(
                [[0, 0, 0, 0],
                 [0, 1, 0, 0],
                 [0, 0, 15, 0],
                 [0, 0, 0, 0]], numpy.int32),
            result.fields['ts_number_of_observations'])

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1,    -1, -1,        -1],
                 [-1, 260.0, -1,        -1],
                 [-1,    -1,  279.924,  -1],
                 [-1,    -1, -1,        -1]], numpy.float32),
            result.fields['tsmean'],
            decimal=3)

        numpy.testing.assert_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, 260.0, -1,    -1],
                 [-1, -1,  280.27,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            result.fields['tsmax'])

        numpy.testing.assert_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, 260.0, -1,    -1],
                 [-1, -1,  279.39,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            result.fields['tsmin'])

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1, -1, -1,    -1],
                 [-1, 0.0, -1,    -1],
                 [-1, -1,  0.0,  -1],
                 [-1, -1, -1,    -1]], numpy.float32),
            result.fields['tsmean_unc_ran'],
            decimal=3)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1,  -1, -1,    -1],
                 [-1, 0.0, -1,    -1],
                 [-1,  -1,  0.0,  -1],
                 [-1,  -1, -1,    -1]], numpy.float32),
            result.fields['tsmean_unc_loc_atm'],
            decimal=3)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1,  -1, -1,    -1],
                 [-1, 0.0, -1,    -1],
                 [-1,  -1,  0.0,  -1],
                 [-1,  -1, -1,    -1]], numpy.float32),
            result.fields['tsmean_unc_loc_sfc'],
            decimal=3)

        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1,  -1, -1,    -1],
                 [-1, 0.0, -1,    -1],
                 [-1,  -1,  0.0,  -1],
                 [-1,  -1, -1,    -1]], numpy.float32),
            result.fields['tsmean_unc_sys'],
            decimal=6)

        # variance should be computed based on valid obs
        # - value taken from spreadsheet
        # note *invalid* for box (1,1) as it has only one valid obs
        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1,  -1, -1,         -1],
                 [-1,  -1, -1,         -1],
                 [-1,  -1,  0.070211,  -1],
                 [-1,  -1, -1,         -1]], numpy.float32),
            result.fields['tsvariance'],
            decimal=6)

        # sampling uncertainty should *not* be zero here
        # - value taken from excel calculation
        # note *invalid* for box (1,1) as it has only one valid obs
        numpy.testing.assert_almost_equal(
            numpy.array(
                [[-1,  -1, -1,         -1],
                 [-1,  -1, -1,         -1],
                 [-1,  -1,  0.110406,  -1],
                 [-1,  -1, -1,         -1]], numpy.float32),
            result.fields['tsmean_unc_spl'],
            decimal=6)
