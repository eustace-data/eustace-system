"""Test the code for processing commmand lines."""
# pylint: disable=duplicate-code,bad-whitespace

import unittest
import tempfile
import numpy
from netCDF4 import Dataset
from .. import command


class TestCommandLine(unittest.TestCase):
    """Test the command module."""

    class PathPatternDate(object):
        """Data class for file name patterns."""

        def __init__(self, path, pattern_lst, pattern_aux, date):
            self.path = path
            self.pattern_lst = pattern_lst
            self.pattern_aux = pattern_aux
            self.date = date

    class AgOptions(object):
        """Data class to simulate a namespace of command line arguments."""

        # Allow short names and many instance attributes, as these are required to
        # simulate command line arguments
        # pylint: disable=invalid-name, too-many-instance-attributes, too-many-locals, too-many-arguments

        def __init__(self, x0, xs, xn, y0, ys, yn, wrapcoords, qc_mask_obs, qc_filter_obs, qc_mask_valid, qc_filter_valid, obsname, maxbinobs, source='', o=None, pattern_lst=None, pattern_aux=None, date=None, path=None, complevel=0):
            self.x0 = x0
            self.xs = xs
            self.xn = xn
            self.y0 = y0
            self.ys = ys
            self.yn = yn
            self.wrapcoords = wrapcoords
            self.qc_mask_obs = qc_mask_obs
            self.qc_filter_obs = qc_filter_obs
            self.qc_mask_valid = qc_mask_valid
            self.qc_filter_valid = qc_filter_valid
            self.obsname = obsname
            self.maxbinobs = maxbinobs
            self.source = source
            self.o = o
            self.pattern_lst = pattern_lst
            self.pattern_aux = pattern_aux
            self.date = date
            self.path = path
            self.complevel = complevel

    # Don't require docstrings for tests if the test names are descriptive enough
    # pylint: disable=missing-docstring

    def test_build_one_filename(self):
        result = command.FilenameList(TestCommandLine.PathPatternDate(
            'bob', 'afile.%Y-%m-%d.%H_%M.nc', 'auxstuff.%Y%m%d.%H%M.nc', '193311031233'))
        self.assertEqual('single', result.time_mode)
        self.assertEqual(1, len(result.filenames))
        self.assertEqual('bob/afile.1933-11-03.12_33.nc',
                         result.filenames[0].filename_lst)
        self.assertEqual('bob/auxstuff.19331103.1233.nc',
                         result.filenames[0].filename_aux)

    def test_build_wholeday_filenames(self):

        result = command.FilenameList(
            TestCommandLine.PathPatternDate('who', 'me.%Y%m%d.%H%M.nc', 'yes.%Y%m%d.%H%M.nc', '20010723'))
        self.assertEqual('day', result.time_mode)
        self.assertEqual(288, len(result.filenames))
        self.assertEqual('who/me.20010723.0000.nc', result.filenames[0].filename_lst)
        self.assertEqual('who/yes.20010723.0000.nc', result.filenames[0].filename_aux)
        self.assertEqual('who/me.20010723.0005.nc', result.filenames[1].filename_lst)
        self.assertEqual('who/yes.20010723.0005.nc', result.filenames[1].filename_aux)
        self.assertEqual('who/me.20010723.0110.nc', result.filenames[14].filename_lst)
        self.assertEqual('who/yes.20010723.0110.nc', result.filenames[14].filename_aux)
        self.assertEqual('who/me.20010723.2350.nc', result.filenames[286].filename_lst)
        self.assertEqual('who/yes.20010723.2350.nc', result.filenames[286].filename_aux)
        self.assertEqual('who/me.20010723.2355.nc', result.filenames[287].filename_lst)
        self.assertEqual('who/yes.20010723.2355.nc', result.filenames[287].filename_aux)

        # causes problems if relying on exception thrown by
        # datetime.strptime(args.date, '%Y%m%d%H%M')
        # to detect single file / whole day mode
        result = command.FilenameList(TestCommandLine.PathPatternDate('a', 'b.%Y%m%d.%H%M.nc', 'c.%Y%m%d.%H%M.nc', '20101102'))
        self.assertEqual('day', result.time_mode)
        self.assertEqual(288, len(result.filenames))
        self.assertEqual('a/b.20101102.0000.nc', result.filenames[0].filename_lst)
        self.assertEqual('a/c.20101102.0000.nc', result.filenames[0].filename_aux)
        self.assertEqual('a/b.20101102.2355.nc', result.filenames[287].filename_lst)
        self.assertEqual('a/c.20101102.2355.nc', result.filenames[287].filename_aux)

    def test_run(self):

        # Lots of local variables required here for tests
        # pylint: disable=too-many-locals, too-many-statements

        X = 250.00 # pylint: disable=invalid-name
        observations = numpy.array(
            [[X,      X,      X,   280.05,    X     ],
             [X,      X,   280.21, 279.43, 280.14   ],
             [279.73, 279.73, 280.27, 279.81, 279.39],
             [280.00, 280.13,    X,   280.08, 279.85],
             [280.07,    X,      X,      X,   279.97]], numpy.float32) # noqa

        qc_cloud = numpy.array(
            [[1, 1, 1, 0, 1],
             [1, 1, 0, 0, 0],
             [0, 0, 0, 0, 0],
             [0, 0, 1, 0, 0],
             [0, 1, 1, 1, 0]], numpy.int32)

        uncertainty_random = numpy.zeros((5, 5), numpy.float32)

        uncertainty_local_atm = numpy.zeros((5, 5), numpy.float32)

        uncertainty_local_sfc = numpy.zeros((5, 5), numpy.float32)

        uncertainty_systematic = numpy.array([0.0], numpy.float32)

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
        r = Dataset(generated_file_lst.name, 'w', 'NETCDF4') # pylint: disable=invalid-name
        r.createDimension('nj', 5)
        r.createDimension('ni', 5)
        lat = r.createVariable('lat', 'f4', ('nj', 'ni'))
        lat[:] = coord_lat
        lon = r.createVariable('lon', 'f4', ('nj', 'ni'))
        lon[:] = coord_lon
        lst = r.createVariable('LST', 'f4', ('nj', 'ni',))
        lst[:] = observations
        qc = r.createVariable('QC', 'i4', ('nj', 'ni',)) # pylint: disable=invalid-name
        qc[:] = qc_cloud
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

        # output file to write (automatically removed when out of scope)
        fileoutput = tempfile.NamedTemporaryFile(
            prefix='satgrid.output.', suffix='.nc')

        # make aggregator with:
        #
        #   latitude grid box boundaries  50.00, 50.25, 50.50, 50.75, 60.00
        #   longitude grid box boundaries -4.00, -3.75, -3.50, -3.25, -3.00
        #
        # so that all of our data falls into box (2, 2)
        #

        options = TestCommandLine.AgOptions(
            y0=50.0, yn=4, ys=0.25,
            x0=-4.0, xn=4, xs=0.25,
            wrapcoords=False,
            qc_mask_obs=64,
            qc_filter_obs=0,
            qc_mask_valid=65,
            qc_filter_valid=0,
            obsname='LST',
            maxbinobs=25,
            o=fileoutput.name,
            pattern_lst=generated_file_lst.name,
            pattern_aux=generated_file_aux.name,
            date='201511251200',
            path='',
            source='made up',
            complevel=4)

        # run it
        command.run(options)

        # re-load output
        result = Dataset(fileoutput.name, 'r', 'NETCDF4')

        # result.variables can have a subscript though pylint doesn't detect this correctly
        # pylint: disable=unsubscriptable-object

        # for v in result.variables:
        #  print v
        #  print result.variables[v][:]

        # check results
        self.assertEqual('made up', result.source)
        self.assertEqual(numpy.int16, result.variables['tsmean'].datatype)
        numpy.testing.assert_equal(
            [[True, True, True, True],
             [True, True, True, True],
             [True, True, False, True],
             [True, True, True, True]],
            result.variables['tsmean'][:].mask)

        numpy.testing.assert_almost_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,       -1],
                 [-1, -1, -1,       -1],
                 [-1, -1,  279.924, -1],
                 [-1, -1, -1,       -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsmean'][:].data,
            decimal=3)

        self.assertEqual(numpy.int16, result.variables['tsmax'].datatype)
        numpy.testing.assert_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  280.27,  -1],
                 [-1, -1, -1,    -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsmax'][:])

        self.assertEqual(numpy.int16, result.variables['tsmin'].datatype)
        numpy.testing.assert_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  279.39,  -1],
                 [-1, -1, -1,    -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsmin'][:])

        self.assertEqual(numpy.int16, result.variables['tsmean_unc_ran'].datatype)
        numpy.testing.assert_almost_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.0,  -1],
                 [-1, -1, -1,    -1]]),
            result.variables['tsmean_unc_ran'][:],
            decimal=3)

        self.assertEqual(numpy.int16, result.variables['tsmean_unc_loc_atm'].datatype)
        numpy.testing.assert_almost_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.0,  -1],
                 [-1, -1, -1,    -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsmean_unc_loc_atm'][:],
            decimal=3)

        self.assertEqual(numpy.int16, result.variables['tsmean_unc_loc_sfc'].datatype)
        numpy.testing.assert_almost_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.0,  -1],
                 [-1, -1, -1,    -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsmean_unc_loc_sfc'][:],
            decimal=3)

        self.assertEqual(numpy.int16, result.variables['tsmean_unc_sys'].datatype)
        numpy.testing.assert_almost_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,    -1],
                 [-1, -1, -1,    -1],
                 [-1, -1,  0.0,  -1],
                 [-1, -1, -1,    -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsmean_unc_sys'][:],
            decimal=6)

        # variance should be computed based on valid obs
        # - value taken from spreadsheet
        # - note precision modified to 3dp here due to int16 storage
        numpy.testing.assert_almost_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,         -1],
                 [-1, -1, -1,         -1],
                 [-1, -1,  0.070211,  -1],
                 [-1, -1, -1,         -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsvariance'][:],
            decimal=3)

        # sampling uncertainty should *not* be zero here
        # - value taken from excel calculation
        # - note precision modified to 3dp here due to int16 storage
        self.assertEqual(numpy.int16, result.variables['tsmean_unc_spl'].datatype)
        numpy.testing.assert_almost_equal(
            numpy.ma.masked_array(
                [[-1, -1, -1,         -1],
                 [-1, -1, -1,         -1],
                 [-1, -1,  0.110406,  -1],
                 [-1, -1, -1,         -1]],
                [[True, True, True, True],
                 [True, True, True, True],
                 [True, True, False, True],
                 [True, True, True, True]]),
            result.variables['tsmean_unc_spl'][:],
            decimal=3)

        # check count (of valid observations)
        self.assertEqual(numpy.int16, result.variables['ts_number_of_observations'].datatype)
        numpy.testing.assert_equal(
            numpy.array(
                [[0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 15, 0],
                 [0, 0, 0, 0]], numpy.int32),
            result.variables['ts_number_of_observations'][:])

        # check count (total number of obs, including ones removed due to
        # cloud)
        self.assertEqual(numpy.int16, result.variables['total_number_of_observations'].datatype)
        numpy.testing.assert_equal(
            numpy.array(
                [[0, 0, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 25, 0],
                 [0, 0, 0, 0]], numpy.int32),
            result.variables['total_number_of_observations'][:])

        # check meta-data

        self.assertEqual('Mean surface temperature (K)', result.variables['tsmean'].long_name)
        self.assertEqual('K', result.variables['tsmean'].units)
        self.assertEqual('surface_temperature', result.variables['tsmean'].standard_name)

        # self.assertEqual('Surface temperature variance (K2)', result.variables['tsvariance'].long_name)
        # self.assertEqual('K2', result.variables['tsvariance'].units)

        self.assertEqual('Minimum surface temperature (K)', result.variables['tsmin'].long_name)
        self.assertEqual('K', result.variables['tsmin'].units)
        self.assertEqual('surface_temperature', result.variables['tsmin'].standard_name)

        self.assertEqual('Maximum surface temperature (K)', result.variables['tsmax'].long_name)
        self.assertEqual('K', result.variables['tsmax'].units)
        self.assertEqual('surface_temperature', result.variables['tsmax'].standard_name)

        self.assertEqual('Random uncertainty in mean surface temperature (K)', result.variables['tsmean_unc_ran'].long_name)
        self.assertEqual('K', result.variables['tsmean_unc_ran'].units)

        self.assertEqual('Locally correlated uncertainty (atm component) in mean surface temperature (K)', result.variables['tsmean_unc_loc_atm'].long_name)
        self.assertEqual('K', result.variables['tsmean_unc_loc_atm'].units)

        self.assertEqual('Locally correlated uncertainty (sfc component) in mean surface temperature (K)', result.variables['tsmean_unc_loc_sfc'].long_name)
        self.assertEqual('K', result.variables['tsmean_unc_loc_sfc'].units)

        self.assertEqual('Systematic uncertainty in mean surface temperature (K)', result.variables['tsmean_unc_sys'].long_name)
        self.assertEqual('K', result.variables['tsmean_unc_sys'].units)

        self.assertEqual('Sampling uncertainty in mean surface temperature (K)', result.variables['tsmean_unc_spl'].long_name)
        self.assertEqual('K', result.variables['tsmean_unc_spl'].units)

        self.assertEqual('Number of observations of surface temperature', result.variables['ts_number_of_observations'].long_name)
        self.assertEqual('1', result.variables['ts_number_of_observations'].units)

        self.assertEqual('Total number of satellite observations in grid box', result.variables['total_number_of_observations'].long_name)
        self.assertEqual('1', result.variables['total_number_of_observations'].units)

        result.close()

    def test_run_compressed(self):
        """Same as test_run but with compression option set."""
        pass
