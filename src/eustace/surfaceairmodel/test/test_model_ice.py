import unittest
import numpy
from netCDF4 import Dataset
import os
from ..model_ice_input_processing import IceSurfaceTemperatureQualityControl
from ..model_ice import ModelIce
from ..model_ice import UNITS
from ..model_ice import ICE_CORRELATED_UNCERTAINTIES
from ..model_ice import ICE_CORRELATED_UNCERTAINTIES_LENGTHS
from ..model_ice import ICE_CORRELATED_UNCERTAINTIES_TIMES
from eustace.analysis.observationsource import ObservationSource
import eustaceconfig

class TestIceSurfaceTemperatureQualityControl(unittest.TestCase):
   
    def test_init(self):

        result = IceSurfaceTemperatureQualityControl(
            daynumber=38891,
            lat=numpy.array([ 2.0, 3.0 ]),
            lon=numpy.array([ 2.0, 3.0, 4.0 ]),
             surface_temperature=numpy.ma.masked_array(
                data=[ [ 3.3, 4.4, 5.5 ], [ 3.2, 4.4, 1.5 ] ], mask=[ [ False, True, False ], [ True, False, False ] ]),
             surface_temperature_min=numpy.ma.masked_array(
                data=[ [ 3.2, 4.3, 5.4 ], [ 3.4, 4.6, 1.9 ] ], mask=[ [ True, True, False ], [ True, False, False ] ]),
             surface_temperature_max=numpy.ma.masked_array(
                data=[ [ 3.9, 4.8, 5.7 ], [ 3.8, 4.8, 2.3 ] ], mask=[ [ False, False, False ], [ True, False, True ] ]),
             surface_temperature_std=numpy.ma.masked_array(
                data=[ [ 0.12, 0.22, 0.012 ], [ 1.44, 1.55, 2.11 ] ], mask=[ [ False, True, False ], [ True, False, False ] ]),
             surface_temperature_3d=numpy.array([ [ [ 1.11, 2.22, 3.33 ],
                                                    [ 4.44, 5.55, 6.66 ] ],
                                                  [ [ 4.32, 5.67, 8.23 ],
                                                    [ 0.11, 0.77, 0.89 ] ],
                                                  [ [ 0.35, 0.78, 0.22 ],
                                                    [ 1.11, 2.22, 3.33 ] ] ]),
             Nobs3d=numpy.array([ [ [ 9, 8, 7 ], 
                                    [ 3, 4, 8 ] ], 
                                  [ [ 1, 4, 5 ],
                                    [ 2, 5, 1 ] ],
                                  [ [ 1, 1, 2 ],
                                    [ 3, 3, 4 ] ] ], dtype=numpy.int64),
             land_mask=numpy.array( [ [ 4, 8, 4 ], [ 1, 2, 3 ] ], dtype=numpy.int64),
             sea_ice_fraction=numpy.array([ [ 0.1, 0.2, 0.3 ], [ 0.3, 0.2, 0.1 ] ]),
             uncorrelated_uncertainty=numpy.ma.masked_array(data=[ [ 28.0, 29.0, 30.0 ], [ 31.0, 32.0, 9.2 ] ]),
             synoptically_correlated_uncertainty=numpy.ma.masked_array(data=[ [ 28.0, 29.0, 30.0 ], [ 31.0, 32.0, 9.3 ] ]),
             large_scale_correlated_uncertainty=numpy.ma.masked_array(data=[ [ 28.0, 29.0, 30.0 ], [ 31.0, 32.0, 9.4 ] ]))

        self.assertEqual(result.hemisphere(), ModelIce.NORTHERN_HEMISPHERE)
        self.assertEqual(38891, result.daynumber)
        numpy.testing.assert_equal(result.surface_temperature.data, [ [ 3.3, 4.4, 5.5 ], [ 3.2, 4.4, 1.5 ] ])
        numpy.testing.assert_equal(result.surface_temperature_min.data, [ [ 3.2, 4.3, 5.4 ], [ 3.4, 4.6, 1.9 ] ])
        numpy.testing.assert_equal(result.surface_temperature_max.data, [ [ 3.9, 4.8, 5.7 ], [ 3.8, 4.8, 2.3 ] ])
        numpy.testing.assert_equal(result.surface_temperature_std.data, [ [ 0.12, 0.22, 0.012 ], [ 1.44, 1.55, 2.11 ] ])
        numpy.testing.assert_equal([ [ [ 1.11, 2.22, 3.33 ],
                                       [ 4.44, 5.55, 6.66 ] ],
                                     [ [ 4.32, 5.67, 8.23 ],
                                       [ 0.11, 0.77, 0.89 ] ],
                                     [ [ 0.35, 0.78, 0.22 ],
                                       [ 1.11, 2.22, 3.33 ] ] ],
                                   result.surface_temperature_3d)
        numpy.testing.assert_equal([ [ [ 9, 8, 7 ], 
                                       [ 3, 4, 8 ] ], 
                                     [ [ 1, 4, 5 ],
                                       [ 2, 5, 1 ] ],
                                     [ [ 1, 1, 2 ],
                                       [ 3, 3, 4 ] ] ],
                                   result.Nobs3d)
        numpy.testing.assert_equal([ [ 4, 8, 4 ], [ 1, 2, 3 ] ], result.land_mask)
        numpy.testing.assert_equal([ [ 0.1, 0.2, 0.3 ], [ 0.3, 0.2, 0.1 ] ], result.sea_ice_fraction)
        numpy.testing.assert_equal([ [ 28.0, 29.0, 30.0 ], [ 31.0, 32.0, 9.2 ] ], result.uncorrelated_uncertainty)
        numpy.testing.assert_equal([ [ 28.0, 29.0, 30.0 ], [ 31.0, 32.0, 9.3 ] ], result.synoptically_correlated_uncertainty)
        numpy.testing.assert_equal([ [ 28.0, 29.0, 30.0 ], [ 31.0, 32.0, 9.4 ] ], result.large_scale_correlated_uncertainty)
                                 
    def test_mask_logic(self):

        icemasks = IceSurfaceTemperatureQualityControl(
            daynumber=40999,
            lat=numpy.array([ -3.0, -2.0 ]),
            lon=numpy.array([ 5.0, 6.0 ]),
            surface_temperature=numpy.ma.masked_array(
                data=[ [ 270.35, 280.16 ], [ 273.1, 280.36 ] ], mask=[ [ False, True ], [ True, False ] ]),
            surface_temperature_min=numpy.ma.masked_array(
                data=[ [ 300.0, 270.0 ], [ 230.0, 265.0   ]  ], mask=[ [ True, True ], [ True, False ] ]),
            surface_temperature_max=numpy.ma.masked_array(
                data=[ [ 290.46, 270.0 ], [ 290.69, 290.36 ] ], mask=[ [ False, False ], [ True, False ] ]),
            surface_temperature_std=numpy.ma.masked_array(
                data=[ [ 0.12, 100.0 ], [ 7.07, 7.0700001 ] ], mask=[ [ False, True ], [ False, False ] ]),
            # surface_temperature_3d=numpy.array(
            #   [ [ [ 270.1, 270.2, 270.3, 270.4, 270.5, 270.6, 270.7, 270.0 ], [ 270.2, 270.3,  1E20,  1E20,  1E20,  1E20, 270.0, 270.1 ] ],
            #     [ [  1E20,  1E20, 270.5, 270.6, 270.7, 271.0,  1E20,  1E20 ], [ 270.4, 270.5, 270.6, 270.7, 270.0, 270.1, 270.2, 270.3 ] ] ] ),
            surface_temperature_3d=numpy.ma.masked_array(data=
                [ [ [  270.1,  270.2 ],
                    [   1E20,  270.4 ] ],
                  [ [ 270.2, 270.3 ],
                    [  1E20, 270.5 ] ],
                  [ [ 270.3, 1E20 ],
                    [ 270.5,  270.6 ] ],
                  [ [ 270.4, 1E20 ],
                    [ 270.6, 270.7 ] ],
                  [ [ 270.5, 1E20 ],
                    [ 270.7, 270.0 ] ],
                  [ [ 270.6, 1E20 ],
                    [ 271.0, 270.1 ] ],
                  [ [  270.7, 270.0 ],
                    [  1E20, 270.2 ] ],
                  [ [ 270.0, 270.1 ],
                    [  1E20, 270.3 ] ] ], 
             mask=
                [ [ [False, False ],
                    [ True, False ] ],
                  [ [False, False ],
                    [ True, False ] ],
                  [ [False,  True ],
                    [False, False ] ],
                  [ [False,  True ],
                    [False, False ] ],
                  [ [False,  True ],
                    [False, False ] ],
                  [ [False,  True ],
                    [False, False ] ],
                  [ [False, False ],
                    [ True, False ] ],
                  [ [False, False ],
                    [ True, False ] ] ]),                                                                         
             # Nobs3d=numpy.array([ [ [ 3, 4, 3, 2, 1, 7, 8, 2 ], [ 3, 4, 0, 0, 0, 0, 8, 2 ] ],
             #                      [ [ 0, 0, 6, 1, 7, 2, 0, 0 ], [ 1, 7, 5 ,2, 6, 2, 1, 1 ] ] ], dtype=numpy.int64),
             Nobs3d = numpy.ma.masked_array( data=[
                    [ [ 3, 3 ],
                      [ 0, 0 ] ],
                    [ [ 4, 4 ],
                      [ 0, 7 ] ],
                    [ [ 3, 0 ],
                      [ 6, 5 ] ],
                    [ [ 2, 0 ],
                      [ 1, 2 ] ],
                    [ [ 1, 0 ],
                      [ 7, 6 ] ],
                    [ [ 7, 0 ],
                      [ 2, 2 ] ],
                    [ [ 8, 8 ],
                      [ 0, 1 ] ],
                    [ [ 2, 2 ],
                      [ 0, 1 ] ] ], mask=numpy.zeros((8,2,2), dtype=numpy.bool)),
             land_mask=numpy.array( [ [ 3, 4 ], [ 1, 0 ] ], dtype=numpy.int64),
             sea_ice_fraction=numpy.ma.masked_array(data=[ [ 0.32, 0.29 ], [ 1.0, 0.24 ] ], mask=[ [ False, False ], [ False, False] ]),
             uncorrelated_uncertainty=numpy.zeros((2,2)),
             synoptically_correlated_uncertainty=numpy.zeros((2,2)),
             large_scale_correlated_uncertainty=numpy.zeros((2,2)))

        self.assertEqual(ModelIce.SOUTHERN_HEMISPHERE, icemasks.hemisphere())

        numpy.testing.assert_equal([ [ True, False ], [ False, True ] ], icemasks.land_ice_mask())
        numpy.testing.assert_equal([ [ False, True ], [ True, True ] ], icemasks.sea_ice_mask())
        numpy.testing.assert_equal([ [ False, False ], [ False, True] ], icemasks.ice_mask())
        numpy.testing.assert_equal([ [ False, True ], [ True, True ] ], icemasks.marginal_ice_zone_mask())
        numpy.testing.assert_equal([ [ False, True ], [ False, False] ], icemasks.observations_day_mask())
        numpy.testing.assert_equal([ [ False, False ], [ True, False] ], icemasks.observations_night_mask())
        numpy.testing.assert_equal([ [ False, True ], [ True, False] ], icemasks.observations_day_or_night_mask())

        # Only masked where there is no mask on input (hence last entry True but first still False even though out of range)
        numpy.testing.assert_equal([ [ False, False ], [ False, True ] ], icemasks.temperature_standard_deviation_mask())

        # Only masked where there is no mask on input (hence last entries True but first still False even though out of range)
        numpy.testing.assert_equal([ [ False, False ], [ False, True ] ], icemasks.temperature_value_mask())

        # Mean of all bins is:
        # [ sum(270.1, 270.2, 270.3, 270.4, 270.5, 270.6, 270.7, 270.0)/8, sum(270.2, 270.3, 270.0, 270.1)/4                             ]
        # [ sum(270.5, 270.6, 270.7, 271.0)/4                            , sum(270.4, 270.5, 270.6, 270.7, 270.0, 270.1, 270.2, 270.3)/8 ]
        # =
        # [ 270.35, 270.15 ]
        # [ 270.70, 270.35 ]
        
        # And must be no more than 10 above this (note entry (1,1) is False because it isn't present in surface_temperature)
        numpy.testing.assert_equal([ [ False, False ], [ False, True ] ], icemasks.combined_bin_comparison_mask())
        
        # Mean of day bins is:
        # [ sum(270.3, 270.4, 270.5, 270.6)/4, 0 / 0                             ]
        # [ sum(270.5, 270.6, 270.7, 271.0)/4, sum(270.6, 270.7, 270.0, 270.1)/4 ]
        # =
        # [ 270.45,   --   ]
        # [ 270.70, 270.35 ]
        #

        # And must be no more than 20 above this (also masked when no data)
        numpy.testing.assert_equal([ [ True, False ], [ False, True ] ], icemasks.day_bin_comparison_mask())

        # Mean of night bins is:
        # [ sum( 270.1, 270.2, 270.7, 270.0)/4, sum( 270.2, 270.3, 270.0, 270.1)/4 ]
        # [                              0 / 0, sum( 270.4, 270.5, 270.2, 270.3)/4 ]
        # =
        # [ 270.25, 270.30 ]
        # [   --  , 270.35 ]
        #

        # Original logic asks if abs(daily minimum - min_of_night_bins) < -20,
        # which is never True
        numpy.testing.assert_equal([ [ False, False ], [ False, False ] ], icemasks.night_bin_comparison_mask())


    def test_local_mean(self):

        # 2D data 6 x 6, generated using:
        # data = 0.1 * numpy.int64( 100 * numpy.random.rand(6,6) )
        # mask = numpy.random.randn(6,6) > 0
        testdata = numpy.ma.masked_array(
            data=[[ 2.5,  4.5, 1E20,  7.9,  7.7, 1E20],
                  [ 3.7,  3.2, 1E20,  2.3,  8.8,  6. ],
                  [ 8.2,  9.8, 1E20, 1E20,  0.5,  5.2],
                  [ 5. , 1E20,  9. ,  1.2, 1E20,  0. ],
                  [1E20,  5.6,  0.1, 1E20, 1E20,  7.9],
                  [1E20,  3.6, 1E20,  1.4, 1E20, 1E20]],
            mask=[[False, False,  True, False, False,  True],
                  [False, False,  True, False, False, False],
                  [False, False,  True,  True, False, False],
                  [False,  True, False, False,  True, False],
                  [ True, False, False,  True,  True, False],
                  [ True, False,  True, False,  True,  True]])

        # Compute
        result = IceSurfaceTemperatureQualityControl.local_mean(2,14, testdata)

        # Should have borders which are all set True
        # And there are four inner computations, e.g. first one is:
        #
        # sum(2.5, 4.5, 7.9, 7.7,
        #     3.7, 3.2, 2.3, 8.8,
        #     8.2, 9.8, 0.5,
        #     5.0, 9.0, 1.2,
        #     5.6, 0.1) / 16 = 5.0
        #
        # Giving:
        #
        # [ 5.0,          4.98125       ]
        # [ 4.10769230769 4.52857142857 ]
        #
        # But the counts of values are:
        # [ 16, 16 ]
        # [ 13, 14 ]
        #
        # So bottom left should be masked according to min obs required.
        #
        #
        numpy.testing.assert_equal(result.mask,
                                   [ [  False, False ],
                                     [  True, False ] ])
        numpy.testing.assert_almost_equal(result.data[0,0], 5.0)
        numpy.testing.assert_almost_equal(result.data[0,1], 4.98125)
        numpy.testing.assert_almost_equal(result.data[1,1], 4.52857142857)

    
    def test_cloud_mask(self):

        # Some constant data with:
        # -- one obvious outlier at (3, 3) which is not flagged
        # -- a bad value at (0, 0) which is flagged as such
        # -- a bad value at (2, 3) which is flagged as such
        # -- a bad value at (3, 2) which is is not flagged but will be via the extra mask below
        testdata = numpy.ma.masked_array(
            data=[ [   0.0, 270.0, 270.0, 270.0, 270.0, 270.0 ],
                   [ 270.0, 270.0, 270.0, 270.0, 270.0, 270.0 ],
                   [ 270.0, 270.0, 270.0, 150.0, 270.0, 270.0 ],
                   [ 270.0, 270.0, 150.0, 150.0, 270.0, 270.0 ],
                   [ 270.0, 270.0, 270.0, 270.0, 270.0, 270.0 ],
                   [ 270.0, 270.0, 270.0, 270.0, 270.0, 270.0 ] ],
            mask=[ [  True, False, False, False, False, False ],
                   [ False, False, False, False, False, False ],
                   [ False, False, False,  True, False, False ],
                   [ False, False, False, False, False, False ],
                   [ False, False, False, False, False, False ],
                   [ False, False, False, False, False, False ] ])
                   
        # The additional mask
        extra_mask = numpy.zeros((6,6), dtype=numpy.bool)
        extra_mask[3,2] = True
        
        # Run the outlier detection
        result = IceSurfaceTemperatureQualityControl.cloud_mask(testdata, extra_mask)
        
        # Should have:
        # -- ignored the first cell
        # -- used the mask at (2, 3) so that this is not flagged
        # -- used the extra mask at (3, 2) so that this is not flagged
        # -- flagged (3,3)
        numpy.testing.assert_equal(
            result,
            [ [ False, False, False, False, False, False ],
              [ False, False, False, False, False, False ],
              [ False, False, False, False, False, False ],
              [ False, False, False,  True, False, False ],
              [ False, False, False, False, False, False ],
              [ False, False, False, False, False, False ] ])

    def test_regression_mask_result(self):
        """Compare with DMI's original output."""

        inputfile = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/AASTI_v1_L3_solartime/2007/01/2007-01-01_sh.nc')
        comparisonfile = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_sh_20070101.nc')

        # Definition of known field names in NetCDF (also match the constructor method params)
        inputfieldnames = [ 'lat',
                            'lon',
                            'surface_temperature', 
                            'surface_temperature_min',
                            'surface_temperature_max',
                            'surface_temperature_std',
                            'surface_temperature_3d',
                            'Nobs3d',
                            'land_mask',
                            'sea_ice_fraction',
                            'uncorrelated_uncertainty',
                            'synoptically_correlated_uncertainty',
                            'large_scale_correlated_uncertainty' ]

        # Load those fields
        inputdata = Dataset(inputfile, 'r')
        inputfields = { name: inputdata.variables[name][:] for name in inputfieldnames }
        inputfields = { name: field[0,:,:] if (field.ndim == 3 and not '3d' in name) else field for name, field in inputfields.iteritems() }

        # Construct
        processor = IceSurfaceTemperatureQualityControl(daynumber=0, **inputfields)

        # Load the reference dataset
        comparisondata = Dataset(comparisonfile, 'r')

        # Check mask info
        numpy.testing.assert_equal(processor.land_ice_mask(), comparisondata.variables['landicemask'][0,:,:] == 0)
        numpy.testing.assert_equal(processor.sea_ice_mask(), comparisondata.variables['seaicemask'][0,:,:] == 0)
        numpy.testing.assert_equal(processor.marginal_ice_zone_mask(), comparisondata.variables['MIZmask'][0,:,:] == 0)
        numpy.testing.assert_equal(processor.quality_controlled_surface_temperature().mask, comparisondata.variables['tas'][0,:,:].mask)
        numpy.testing.assert_equal(processor.quality_controlled_surface_temperature_max().mask, comparisondata.variables['tasmax'][0,:,:].mask)
        numpy.testing.assert_equal(processor.quality_controlled_surface_temperature_min().mask, comparisondata.variables['tasmin'][0,:,:].mask)

        #-- debugging tmin
        # badfalse = numpy.logical_and(numpy.logical_not(processor.quality_controlled_surface_temperature_min().mask), comparisondata.variables['tasmin'][0,:,:].mask)
        # print numpy.nonzero(badfalse)
        # badtrue = numpy.logical_and(processor.quality_controlled_surface_temperature_min().mask, numpy.logical_not(comparisondata.variables['tasmin'][0,:,:].mask))
        # print numpy.nonzero(badtrue)

        #-- debugging tmax
        # badfalse = numpy.logical_and(numpy.logical_not(processor.quality_controlled_surface_temperature_max().mask), comparisondata.variables['tasmax'][0,:,:].mask)
        # print 'badfalse: ', numpy.nonzero(badfalse)
        # badtrue = numpy.logical_and(processor.quality_controlled_surface_temperature_max().mask, numpy.logical_not(comparisondata.variables['tasmax'][0,:,:].mask))
        # print 'badtrue: ', numpy.nonzero(badtrue)

    def test_regression_mask_result2(self):
        """Compare with DMI's original output (again)."""

        inputfile = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/AASTI_v1_L3_solartime/2007/02/2007-02-23_sh.nc')
        comparisonfile = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_sh_20070223.nc')

        # Definition of known field names in NetCDF (also match the constructor method params)
        inputfieldnames = [ 'lat',
                            'lon',
                            'surface_temperature', 
                            'surface_temperature_min',
                            'surface_temperature_max',
                            'surface_temperature_std',
                            'surface_temperature_3d',
                            'Nobs3d',
                            'land_mask',
                            'sea_ice_fraction',
                            'uncorrelated_uncertainty',
                            'synoptically_correlated_uncertainty',
                            'large_scale_correlated_uncertainty' ]

        # Load those fields
        inputdata = Dataset(inputfile, 'r')
        inputfields = { name: inputdata.variables[name][:] for name in inputfieldnames }
        inputfields = { name: field[0,:,:] if (field.ndim == 3 and not '3d' in name) else field for name, field in inputfields.iteritems() }

        # Construct
        processor = IceSurfaceTemperatureQualityControl(daynumber=0, **inputfields)

        # Load the reference dataset
        comparisondata = Dataset(comparisonfile, 'r')

        # Check mask info
        numpy.testing.assert_equal(processor.land_ice_mask(), comparisondata.variables['landicemask'][0,:,:] == 0)
        numpy.testing.assert_equal(processor.sea_ice_mask(), comparisondata.variables['seaicemask'][0,:,:] == 0)
        numpy.testing.assert_equal(processor.marginal_ice_zone_mask(), comparisondata.variables['MIZmask'][0,:,:] == 0)
        numpy.testing.assert_equal(processor.quality_controlled_surface_temperature().mask, comparisondata.variables['tas'][0,:,:].mask)
        numpy.testing.assert_equal(processor.quality_controlled_surface_temperature_max().mask, comparisondata.variables['tasmax'][0,:,:].mask)
        numpy.testing.assert_equal(processor.quality_controlled_surface_temperature_min().mask, comparisondata.variables['tasmin'][0,:,:].mask)

    def test_matlab_exact(self):
        """A close representation of original DMI Matlab code - useful for diagnosing any discrepancies."""

        inputfile = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/AASTI_v1_L3_solartime/2007/01/2007-01-01_sh.nc')
        inputdata = Dataset(inputfile, 'r')

        st = inputdata.variables['surface_temperature'][0,:,:]
        stmin = inputdata.variables['surface_temperature_min'][0,:,:]
        stmax = inputdata.variables['surface_temperature_max'][0,:,:]
        stSTD = inputdata.variables['surface_temperature_std'][0,:,:]
        st_bin = inputdata.variables['surface_temperature_3d'][:,:,:]
        Nobs_bin = inputdata.variables['Nobs3d'][:,:,:]
        iceconc =  inputdata.variables['sea_ice_fraction'][0,:,:]
        LandIce_mask = (inputdata.variables['land_mask'][0,:,:] != 1) & (inputdata.variables['land_mask'][0,:,:] != 4)
        SeaIce_mask = (iceconc.data < 0.3) | iceconc.mask
            
        mask_st = numpy.zeros(st.shape, dtype=numpy.bool)
        mask_stmin = numpy.zeros(st.shape, dtype=numpy.bool)
        mask_stmax = numpy.zeros(st.shape, dtype=numpy.bool)

        mask_st |= st.mask
        mask_st |= (LandIce_mask & SeaIce_mask)

        mask_stmin |= stmin.mask
        mask_stmin |= (LandIce_mask & SeaIce_mask)

        mask_stmax |= stmax.mask
        mask_stmax |= (LandIce_mask & SeaIce_mask)

        bd = (stSTD.data > 7.07) & ~stSTD.mask
        mask_st |= bd
        mask_stmin |= bd
        mask_stmax |= bd

        bd = (st.data > 273.15 + 5) & ~st.mask
        mask_st |= bd
        mask_stmin |= bd
        mask_stmax |= bd
            
        nightbins = [0, 1, 6, 7]
        daybins   = [2, 3, 4, 5]
        Nobs_night = numpy.sum(Nobs_bin[nightbins,:,:], axis=0)
        Nobs_day   = numpy.sum(Nobs_bin[daybins,:,:], axis=0)

        mask_st |= (Nobs_night.data == 0) & ~Nobs_night.mask
        mask_st |= (Nobs_day == 0) & ~Nobs_day.mask
        mask_stmin |= (Nobs_night == 0)
        mask_stmax |= (Nobs_day == 0)
                        
        bin_mean_all = numpy.ma.mean(st_bin, axis=0)
        bin_mean_day = numpy.ma.mean(st_bin[daybins,:,:], axis=0)
        bin_mean_night = numpy.ma.mean(st_bin[nightbins,:,:], axis=0)
            
        mask_st |= (numpy.abs(st - bin_mean_all) > 10)
        mask_stmax |= (numpy.abs(stmax-bin_mean_day) > 20)

        testmask = (mask_stmax | SeaIce_mask)[96-2:96+3,1181-2:1181+3]
        testmask[2,2] = True
        square = numpy.ma.masked_array(data=stmax.data[96-2:96+3,1181-2:1181+3], mask=testmask)        
        diff = square.mean()- stmax.data[96,1181]
        self.assertTrue(diff < 10.0)


class TestModelIce(unittest.TestCase):

    def test_load_model_parameters_hemisphere(self):

        # File for test data
        filename = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ice/m_values_sh_TskinSeason_MS9.mat')        

        # Load it
        components = ModelIce.load_hemisphere_model_parameters(filename)
        
        # Comparison with numbers loaded into Octave
        # octave --no-gui
        # load('/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ice/m_values_sh_TskinSeason_MS9.mat')
        # format('long')
        # M{1,:}
        # M{2,:}
        numpy.testing.assert_almost_equal(components['LandIce'][   'tas'], [ 5.703469823462727 , 1.042718691887182 ,-0.422019294681218 ,-0.217379801228344  ])
        numpy.testing.assert_almost_equal(components['LandIce']['tasmin'], [ 4.6226022258302617, 1.0117634736301606, 0.1978928186064804,-0.0311071036770317 ])
        numpy.testing.assert_almost_equal(components['LandIce']['tasmax'], [ 5.081021982055520 , 1.012180518037014 ,-0.912474714994831 ,-0.301923308096037  ])
        numpy.testing.assert_almost_equal(components[ 'SeaIce'][   'tas'], [ 1.413531017232533 , 0.866334674600820 , 0.961263405623792 , 0.758524683645621  ])
        numpy.testing.assert_almost_equal(components[ 'SeaIce']['tasmin'], [ 1.923594466491028 , 0.777119766100686 , 2.375391786874733 , 0.741157664853369  ])
        numpy.testing.assert_almost_equal(components[ 'SeaIce']['tasmax'], [ 0.473114501062541 , 0.880700952893396 , 1.330191023961684 , 0.406837835681902  ])
       
    def test_load_model_parameters(self):

        # File for test data
        modelfilepattern = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ice/m_values_{hemisphere}_TskinSeason_MS9.mat')
        filenames = { hemisphere: modelfilepattern.format(hemisphere=hemisphere) for hemisphere in [ ModelIce.SOUTHERN_HEMISPHERE, ModelIce.NORTHERN_HEMISPHERE ] }

        # Load it
        model = ModelIce()
        model.load_model_parameters(filenames)
        
        # Comparison with numbers loaded into Octave
        # octave --no-gui
        # load('/gws/nopw/j04/eustace/data/internal/surfaceair_model_parameters/ice/m_values_sh_TskinSeason_MS9.mat')
        # format('long')
        # M{1,:}
        # M{2,:}
        numpy.testing.assert_almost_equal(model.components['sh']['LandIce'][   'tas'], [ 5.703469823462727 , 1.042718691887182 ,-0.422019294681218 ,-0.217379801228344  ])
        numpy.testing.assert_almost_equal(model.components['sh']['LandIce']['tasmin'], [ 4.6226022258302617, 1.0117634736301606, 0.1978928186064804,-0.0311071036770317 ])
        numpy.testing.assert_almost_equal(model.components['sh']['LandIce']['tasmax'], [ 5.081021982055520 , 1.012180518037014 ,-0.912474714994831 ,-0.301923308096037  ])
        numpy.testing.assert_almost_equal(model.components['sh'][ 'SeaIce'][   'tas'], [ 1.413531017232533 , 0.866334674600820 , 0.961263405623792 , 0.758524683645621  ])
        numpy.testing.assert_almost_equal(model.components['sh'][ 'SeaIce']['tasmin'], [ 1.923594466491028 , 0.777119766100686 , 2.375391786874733 , 0.741157664853369  ])
        numpy.testing.assert_almost_equal(model.components['sh'][ 'SeaIce']['tasmax'], [ 0.473114501062541 , 0.880700952893396 , 1.330191023961684 , 0.406837835681902  ])
        
        #Testing correct values of correlation ranges        
        indexes=[2, 7, 12] 
        for i,j in enumerate(indexes):	
	  self.assertEqual(model.ANCILLARY_VARIABLES[j].name,ICE_CORRELATED_UNCERTAINTIES[i])
	  self.assertEqual(model.ANCILLARY_VARIABLES[j].length_scale,ICE_CORRELATED_UNCERTAINTIES_LENGTHS[i])
	  self.assertEqual(model.ANCILLARY_VARIABLES[j].time_scale,ICE_CORRELATED_UNCERTAINTIES_TIMES[i])
	  self.assertEqual(model.ANCILLARY_VARIABLES[j].length_scale_units,UNITS[0])
	  self.assertEqual(model.ANCILLARY_VARIABLES[j].time_scale_units,UNITS[1])

    def test_process_hemisphere(self):

        # Files used for test
        modelfiles = { hemisphere: os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ice/m_values_{hemisphere}_TskinSeason_MS9.mat'.format(hemisphere=hemisphere)) for hemisphere in [ ModelIce.SOUTHERN_HEMISPHERE, ModelIce.NORTHERN_HEMISPHERE ] }
        inputfile = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/AASTI_v1_L3_solartime/2007/01/2007-01-01_sh.nc')
        comparisonfile = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_sh_20070101.nc')

        # The daynumber of 2007-01-01 compared with EUSTACE epoch 1850-01-01
        daynumber = 57343

        # Definition of known field names in NetCDF (also match the constructor method params)
        inputfieldnames = [ 'lat',
                            'lon',
                            'surface_temperature', 
                            'surface_temperature_min',
                            'surface_temperature_max',
                            'surface_temperature_std',
                            'surface_temperature_3d',
                            'Nobs3d',
                            'land_mask',
                            'sea_ice_fraction',
                            'uncorrelated_uncertainty',
                            'synoptically_correlated_uncertainty',
                            'large_scale_correlated_uncertainty' ]

        # Load those fields
        inputdata = Dataset(inputfile, 'r')
        inputfields = { name: inputdata.variables[name][:] for name in inputfieldnames }
        inputfields = { name: field[0,:,:] if (field.ndim == 3 and not '3d' in name) else field for name, field in inputfields.iteritems() }

        # Construct
        inputqc = IceSurfaceTemperatureQualityControl(daynumber, **inputfields)        

        # Load the reference dataset
        comparisondata = Dataset(comparisonfile, 'r')

        # Apply model
        model = ModelIce(modelfiles)
        results = model.process_hemisphere(inputqc)

        # Check validity fields
        numpy.testing.assert_equal(results['tas'].mask, comparisondata.variables['tas'][0,:,:].mask)
        numpy.testing.assert_equal(results['tasmin'].mask, comparisondata.variables['tasmin'][0,:,:].mask)
        numpy.testing.assert_equal(results['tasmax'].mask, comparisondata.variables['tasmax'][0,:,:].mask)

        # Do comparisons only where valid data
        tasvalid = numpy.nonzero(~results['tas'].mask)
        tasminvalid = numpy.nonzero(~results['tasmin'].mask)
        tasmaxvalid = numpy.nonzero(~results['tasmax'].mask)

        # Compare observation fields
        numpy.testing.assert_almost_equal(results['tas'][tasvalid], comparisondata.variables['tas'][0,:,:].data[tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin'][tasminvalid], comparisondata.variables['tasmin'][0,:,:].data[tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax'][tasmaxvalid], comparisondata.variables['tasmax'][0,:,:].data[tasmaxvalid], decimal=4)

        # Uncertainty fields should match masks
        numpy.testing.assert_equal(results['tas_unc_rand'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_rand'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_rand'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_sys'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_sys'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_sys'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_cloud'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_cloud'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_cloud'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_corr_local'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_corr_local'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_corr_local'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_no_cloud'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_no_cloud'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_no_cloud'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tasuncertainty'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasminuncertainty'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmaxuncertainty'].mask, results['tasmax'].mask)

        # Uncertainty fields should match where observations valid
        numpy.testing.assert_almost_equal(results['tas_unc_rand'].data[tasvalid], comparisondata.variables['RU'][0,:,:].data[tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_rand'].data[tasminvalid], comparisondata.variables['RUmin'][0,:,:].data[tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_rand'].data[tasmaxvalid], comparisondata.variables['RUmax'][0,:,:].data[tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_corr_local'].data[tasvalid], comparisondata.variables['SSU'][0,:,:].data[tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_corr_local'].data[tasminvalid], comparisondata.variables['SSUmin'][0,:,:].data[tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_corr_local'].data[tasmaxvalid], comparisondata.variables['SSUmax'][0,:,:].data[tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_sys'].data[tasvalid], comparisondata.variables['LSU'][0,:,:].data[tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_sys'].data[tasminvalid], comparisondata.variables['LSUmin'][0,:,:].data[tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_sys'].data[tasmaxvalid], comparisondata.variables['LSUmax'][0,:,:].data[tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_cloud'].data[tasvalid], comparisondata.variables['CU'][0,:,:].data[tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_cloud'].data[tasminvalid], comparisondata.variables['CUmin'][0,:,:].data[tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_cloud'].data[tasmaxvalid], comparisondata.variables['CUmax'][0,:,:].data[tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_no_cloud'].data[tasvalid], comparisondata.variables['TotU'][0,:,:].data[tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_no_cloud'].data[tasminvalid], comparisondata.variables['TotUmin'][0,:,:].data[tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_no_cloud'].data[tasmaxvalid], comparisondata.variables['TotUmax'][0,:,:].data[tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasuncertainty'].data[tasvalid], comparisondata.variables['TotUC'][0,:,:].data[tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasminuncertainty'].data[tasminvalid], comparisondata.variables['TotUCmin'][0,:,:].data[tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmaxuncertainty'].data[tasmaxvalid], comparisondata.variables['TotUCmax'][0,:,:].data[tasmaxvalid], decimal=4)

    def test_process_global(self):

        # Files used for test
        modelfilepattern = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/surfaceair_model_parameters/ice/m_values_{hemisphere}_TskinSeason_MS9.mat')
        inputfilepattern = os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/incoming/AASTI_v1_L3_solartime/2007/01/2007-01-01_{hemisphere}.nc')
        modelfiles = { hemisphere: modelfilepattern.format(hemisphere=hemisphere) for hemisphere in [ ModelIce.SOUTHERN_HEMISPHERE, ModelIce.NORTHERN_HEMISPHERE ] }

        # The daynumber of 2007-01-01 compared with EUSTACE epoch 1850-01-01
        daynumber = 57343

        # Definition of known field names in NetCDF (also match the constructor method params)
        inputfieldnames = [ 'lat',
                            'lon',
                            'surface_temperature', 
                            'surface_temperature_min',
                            'surface_temperature_max',
                            'surface_temperature_std',
                            'surface_temperature_3d',
                            'Nobs3d',
                            'land_mask',
                            'sea_ice_fraction',
                            'uncorrelated_uncertainty',
                            'synoptically_correlated_uncertainty',
                            'large_scale_correlated_uncertainty' ]

        # Input dictionary
        inputs = { }

        # Load those fields
        for hemisphere in [ ModelIce.SOUTHERN_HEMISPHERE, ModelIce.NORTHERN_HEMISPHERE ]:

            # filename for this input
            inputfilename = inputfilepattern.format(hemisphere=hemisphere)

            # Load
            inputdata = Dataset(inputfilename, 'r')

            # Get fields
            inputfields = { name: inputdata.variables[name][:] for name in inputfieldnames }

            # Reduce to 2D those with time dimension
            inputfields = { name: field[0,:,:] if (field.ndim == 3 and not '3d' in name) else field for name, field in inputfields.iteritems() }            

            # Construct
            inputs[hemisphere] = IceSurfaceTemperatureQualityControl(daynumber, **inputfields)

        # Apply model
        results = ModelIce(modelfiles).process_global(inputs[ModelIce.SOUTHERN_HEMISPHERE], inputs[ModelIce.NORTHERN_HEMISPHERE])

        # Data to compare
        ref_sh = Dataset(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_sh_20070101.nc'), 'r')
        ref_nh = Dataset(os.path.join(eustaceconfig.WORKSPACE_PATH, 'data/internal/iat_from_ist/EUSTACE_satstace_v1/2007/eustace_satstace_ice_nh_20070101.nc'), 'r')

        # Check validity fields
        numpy.testing.assert_equal(results['tas'][0:160,:].mask, ref_sh.variables['tas'][0,:,:].mask)
        numpy.testing.assert_equal(results['tasmin'][0:160,:].mask, ref_sh.variables['tasmin'][0,:,:].mask)
        numpy.testing.assert_equal(results['tasmax'][0:160,:].mask, ref_sh.variables['tasmax'][0,:,:].mask)
        self.assertTrue(results['tas'][161:560,:].mask.all())
        self.assertTrue(results['tasmin'][161:560,:].mask.all())
        self.assertTrue(results['tasmax'][161:560,:].mask.all())
        numpy.testing.assert_equal(results['tas'][560:720,:].mask, ref_nh.variables['tas'][0,:,:].mask)
        numpy.testing.assert_equal(results['tasmin'][560:720,:].mask, ref_nh.variables['tasmin'][0,:,:].mask)
        numpy.testing.assert_equal(results['tasmax'][560:720,:].mask, ref_nh.variables['tasmax'][0,:,:].mask)

        # Do comparisons only where valid data
        # - these are valid entries in whole field
        tasvalid = numpy.nonzero(~results['tas'].mask)
        tasminvalid = numpy.nonzero(~results['tasmin'].mask)
        tasmaxvalid = numpy.nonzero(~results['tasmax'].mask)
        # - and these are filtered to include only southern hemisphere
        sh_tasvalid = (tasvalid[0][tasvalid[0] < 160], tasvalid[1][tasvalid[0] < 160])
        sh_tasminvalid = (tasminvalid[0][tasminvalid[0] < 160], tasminvalid[1][tasminvalid[0] < 160])
        sh_tasmaxvalid = (tasmaxvalid[0][tasmaxvalid[0] < 160], tasmaxvalid[1][tasmaxvalid[0] < 160])
        # - and the corresponding indices in reference data
        ref_sh_tasvalid = sh_tasvalid
        ref_sh_tasminvalid = sh_tasminvalid
        ref_sh_tasmaxvalid = sh_tasmaxvalid
        # - and these are filtered to include only northern hemisphere
        nh_tasvalid = (tasvalid[0][tasvalid[0] >= 560], tasvalid[1][tasvalid[0] >= 560])
        nh_tasminvalid = (tasminvalid[0][tasminvalid[0] >= 560], tasminvalid[1][tasminvalid[0] >= 560])
        nh_tasmaxvalid = (tasmaxvalid[0][tasmaxvalid[0] >= 560], tasmaxvalid[1][tasmaxvalid[0] >= 560])
        # - and the corresponding indices in reference data (need to offset longitudes)
        ref_nh_tasvalid = (nh_tasvalid[0]-560, nh_tasvalid[1])
        ref_nh_tasminvalid = (nh_tasminvalid[0]-560, nh_tasminvalid[1])
        ref_nh_tasmaxvalid = (nh_tasmaxvalid[0]-560, nh_tasmaxvalid[1])

        # Compare observation fields
        numpy.testing.assert_almost_equal(results['tas'][sh_tasvalid], ref_sh.variables['tas'][0,:,:].data[ref_sh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin'][sh_tasminvalid], ref_sh.variables['tasmin'][0,:,:].data[ref_sh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax'][sh_tasmaxvalid], ref_sh.variables['tasmax'][0,:,:].data[ref_sh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas'][nh_tasvalid], ref_nh.variables['tas'][0,:,:].data[ref_nh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin'][nh_tasminvalid], ref_nh.variables['tasmin'][0,:,:].data[ref_nh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax'][nh_tasmaxvalid], ref_nh.variables['tasmax'][0,:,:].data[ref_nh_tasmaxvalid], decimal=4)

        # Uncertainty fields should match masks
        numpy.testing.assert_equal(results['tas_unc_rand'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_rand'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_rand'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_sys'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_sys'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_sys'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_cloud'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_cloud'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_cloud'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_corr_local'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_corr_local'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_corr_local'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tas_unc_no_cloud'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasmin_unc_no_cloud'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmax_unc_no_cloud'].mask, results['tasmax'].mask)
        numpy.testing.assert_equal(results['tasuncertainty'].mask, results['tas'].mask)
        numpy.testing.assert_equal(results['tasminuncertainty'].mask, results['tasmin'].mask)
        numpy.testing.assert_equal(results['tasmaxuncertainty'].mask, results['tasmax'].mask)

        # Uncertainty fields should match where observations valid (southern hemisphere)
        numpy.testing.assert_almost_equal(results['tas_unc_rand'].data[sh_tasvalid], ref_sh.variables['RU'][0,:,:].data[ref_sh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_rand'].data[sh_tasminvalid], ref_sh.variables['RUmin'][0,:,:].data[ref_sh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_rand'].data[sh_tasmaxvalid], ref_sh.variables['RUmax'][0,:,:].data[ref_sh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_corr_local'].data[sh_tasvalid], ref_sh.variables['SSU'][0,:,:].data[ref_sh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_corr_local'].data[sh_tasminvalid], ref_sh.variables['SSUmin'][0,:,:].data[ref_sh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_corr_local'].data[sh_tasmaxvalid], ref_sh.variables['SSUmax'][0,:,:].data[ref_sh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_sys'].data[sh_tasvalid], ref_sh.variables['LSU'][0,:,:].data[ref_sh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_sys'].data[sh_tasminvalid], ref_sh.variables['LSUmin'][0,:,:].data[ref_sh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_sys'].data[sh_tasmaxvalid], ref_sh.variables['LSUmax'][0,:,:].data[ref_sh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_cloud'].data[sh_tasvalid], ref_sh.variables['CU'][0,:,:].data[ref_sh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_cloud'].data[sh_tasminvalid], ref_sh.variables['CUmin'][0,:,:].data[ref_sh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_cloud'].data[sh_tasmaxvalid], ref_sh.variables['CUmax'][0,:,:].data[ref_sh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_no_cloud'].data[sh_tasvalid], ref_sh.variables['TotU'][0,:,:].data[ref_sh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_no_cloud'].data[sh_tasminvalid], ref_sh.variables['TotUmin'][0,:,:].data[ref_sh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_no_cloud'].data[sh_tasmaxvalid], ref_sh.variables['TotUmax'][0,:,:].data[ref_sh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasuncertainty'].data[sh_tasvalid], ref_sh.variables['TotUC'][0,:,:].data[ref_sh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasminuncertainty'].data[sh_tasminvalid], ref_sh.variables['TotUCmin'][0,:,:].data[ref_sh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmaxuncertainty'].data[sh_tasmaxvalid], ref_sh.variables['TotUCmax'][0,:,:].data[ref_sh_tasmaxvalid], decimal=4)

        # Uncertainty fields should match where observations valid (northern hemisphere)
        numpy.testing.assert_almost_equal(results['tas_unc_rand'].data[nh_tasvalid], ref_nh.variables['RU'][0,:,:].data[ref_nh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_rand'].data[nh_tasminvalid], ref_nh.variables['RUmin'][0,:,:].data[ref_nh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_rand'].data[nh_tasmaxvalid], ref_nh.variables['RUmax'][0,:,:].data[ref_nh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_corr_local'].data[nh_tasvalid], ref_nh.variables['SSU'][0,:,:].data[ref_nh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_corr_local'].data[nh_tasminvalid], ref_nh.variables['SSUmin'][0,:,:].data[ref_nh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_corr_local'].data[nh_tasmaxvalid], ref_nh.variables['SSUmax'][0,:,:].data[ref_nh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_sys'].data[nh_tasvalid], ref_nh.variables['LSU'][0,:,:].data[ref_nh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_sys'].data[nh_tasminvalid], ref_nh.variables['LSUmin'][0,:,:].data[ref_nh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_sys'].data[nh_tasmaxvalid], ref_nh.variables['LSUmax'][0,:,:].data[ref_nh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_cloud'].data[nh_tasvalid], ref_nh.variables['CU'][0,:,:].data[ref_nh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_cloud'].data[nh_tasminvalid], ref_nh.variables['CUmin'][0,:,:].data[ref_nh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_cloud'].data[nh_tasmaxvalid], ref_nh.variables['CUmax'][0,:,:].data[ref_nh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tas_unc_no_cloud'].data[nh_tasvalid], ref_nh.variables['TotU'][0,:,:].data[ref_nh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmin_unc_no_cloud'].data[nh_tasminvalid], ref_nh.variables['TotUmin'][0,:,:].data[ref_nh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmax_unc_no_cloud'].data[nh_tasmaxvalid], ref_nh.variables['TotUmax'][0,:,:].data[ref_nh_tasmaxvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasuncertainty'].data[nh_tasvalid], ref_nh.variables['TotUC'][0,:,:].data[ref_nh_tasvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasminuncertainty'].data[nh_tasminvalid], ref_nh.variables['TotUCmin'][0,:,:].data[ref_nh_tasminvalid], decimal=4)
        numpy.testing.assert_almost_equal(results['tasmaxuncertainty'].data[nh_tasmaxvalid], ref_nh.variables['TotUCmax'][0,:,:].data[ref_nh_tasmaxvalid], decimal=4)
