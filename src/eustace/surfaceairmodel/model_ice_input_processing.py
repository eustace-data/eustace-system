# Ice model pre-processing for quality control
# Ported from DMI Matlab code

import numpy
import scipy.signal
import scipy.io
from model_ice import ModelIce

class IceSurfaceTemperatureQualityControl(object):
    """Apply quality control to input surface air temperatures and make results available for processing."""

    # Quality control parameters
    MIN_SEA_ICE_CONC = 0.3
    MIZ_UPPER_THRESHOLD = 0.85
    TEMPERATURE_STANDARD_DEVIATION_MAX = 7.07  # ( expected daily range up to 20C, std(10*sin(2pi*t/T)) )
    TEMPERATURE_MAX = 278.15 # Exclude ice temperatures measured at greater than 5C
    NB_TOL = 10
    COMBINED_BIN_MAX = 10.0
    DAY_BIN_MAX = 20.0
    NIGHT_BIN_MIN = -20.0
    CLOUD_MASK_HALFWIDTH = 2
    CLOUD_MASK_MIN_COUNT = 2

    # Definition of expected input bins
    NIGHTBINS = [ 0, 1, 6, 7 ]
    DAYBINS = [ 2, 3, 4, 5 ]

    def __init__(self,
                 daynumber,
                 lat,
                 lon,
                 surface_temperature,
                 surface_temperature_min,
                 surface_temperature_max,
                 surface_temperature_std,
                 surface_temperature_3d,
                 Nobs3d,
                 land_mask,
                 sea_ice_fraction,
                 uncorrelated_uncertainty,
                 synoptically_correlated_uncertainty,
                 large_scale_correlated_uncertainty):
        """Construct from relevant NetCDF variables."""

        # Check expected dimensions
        if lat.ndim != 1:
            raise ValueError('latitude not 1-dmensional')
        if lon.ndim != 1:
            raise ValueError('longitude not 1-dmensional')
        if surface_temperature.ndim != 2:
            raise ValueError('surface_temperature not 2-dmensional')
        if Nobs3d.ndim != 3:
            raise ValueError('Nobs3d not 3-dmensional')
        if surface_temperature_3d.ndim != 3:
            raise ValueError('surface_temperature_3d not 3-dmensional')

        # Check sizes are consistent
        if surface_temperature.shape != (lat.shape[0], lon.shape[0]):
            raise ValueError('surface_temperature_min matrix shape not consistent with latitude and longitude dimensions')
        if surface_temperature.shape != surface_temperature_min.shape:
            raise ValueError('surface_temperature_min matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != surface_temperature_max.shape:
            raise ValueError('surface_temperature_max matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != surface_temperature_std.shape:
            raise ValueError('surface_temperature_std matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != surface_temperature_3d.shape[1:3]:
            raise ValueError('surface_temperature_3d matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != Nobs3d.shape[1:3]:
            raise ValueError('Nobs3d matrix shape not consistent with surface_temperature')
        if land_mask.shape != land_mask.shape:
            raise ValueError('land_mask matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != sea_ice_fraction.shape:
            raise ValueError('sea_ice_fraction matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != uncorrelated_uncertainty.shape:
            raise ValueError('uncorrelated_uncertainty matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != large_scale_correlated_uncertainty.shape:
            raise ValueError(' matrix shape not consistent with surface_temperature')
        if surface_temperature.shape != synoptically_correlated_uncertainty.shape:
            raise ValueError(' matrix shape not consistent with surface_temperature')

        # Store
        self.daynumber = daynumber
        self.latitude = lat
        self.longitude = lon
        self.surface_temperature = surface_temperature
        self.surface_temperature_min = surface_temperature_min
        self.surface_temperature_max = surface_temperature_max
        self.surface_temperature_3d = surface_temperature_3d
        self.surface_temperature_std = surface_temperature_std
        self.Nobs3d = Nobs3d
        self.land_mask = land_mask
        self.sea_ice_fraction = sea_ice_fraction
        self.uncorrelated_uncertainty = uncorrelated_uncertainty
        self.synoptically_correlated_uncertainty = synoptically_correlated_uncertainty
        self.large_scale_correlated_uncertainty = large_scale_correlated_uncertainty

    def hemisphere(self):
        """Use latitude information to determine hemisphere."""

        return ModelIce.NORTHERN_HEMISPHERE if self.latitude[0] > 0.0 else ModelIce.SOUTHERN_HEMISPHERE

    def land_ice_mask(self):
        """2D mask which is True when there is no land ice, false otherwise."""

        return ((self.land_mask != 1) & (self.land_mask != 4))

    def sea_ice_concentration_mask(self):
        """2D mask which is True when ice conentration is unknown or less than minimum sea ice threshold, False otherwise."""

        return (self.sea_ice_fraction < IceSurfaceTemperatureQualityControl.MIN_SEA_ICE_CONC).filled(True)

    def sea_ice_mask(self):
        """2D mask which is True when over land or when ice conentration is less than minimum sea ice threshold, False otherwise."""

        return self.sea_ice_concentration_mask() | ~self.land_ice_mask()

    def ice_mask(self):
        """2D mask which is True when neither land nor sea ice is present, False otherwise."""

        return self.land_ice_mask() & self.sea_ice_concentration_mask()

    def marginal_ice_zone_mask(self):
        """2D mask which is False where ice concentration is in in the marginal ice zone, True when outside it."""

        return self.sea_ice_mask() | (self.sea_ice_fraction >= IceSurfaceTemperatureQualityControl.MIZ_UPPER_THRESHOLD)

    def temperature_standard_deviation_mask(self):
        """2D mask which is True when standard deviation is defined and out of range, False when in range or unknown."""
        
        return (self.surface_temperature_std > IceSurfaceTemperatureQualityControl.TEMPERATURE_STANDARD_DEVIATION_MAX).filled(False)

    def temperature_value_mask(self):
        """2D mask which is True when temperature is out of range, False when in range or unknown."""

        return (self.surface_temperature > IceSurfaceTemperatureQualityControl.TEMPERATURE_MAX).filled(False)

    def observations_day_mask(self):
        """2D mask which is True where there are zero daytime observations  or observation count unknown."""

        num_obs_day = numpy.sum(self.Nobs3d[IceSurfaceTemperatureQualityControl.DAYBINS,:,:], axis=0)
        return (num_obs_day == 0).filled(True)

    def observations_night_mask(self):
        """2D mask which is True where there are zero nighttime observations or observation count unknown."""

        num_obs_night = numpy.sum(self.Nobs3d[IceSurfaceTemperatureQualityControl.NIGHTBINS,:,:], axis=0)
        return (num_obs_night == 0).filled(True)

    def observations_day_or_night_mask(self):
        """2D mask which is True if there are zero daytime observations or zero nighttime observations or counts are unknown."""
        
        return self.observations_day_mask() | self.observations_night_mask()

    def combined_bin_comparison_mask(self):
        """2D mask which is True if surface temperature far from mean bin value, False when in range or unknown."""

        bin_mean_all = self.surface_temperature_3d.mean(axis=0)
        comparison = numpy.abs(self.surface_temperature - bin_mean_all) > IceSurfaceTemperatureQualityControl.COMBINED_BIN_MAX
        return comparison.filled(False)

    def day_bin_comparison_mask(self):
        """2D mask which is True if maximum surface temperature far from mean daytime bin value, False when in range or unknown."""

        bin_mean_day = numpy.mean(self.surface_temperature_3d[IceSurfaceTemperatureQualityControl.DAYBINS,:,:], axis=0)
        comparison = numpy.abs(self.surface_temperature_max - bin_mean_day) > IceSurfaceTemperatureQualityControl.DAY_BIN_MAX
        return comparison.filled(False)

    def night_bin_comparison_mask(self):
        """2D mask which is True if minimum surface temperature far from mean nighttime bin value."""

        bin_mean_night = numpy.mean(self.surface_temperature_3d[IceSurfaceTemperatureQualityControl.NIGHTBINS,:,:], axis=0)
        comparison = numpy.abs(self.surface_temperature_min - bin_mean_night) < IceSurfaceTemperatureQualityControl.NIGHT_BIN_MIN
        return comparison.filled(False)

    def filter_outliers(self, variable, threshold_mask):
        """Cloud outlier detection."""

        # Cloud mask from outlier detection
        sea_ice_outliers = IceSurfaceTemperatureQualityControl.cloud_mask(variable, threshold_mask | self.sea_ice_concentration_mask())
        land_ice_outliers = IceSurfaceTemperatureQualityControl.cloud_mask(variable, threshold_mask | self.land_ice_mask())

        # Return result
        return numpy.ma.masked_array(data=variable.data, mask=(variable.mask | self.ice_mask() | sea_ice_outliers | land_ice_outliers | threshold_mask))

    def quality_controlled_surface_temperature(self):
        """Combined mask for surface temperature, False where valid, True otherwise."""

        # Mask from data threshold logic
        threshold_mask = \
            self.temperature_standard_deviation_mask() | \
            self.temperature_value_mask() | \
            self.observations_day_or_night_mask() | \
            self.combined_bin_comparison_mask()

        return self.filter_outliers(self.surface_temperature, threshold_mask)

    def quality_controlled_surface_temperature_min(self):
        """Combined mask for surface temperature, False where valid, True otherwise."""

        # Mask from data threshold logic
        threshold_mask = \
            self.temperature_standard_deviation_mask() | \
            self.temperature_value_mask() | \
            self.observations_night_mask() | \
            self.night_bin_comparison_mask()

        return self.filter_outliers(self.surface_temperature_min, threshold_mask)

    def quality_controlled_surface_temperature_max(self):
        """Combined mask for surface temperature, False where valid, True otherwise."""

        # Mask from data threshold logic
        threshold_mask = \
            self.temperature_standard_deviation_mask() | \
            self.temperature_value_mask() | \
            self.observations_day_mask() | \
            self.day_bin_comparison_mask()

        return self.filter_outliers(self.surface_temperature_max, threshold_mask)

    def uncertainty_uncorrelated(self):
        """Uncorrelated uncertainty component."""

        return self.uncorrelated_uncertainty

    def uncertainty_synoptically_correlated(self):
        """Locally correlated uncertainty component 0 (synoptic scale)."""

        return self.synoptically_correlated_uncertainty

    def uncertainty_large_scale_correlated(self):
        """Locally correlated uncertainty component 1 (large scale)."""

        return self.large_scale_correlated_uncertainty

    @staticmethod
    def cloud_mask(variable, icemask):

        # Apply icemask to variable
        variable_remasked = numpy.ma.array(data=variable.data, mask=(variable.mask | icemask))

        # Parameters to apply
        halfwidth = IceSurfaceTemperatureQualityControl.CLOUD_MASK_HALFWIDTH
        mincount = IceSurfaceTemperatureQualityControl.CLOUD_MASK_MIN_COUNT

        # Compute local mean
        # (masked array which is masked where fewer than MIN_COUNT pixels valid in local region)
        localmean = IceSurfaceTemperatureQualityControl.local_mean(halfwidth, mincount, variable_remasked)

        # Evaluate how much less than local mean this is (where we have MIN_COUNT information to tell us)
        amount_below_mean = localmean - variable_remasked[halfwidth:-halfwidth,halfwidth:-halfwidth]

        # Threshold this for those pixels where we have info
        central_result = (amount_below_mean.data > IceSurfaceTemperatureQualityControl.NB_TOL) & ~amount_below_mean.mask

        # Pad to make same size as original input
        expanded_result = numpy.pad(central_result, ((halfwidth,halfwidth), (halfwidth,halfwidth)), mode='constant', constant_values=False)

        return expanded_result

    @staticmethod
    def local_mean(halfwidth, mincount, x):
        """
        Compute local mean of masked array x over square regions of side length 2.halfwidth + 1.

        Returns a masked array. Elements with local count fewer than mincount are masked in output.
        """

        # Put zeros where mask is True
        xfilled = x.filled(0.0)

        # Kernel for sum is just ones over square of required size but zero in centre
        kernel = numpy.ones((2*halfwidth+1, 2*halfwidth+1), dtype=numpy.int64)
        kernel[halfwidth,halfwidth] = 0

        # Compute the sum of unmasked elements over kernel area
        # - will return central portion only
        s = scipy.signal.convolve(xfilled, kernel, 'valid')

        # Count unmasked elements over kernel area
        n = scipy.signal.convolve(x.mask == False, kernel, 'valid')

        # Masked version of s using threshold on n
        s_masked = numpy.ma.masked_array(data=s, mask=(n < mincount))

        # Result
        return s_masked / n
