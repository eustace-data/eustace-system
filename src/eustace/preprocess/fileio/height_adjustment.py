"""Adjust observations to a standard altitude surface"""

import numpy

from netCDF4 import Dataset

from eustace.analysis.advanced_standard.elements.geography_based import GeographyBasedCovariateFunction

LAPSE_RATE = -6.5            # Lapse rate used for the adjustment. Approximately the environmental lapse rate in K km^-1.
LAPSE_RATE_UNCERTAINTY = 1.1 # 1 sigma uncertainty in lapse rate in K km^-1.

class AltitudeAdjustment(object):
    """Lapse rate based altitude adjustment to an interpolated reference surface."""

    def __init__(self,  filename, latitude_label, longitude_label, covariate_label, rescale_factor):
        """Partially initialise the class with (non) optional argument to prevent EUMOPPS from trying store the elevation map in """
        self.filename = filename
        self.latitude_label = latitude_label
        self.longitude_label = longitude_label
        self.covariate_label = covariate_label
        self.rescale_factor = rescale_factor

    def load_reference_map(self):
        """Load the reference elevation map data from a specific covariate file"""

        covariate_file = Dataset(self.filename)
        self.latitude = covariate_file.variables[self.latitude_label][:]
        self.longitude = covariate_file.variables[self.longitude_label][:]
        self.covariate = covariate_file.variables[self.covariate_label][:]
        covariate_file.close()

    def reference_altitudes(self, locations):
        """Get reference altitudes by interpolation of elevation map"""
        return GeographyBasedCovariateFunction(self.latitude, self.longitude, self.covariate, self.rescale_factor).compute(locations).ravel()
        
    def adjustment(self, altitudes, locations):
        """Compute the additive adjustment factor by scaling the altitude difference from reference by the negative of the lapse rate"""        
        return -LAPSE_RATE * (altitudes - self.reference_altitudes(locations))
        
    def adjustment_uncertainty(self, altitudes, locations):
        """Compute the additive adjustment uncertainty by scaling the altitude difference from reference by the lapse rate"""        
        return LAPSE_RATE_UNCERTAINTY * numpy.abs((altitudes - self.reference_altitudes(locations)))


