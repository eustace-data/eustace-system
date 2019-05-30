"""
Surface-air ice data
------------------------
Read ice surface input apply quality control ready for input to model.
Read surface-air model output for ice surface in original WP1 format.
"""

from cubemap import ObservationSourceCubeMapGrid
from cubemap import ObservationMap
from eustace.analysis.observationsource import ObservationSource
from eustace.timeutils.epoch import days_since_epoch
from ..model_ice_input_processing import IceSurfaceTemperatureQualityControl
from ..model_ice import ModelIce
import numpy
import netCDF4
from datetime import datetime

class ObservationSourceIceDMI(ObservationSourceCubeMapGrid):
    """Provide ObservationSource interface to ice surface data loaded from Iris."""

    OBSERVATIONMAPS = { 
        ObservationSource.TMEAN: ObservationMap('tas', 'RU', [ 'SSU' ], [ 500000.0 ]),
        ObservationSource.TMAX: ObservationMap('tasmax', 'RUmax', [ 'SSUmax' ], [ 500000.0 ]),
        ObservationSource.TMIN: ObservationMap('tasmin', 'RUmin', [ 'SSUmin' ], [ 500000.0 ])
    }

    def __init__(self, filename):
        super(ObservationSourceIceDMI, self).__init__(ObservationSourceIceDMI.OBSERVATIONMAPS, filename)


class ObservationSourceIceDMICombined(ObservationSource):
    """Provide ObservationSource interface to a pair of files (for northern and southern hemisphere combined)."""

    OBSERVATIONMAPS = ObservationSourceIceDMI.OBSERVATIONMAPS

    def __init__(self, filenames):
        """Load with list of filenames (should be two for NH and SH)."""

        if len(filenames) != 2:
            raise ValueError('ObservationSourceIceDMICombined expects exactly two input filenames (got {0})'.format(len(filenames)))

        # Get individual sources
        self.sources = [ ObservationSourceIceDMI(filename) for filename in filenames ]

        # Retrieve measurements item for each source
        measurements = [ source.get_field_by_name(ObservationSourceIceDMICombined.OBSERVATIONMAPS[ObservationSource.TMEAN].measurement) for source in self.sources ]
        longitudes = [ m.coords('longitude')[0].points.ravel() for m in measurements ]
        latitudes = [ m.coords('latitude')[0].points.ravel() for m in measurements ]

        # Longitudes should simply match
        diff_longitude = numpy.max(numpy.abs(longitudes[1] - longitudes[0]))
        if diff_longitude > 0.000001:
            raise ValueError('ObservationSourceIceDMICombined longitudes mismatch')

        # Latitudes should be at same resolution for both (and same sign of direction of increase)
        latitude_resolution = latitudes[0][1] - latitudes[0][0]
        if latitude_resolution < 0.0:
            raise ValueError('ObservationSourceIceDMICombined expect latitudes increasing (but are decreasing)')
        if latitudes[1][0] <= latitudes[0][-1]:
            raise ValueError('ObservationSourceIceDMICombined expect second file to have latitudes greater than first (but does not)')
        diff_latitude0 = numpy.max(numpy.abs((latitudes[0][1:] - latitudes[0][0:-1]) - latitude_resolution))
        if diff_latitude0 > 0.000001:
            raise ValueError('ObservationSourceIceDMICombined latitude resolution non-constant')
        diff_latitude1 = numpy.max(numpy.abs((latitudes[1][1:] - latitudes[1][0:-1]) - latitude_resolution))
        if diff_latitude1 > 0.000001:
            raise ValueError('ObservationSourceIceDMICombined latitude resolution mismatch')

        # These are the axes to be used for lookup
        self.longitude_axis_points = longitudes[0]
        self.latitude_axis_points = numpy.arange(latitudes[0][0], latitudes[1][-1] + 0.000001, step=latitude_resolution)

        # The difference in latitude increments between the last latitude in the first file and the first latitude in the second file
        latitude_increments_between_files = int( ((latitudes[1][0] - latitudes[0][-1]) / latitude_resolution) + 0.000001 )

        # The offset that must be added to locations for new lookup
        self.locationoffsets = [ 0, (latitudes[0].shape[0] - 1 + latitude_increments_between_files) * self.longitude_axis_points.shape[0] ]


    def observables(self):
        """The names of variables estimated from this source."""
        return ObservationSourceIceDMICombined.OBSERVATIONMAPS.keys()

    def number_of_observations(self):
        """Total number of observations."""
        
        # Sum over multiple sources
        return sum( [ source.number_of_observations() for source in self.sources ] )

    def observation_location_lookup(self):
        """NumPy array in which column number corresponds to location id (grid box index in this case),
           and rows are latitude and longitude."""

        mesh = numpy.meshgrid(self.latitude_axis_points, self.longitude_axis_points, indexing='ij')    
        return numpy.vstack([mesh[0].ravel(), mesh[1].ravel()])


    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""

        # Get first source obs
        obs = self.sources[0].observations(observable)

        # Append any more
        for index in range(1, len(self.sources)):
            obs.append(self.sources[index].observations(observable), self.locationoffsets[index])

        return obs

    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty."""

        return ObservationSourceIceDMICombined.OBSERVATIONMAPS[observable].localcorrelationscale


class IceSurfaceTemperatureQualityControlNetCDF(IceSurfaceTemperatureQualityControl):
    """Wrapper to read input to ice model via quality control class."""

    # Definition of known field names in NetCDF (also match the constructor method params in base class)
    FIELDNAMES = [ 'lat',
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
    
    def __init__(self, filename):
        """Load from filename."""

        # Load file
        inputdata = netCDF4.Dataset(filename, 'r')

        # Retrieve fields
        inputfields = { name: inputdata.variables[name][:] for name in IceSurfaceTemperatureQualityControlNetCDF.FIELDNAMES }

        # Reduce to 2D those with time dimension
        inputfields = { name: field[0,:,:] if (field.ndim == 3 and not '3d' in name) else field for name, field in inputfields.iteritems() }            

        # Get time
        timevariable = inputdata.variables['time']
        fieldtime = netCDF4.num2date(timevariable[0], units=timevariable.units)

        # Convert to days since EUSTACE epoch
        daynumber = int(days_since_epoch(fieldtime))

        # Construct
        super(IceSurfaceTemperatureQualityControlNetCDF, self).__init__(daynumber, **inputfields)

