"""
Observation source interface to iris cube structure
---------------------------------------------------
Map an iris cube as an observation source.
"""

from iris import load
from iris import Constraint
from iris.cube import CubeList
import numpy

from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations

class ObservationMap(object):
    """Character strings defining the mapping to file."""

    def __init__(self, measurement, uncorrelatederror, locallycorrelatederror, localcorrelationscale):
        self.measurement = measurement
        self.uncorrelatederror = uncorrelatederror
        self.locallycorrelatederror = locallycorrelatederror
        self.localcorrelationscale = localcorrelationscale

class ObservationSourceCubeMapGrid(ObservationSource, CubeList):
    """Provide ObservationSource interface to ice surface data loaded from Iris."""

    def __new__(cls, *args, **kwargs):
        """"This is needed because CubeList overrides __new__ without allowing for variable args and kwargs,
            so we have problems with subclassing unless we explicitly allow it again here."""

        return CubeList.__new__(cls)

    def __init__(self, obsmaps, filename):
        """Construct with specified observation maps."""

        super(ObservationSourceCubeMapGrid, self).__init__()
        self.obsmaps = obsmaps
        self.extend( load(filename) )

    def observables(self):
        """The names of variables estimated from this source."""
        return self.obsmaps.keys()

    def number_of_observations(self):
        """Total number of observations."""

        measurement = self.get_field_by_name(self.obsmaps.values()[0].measurement)
        return numpy.prod(measurement.shape)

    def observation_location_lookup(self):
        """NumPy array in which column number corresponds to location id (grid box index in this case),
           and rows are latitude and longitude."""

        measurement = self.get_field_by_name(self.obsmaps.values()[0].measurement)
        longitude = measurement.coords('longitude')
        latitude = measurement.coords('latitude')
        mesh = numpy.meshgrid(numpy.float64(latitude[0].points.ravel()), numpy.float64(longitude[0].points.ravel()), indexing='ij')
        return numpy.vstack([mesh[0].ravel(), mesh[1].ravel()])

    def get_field_by_name(self, searchname):
        """Helper to retrieve field using var_name. Raises ValueError if not found."""

        results =  self.extract(Constraint(cube_func=lambda cube: cube.var_name==searchname))
        if len(results) != 1:
            raise ValueError('No unique field found with var_name: {0}'.format(searchname))
        return results[0]

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""

        obsmap = self.obsmaps[observable]
        measurement = self.get_field_by_name(obsmap.measurement)
        location = numpy.array(range(self.number_of_observations()), numpy.uint64) # location ids are just grid indices
        time = measurement.coords('time')[0].points[0]
        uncorrelatederror = self.get_field_by_name(obsmap.uncorrelatederror).data
        locallycorrelatederror = [ self.get_field_by_name(mapname).data.data.ravel() for mapname in obsmap.locallycorrelatederror ]
        localcorrelationscale = obsmap.localcorrelationscale
        return Observations(
            measurement.data.mask.ravel(), 
            time, 
            location, 
            measurement.data.data.ravel(), 
            uncorrelatederror.data.ravel(), 
            locallycorrelatederror)

    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty."""

        return self.obsmaps[observable].localcorrelationscale
