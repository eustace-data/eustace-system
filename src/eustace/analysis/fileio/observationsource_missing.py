"""A dummy observation source for use is no data is available"""

import numpy
from ..observationsource import ObservationSource
from ..observationsource import Observations

class ObservationSourceMissing(ObservationSource):
    """Define an observation source."""

    def __init__(self, observables):
        """Empty constructor in bass class"""
        
        self.missing_observables = observables

    def observables(self):
        """The names of variables estimated from this source."""
        return self.missing_observables

    def observation_location_lookup(self):
        """NumPy array in which column number corresponds to location id and rows are latitude and longitude."""
        return numpy.zeros( (0,0) )

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""
        
        if self.missing_observable == observable:
            # return as Observations object
            return ObservationsMissing(None)
        else:
            raise ValueError('Unknown observable: '+observable)

    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty."""
        return numpy.zeros( (0,0) )

class ObservationsMissing(Observations):
    
    def __init__(self, time):
        #"""Build observation set from NumPy arrays."""
        self.mask = numpy.ones(shape=(0), dtype=numpy.bool)
        self.time = time
        self.location = numpy.empty(shape=(0))
        self.measurement = numpy.empty(shape=(0))
        self.uncorrelatederror = numpy.empty(shape=(0))
        self.locallycorrelatederror = numpy.empty(shape=(0))

    def number_of_observations(self):
        
        return 0