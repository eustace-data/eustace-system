"""
In-situ land source
-------------------
Observation break points source interface to output of in-situ land preprocessing.
"""

import numpy
import netCDF4
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import ObservationsBreakPoints
from datetime import datetime
from eustace.timeutils.epoch import days_since_epoch
from eustace.outputformats import definitions

class ObservationBreakPointSourceInsituLand(ObservationSource):
    """Provide observation break points source interface to output of pre-processed in-situ land data.
       Optionally restrict to a single day."""

    # Map observation source keys to file variable definitions
    OBSERVATIONMAPS = {	'merged_break' :         {  'dimension':definitions.DIMENSION_NAME_MERGED_BREAK,
                                                    'break_time': definitions.MERGED_BREAK_TIME,
                                                    'break_station': definitions.MERGED_BREAK_STATION,
                                                    'break_likelihood': definitions.MERGED_BREAK_LIKELIHOOD,
                                                    'detection_feasibility': definitions.DETECTION_FEASIBILITY}}

    def __init__(self, filename, station_index=None):
        """Load and parse specified filename according to specified format instance."""

        super(ObservationSource, self).__init__()
        self.dataset = netCDF4.Dataset(filename, 'r')
        self.station_index = station_index

    def observables(self):
        """The names of variables estimated from this source."""
        return ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS.keys()

    def number_of_observations(self, observable):
        """Total number of observed breaking points."""

        num_break_points = self.dataset.dimensions[ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS[observable]['dimension']].size
        if self.station_index is None:
            return num_break_points
        else:
            filtered_num_break_points = numpy.where(self.dataset.variables[ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS[observable]['break_station'].name][:]==self.station_index)[0].size
            return filtered_num_break_points

    def observation_location_lookup(self):
        """NumPy array of 3D coordinates of latitude, longitude, reference time."""

        # station coordinates
        longitude = self.dataset.variables[definitions.LONGITUDE.name][:]
        latitude = self.dataset.variables[definitions.LATITUDE.name][:]

        # allow for longitudes expressed as [0,360) instead of [-180,180)
        wrapped_longitudes = numpy.nonzero(longitude >= numpy.float(180.0))
        longitude[wrapped_longitudes] -= 360.0

        # combined coordinate
        location = numpy.vstack( [ latitude, longitude ] )

        return location

    def observations(self,observable='merged_break'):
        """Retrieve observations break points for specified observable quantity."""

        # generate filtering indices if station_index is generated
        if self.station_index is not None:
            station_mask = self.dataset.variables[ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS[observable]['break_station'].name][:] == self.station_index

        break_likelihood = self.dataset.variables[ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS[observable]['break_likelihood'].name]

        if hasattr(break_likelihood[:],'mask'):
            likelihood_mask = ~break_likelihood[:].mask
        else:
            likelihood_mask = numpy.broadcast_to(numpy.array([True]),break_likelihood[:].shape)

        parameters = {}
        keys = ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS[observable].keys()
        
        # we do not need the information about the dimension
        keys.remove('dimension')

        # Retrieve measurement and variable info
        for key in keys:
            
            if key == 'break_station':
                variable = self.dataset.variables[ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS[observable][key].name]
            else:
                variable = self.dataset.variables[ObservationBreakPointSourceInsituLand.OBSERVATIONMAPS[observable][key].name]
            if hasattr(variable[:],'mask'):
                measurement = variable[:].data
            else:
                measurement = variable[:]

            if self.station_index is not None:
                if key == 'detection_feasibility':
                    measurement = numpy.array([measurement[self.station_index]])
                else:
                    final_mask = numpy.logical_and(station_mask, likelihood_mask)
                    measurement = measurement[final_mask]
            else:
                if (key != 'detection_feasibility'):
                    measurement =  measurement[likelihood_mask]

            parameters[key] = measurement
        return ObservationsBreakPoints(**parameters)

    def local_correlation_length_scale(self, observable):
        pass
