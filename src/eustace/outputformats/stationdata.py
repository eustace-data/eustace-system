"""Represent station data."""

__version__ = "$Revision: 550 $"
__author__ = "Joel R. Mitchelson"

import json
import numpy
import numpy.ma
import definitions


class StationDescription(object):
    """Static information describing a measurement station."""

    def __init__(self, name, latitude, longitude, elevation):
        self.name = name
        self.latitude = latitude
        self.longitude = longitude
        self.elevation = elevation

    def __repr__(self):
        return json.dumps(self.__dict__)


class StationDescriptionArray(object):
    """Array of station descriptions."""

    def __init__(self, stations):
        """Build array and create station indices."""
        self.stations = stations
        self.station_index = {station.name: index for (
            index, station) in enumerate(stations)}

    def get_number_of_stations(self):
        """Number of stations in this array."""
        return len(self.stations)

    def get_station_names(self):
        """Retrieve list of station names."""
        return [station.name for station in self.stations]

    def get_station_index(self, stationname):
        """Retrieve the corresponding one-based ID for the given name."""
        return self.station_index[stationname] if stationname in self.station_index else None


class StationTemperatureSample(object):
    """Represent one day of temperature measurements from one station."""


    def __init__(self, station_name, time, tasmin, tasmin_qc, tasmax, tasmax_qc):

        # Data class - allow many arguments in constructor
        # pylint: disable=too-many-arguments

        self.station_name = station_name
        self.time = time
        self.tasmin = tasmin
        self.tasmin_qc = tasmin_qc
        self.tasmax = tasmax
        self.tasmax_qc = tasmax_qc


class StationBreakDetectionQuality(object):
    """Represent the quality of break-detection."""

    NOT_POSSIBLE = 0
    POSSIBLE_BUT_UNRELIABLE = 1
    RELIABLE = 2

    def __init__(self, detection_time, tasmin_detection_qc, tasmax_detection_qc):
        self.detection_time = detection_time
        self.tasmin_detection_qc = tasmin_detection_qc
        self.tasmax_detection_qc = tasmax_detection_qc


class StationBreak(object):
    """Represent a break in station data."""

    def __init__(self, station_name, amplitude, time_bounds, time_affected_bounds):
        self.station_name = station_name
        self.amplitude = amplitude
        self.time_bounds = time_bounds
        self.time_affected_bounds = time_affected_bounds

    def get_station_name(self):
        """Retrieve name of station to which break applies."""
        return self.station_name


class MeasurementTime(object):
    """Utility class for finding unique measurement times."""

    # Measurement times within this tolerance are considered identical for the
    # purposes of assembling data onto a common time axis
    RESOLUTION = 0.0005

    def __init__(self, time):
        """Initialise to specified time value."""

        self.time = time

    def __hash__(self):
        """Use value rounded to integer as hash."""

        return int(self.time)

    def __eq__(self, other):
        """Equality if equal up to MeasurementTime.RESOLUTION."""

        return abs(other.time - self.time) < MeasurementTime.RESOLUTION


class StationData(object):
    """Represent temperature time series and break detection."""

    def __init__(self, descriptions):
        self.descriptions = descriptions
        self.temperature_samples = []
        self.break_detection_quality = []
        self.breaks = {definitions.TASMIN.name: [],
                       definitions.TASMAX.name: []}

    def add_temperature_sample(self, sample):
        """Add temperature measurement."""

        self.temperature_samples.append(sample)

    def add_break_detection_period(self, quality):
        """Add a period of break detection quality information."""

        self.break_detection_quality.append(quality)

    def add_break(self, variable, stationbreak):
        """Add break information in the specified variable."""

        self.breaks[variable.name].append(stationbreak)

    def get_unique_temperature_sample_times(self):
        """Array of unique times based on the MeasurementTime.RESOLUTION."""

        result = [m.time for m in set(
            [MeasurementTime(sample.time) for sample in self.temperature_samples])]
        result.sort()
        return result

    def get_temperature_samples(self, time_axis, variable):
        """Get a temperature variable along the specified time axis, as a NumPy float32 masked array.
           Any missing entries have the mask flag set."""

        # make result of dimensions (time points, number of stations) as 32-bit float with mask
        # initialise mask to 1 (meaning masked out) and data to fill value of
        # variable
        dimensions = (len(time_axis),
                      self.descriptions.get_number_of_stations())
        data = numpy.empty(dimensions, numpy.float32)
        mask = numpy.empty(dimensions, numpy.bool_)
        data.fill(variable.fill_value)
        mask.fill(True)

        # use measurement time class for comparison of time values
        measurement_times = [MeasurementTime(t) for t in time_axis]

        # extract measurements along time axis and put into result
        # - set mask to zero where found
        for sample in self.temperature_samples:
            station_index = self.descriptions.get_station_index(
                sample.station_name)
            time_index = measurement_times.index(MeasurementTime(sample.time))
            value = sample.__dict__[variable.name]
            data[time_index][station_index] = value
            mask[time_index][station_index] = False

        # form masked result
        return numpy.ma.array(data, mask=mask)

    def get_break_detection_quality_variable(self, variable):
        """Get one aspect of break detection quality as a NumPy array."""

        return numpy.array([quality.__dict__[variable.name] for quality in self.break_detection_quality], variable.dtype)

    def get_break_variable(self, variable):
        """Get one aspect of breaks as a NumPy array."""

        # A variable name like tasmin_break_amplitude gives us
        # a measurement name of 'tasmin' and paramter name of 'amplitude'
        measurement_name, parameter_name = variable.name.split('_break_')

        # stations are a special case - must retrieve as a name first
        # then convert to indices
        if parameter_name == 'station':
            parameter_name = 'station_name'

        # get values
        values = [breakdata.__dict__[parameter_name]
                  for breakdata in self.breaks[measurement_name]]

        # if station names were retrieved, convert now to station indices
        if parameter_name == 'station_name':
            values = [self.descriptions.get_station_index(
                name) for name in values]

        # convert to NumPy type and return
        return numpy.array(values, variable.dtype)
