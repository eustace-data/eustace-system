"""Build EUSTACE-compliant station data files."""

__version__ = "$Revision: 552 $"
__author__ = "Joel R. Mitchelson"

import os
import numpy
from filebuilder import FileBuilder
from eustace.outputformats import definitions
from eustace.outputformats.outputvariable import OutputVariable
from eustace.timeutils.epoch import days_since_epoch

class FileBuilderStationData(object):
    """Build a NetCDF file with EUSTACE specification for station data."""

    @staticmethod
    def create_component_builder(filenamepattern, titlesuffix, directory, region, version, institution, comment, history, source):
        """Helper to create component file."""
        pathname = os.path.join(directory, filenamepattern.format(
            region=region.lower(), version=version))
        title = definitions.TITLE_PREFIX + titlesuffix
        builder = FileBuilder()
        builder.create(pathname, title, institution, comment, history, source)
        return builder

    def __init__(self, directory, region, version, institution, comment, history, source):
        """Create with specified attributes."""

        super(FileBuilderStationData, self).__init__()

        self.temperature = FileBuilderStationData.create_component_builder(
            definitions.STATIONS_TEMPERATURE_FILENAME_PATTERN, definitions.STATIONS_TEMPERATURE_TITLE_SUFFIX,
            directory, region, version, institution, comment, history, source)

        self.status = FileBuilderStationData.create_component_builder(
            definitions.STATIONS_STATUS_FILENAME_PATTERN, definitions.STATIONS_STATUS_TITLE_SUFFIX,
            directory, region, version, institution, comment, history, source)

        # flag temperature file as time series
        self.temperature.nc.__setattr__(definitions.STATIONS_TEMPERATURE_FEATURE_TYPE_NAME, definitions.STATIONS_TEMPERATURE_FEATURE_TYPE_VALUE)

        # always have a 2D dimension for time ranges
        self.add_status_dimension(definitions.DIMENSION_NAME_BOUNDS, 2)

        # empty stations item
        self.stations = None

    def set_temperature_dimensions(self, stations, time_axis):
        """Create station and time axis used for temperature variables. Must do this before adding temperature variables."""

        # store station info
        self.stations = stations

        # Station dimension (known number of stations)
        self.temperature.add_dimension(definitions.DIMENSION_NAME_STATION, self.stations.get_number_of_stations())

        # Length of strings for station names
        self.temperature.add_dimension(definitions.DIMENSION_NAME_NAME_STRLEN, definitions.STATIONS_NAME_STRLEN)

        # NumPy character array of station names
        station_names = self.stations.get_station_names()
        numpy_station_names = numpy.empty(
            (len(station_names), definitions.STATIONS_NAME_STRLEN), definitions.STATION_NAME.dtype)
        for station_index, station_name in enumerate(station_names):
            numpy_station_name = numpy.fromstring(
                station_name, dtype=definitions.STATION_NAME.dtype)
            numpy_station_name.resize((definitions.STATIONS_NAME_STRLEN,))
            numpy_station_names[station_index] = numpy_station_name

        # set station names
        self.temperature.add_variable(definitions.STATION_NAME, [definitions.DIMENSION_NAME_STATION, definitions.DIMENSION_NAME_NAME_STRLEN], numpy_station_names)

        # Make time dimension of unlimited extent
        self.temperature.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_TIME,
            variable=definitions.TIME,
            values=time_axis,
            unlimited=True)

        # always have time bounds too (1 day)
        self.temperature.add_dimension(definitions.DIMENSION_NAME_BOUNDS, 2)
        self.temperature.add_variable(definitions.TIME_BOUNDS, [ definitions.DIMENSION_NAME_TIME, definitions.DIMENSION_NAME_BOUNDS ],
            numpy.transpose(numpy.vstack((time_axis, time_axis+1.0))))

    def add_temperature_variable(self, variable, values):
        """Add the specififed variable with given initial values."""

        # Temperature data dependent on time, with values for all stations
        # given at each time
        self.temperature.add_variable(variable, [definitions.DIMENSION_NAME_TIME, definitions.DIMENSION_NAME_STATION], values)

    def add_status_dimension(self, index_dimension_name, size=None):
        """Add a dimension to the status file (such as break detection index)."""

        self.status.add_dimension(index_dimension_name, size)

    def add_status_variable(self, index_dimension_name, variable, values, bounds=False):
        """Add status variable (dimension must already exist)."""

        # Always include index dimension (like break index)
        dimension_names = [index_dimension_name]

        # Optionally add bounds
        if bounds:
            dimension_names.append(definitions.DIMENSION_NAME_BOUNDS)

        # Create
        self.status.add_variable(variable, dimension_names, values)

    def save_and_close(self):
        """Close the files and store results."""

        self.temperature.save_and_close()
        self.status.save_and_close()
