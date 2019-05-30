"""Example output."""

from eustace.outputformats import definitions
from eustace.outputformats.filebuilder import FileBuilder
from eustace.outputformats.outputvariable import OutputVariable
from eustace.timeutils.epoch import days_since_epoch
import numpy
from netCDF4 import default_fillvals

TAS_ANOMALY = OutputVariable(
    name='tas',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    standard_name='air_temperature',
    long_name='near_surface_temperature_anomaly',
    units='K')

class FileBuilderHadCRUT4ExampleOutput(FileBuilder):
    """Build a NetCDF file a bit like EUSTACE output format."""

    def __init__(self, pathname, outputstructure):
        """Create with specified output structure."""

        super(FileBuilder, self).__init__()

        # get EUSTACE daynumber
        daynumber = days_since_epoch(outputstructure.time_datetime())

        # make the file
        self.create(pathname, title='Test output', institution='', comment='Test output', history='', source='')

        # time variable
        time_variable = OutputVariable(
            name='time',
            dtype=numpy.float32,
            fill_value=None,
            standard_name='time',
            long_name='Time',
            units=definitions.TIME_UNITS_DAYS_SINCE_EPOCH,
            calendar='gregorian',
            axis='T')

        # set time in days since the EUSTACE epoch [and set UNLIMITED]
        self.add_dimension_and_variable(
            dimensionname='time',
            variable=time_variable,
            values=numpy.array([daynumber], numpy.float32),
            unlimited=True)

        # global latitude axis
        self.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_LATITUDE,
            variable=definitions.LATITUDE,
            values=outputstructure.latitudes)

        # global longitude axis
        self.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_LONGITUDE,
            variable=definitions.LONGITUDE,
            values=outputstructure.longitudes)

    def add_global_field(self, variable, values):
        """Add a field (must be a global field)."""

        self.add_variable(
            variable,
            (definitions.DIMENSION_NAME_TIME, definitions.DIMENSION_NAME_LATITUDE, definitions.DIMENSION_NAME_LONGITUDE),
            values)

    def add_uncertainty_parameter(self, name, long_name, values):
        """Add uncertainty parameter (must have global coverage)."""

        self.add_variable(
            OutputVariable(
                name=name,
                dtype=numpy.int16,
                fill_value=definitions.TEMPERATURE_UNCERTAINTY_FILL_VALUE,
                scale_factor=definitions.TEMPERATURE_UNCERTAINTY_SCALE_FACTOR,
                add_offset=definitions.TEMPERATURE_UNCERTAINTY_ADD_OFFSET,
                long_name=long_name,
                units='K'),
            (definitions.DIMENSION_NAME_TIME, definitions.DIMENSION_NAME_LATITUDE, definitions.DIMENSION_NAME_LONGITUDE),
            values)


