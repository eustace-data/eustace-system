"""Specification for global field verification."""

from specification import VariableSpecification
from specification import FileSpecification
import numpy


class SpecificationGlobalField(FileSpecification):
    """Define the specification against which global field files are checked."""

    def __init__(self):

        # Initialise attributes and empty variable/dimension structures
        super(SpecificationGlobalField, self).__init__()

        # Set required dimensions
        self.dimensions['latitude'] = 720
        self.dimensions['longitude'] = 1440
        self.dimensions['time'] = None
        self.dimensions['bounds'] = 2

        # Set variables

        self.variables.append(
            VariableSpecification('time', numpy.float32)
            .add_metadata('standard_name', 'time')
            .add_metadata('long_name', 'Time at zero longitude')
            .add_metadata('units', 'days since 1850-01-01T00:00:00Z')
            .add_metadata('calendar', 'gregorian')
            .add_metadata('axis', 'T')
            .add_metadata('bounds', 'timebounds')
            .add_metadata('ancillary_variables', 'timeoffset'))

        self.variables.append(
            VariableSpecification('timebounds', numpy.float32))

        self.variables.append(
            VariableSpecification('timeoffset', numpy.float32)
            .add_metadata('long_name', 'Local time offset from UTC (days)')
            .add_metadata('units', 'days'))

        self.variables.append(
            VariableSpecification('latitude', numpy.float32)
            .add_metadata('standard_name', 'latitude')
            .add_metadata('long_name', 'Latitude (deg)')
            .add_metadata('units', 'degrees_north')
            .add_metadata('axis', 'Y'))

        self.variables.append(
            VariableSpecification('longitude', numpy.float32)
            .add_metadata('standard_name', 'longitude')
            .add_metadata('long_name', 'Longitude (deg)')
            .add_metadata('units', 'degrees_east')
            .add_metadata('axis', 'X'))

        self.variables.append(
            VariableSpecification('tas', numpy.int16)
            .add_metadata('standard_name', 'air_temperature')
            .add_metadata('long_name', 'Average daily surface air temperature')
            .add_metadata('units', 'K')
            .add_metadata('cell_methods', 'time: mean'))

        self.variables.append(
            VariableSpecification('tasmin', numpy.int16)
            .add_metadata('standard_name', 'air_temperature')
            .add_metadata('long_name', 'Minimum daily surface air temperature')
            .add_metadata('units', 'K')
            .add_metadata('cell_methods', 'time: minimum'))

        self.variables.append(
            VariableSpecification('tasmax', numpy.int16)
            .add_metadata('standard_name', 'air_temperature')
            .add_metadata('long_name', 'Maximum daily surface air temperature')
            .add_metadata('units', 'K')
            .add_metadata('cell_methods', 'time: maximum'))
