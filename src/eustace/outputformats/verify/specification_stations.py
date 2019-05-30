"""Specification for station data verification."""

from specification import VariableSpecification
from specification import FileSpecification
import numpy


class SpecificationStationTemperature(FileSpecification):
    """Define the specification against which station temperature files are checked."""

    def __init__(self):

        # Set attributes and empty variables and dimensions
        super(SpecificationStationTemperature, self).__init__()

        # Set variables

        self.variables.append(
            VariableSpecification('time', numpy.float32)
            .add_metadata('standard_name', 'time')
            .add_metadata('long_name', 'Time at zero longitude')
            .add_metadata('units', 'days since 1850-01-01T00:00:00Z')
            .add_metadata('calendar', 'gregorian'))

        self.variables.append(
            VariableSpecification('latitude', numpy.float32)
            .add_metadata('standard_name', 'latitude')
            .add_metadata('long_name', 'Latitude (deg)')
            .add_metadata('units', 'degrees_north'))

        self.variables.append(
            VariableSpecification('longitude', numpy.float32)
            .add_metadata('standard_name', 'longitude')
            .add_metadata('long_name', 'Longitude (deg)')
            .add_metadata('units', 'degrees_east'))

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

        self.variables.append(
            VariableSpecification('station_name', numpy.dtype('S1'))
            .add_metadata('long_name', 'Station name')
            .add_metadata('cf_role', 'timeseries_id'))

        self.variables.append(
            VariableSpecification('elevation', numpy.float32)
            .add_metadata('long_name', 'Height above the geoid (m)')
            .add_metadata('standard_name', 'surface_altitude')
            .add_metadata('units', 'm'))

        self.variables.append(
            VariableSpecification('tasmin_qc', numpy.uint8)
            .add_metadata('long_name', 'Quality control flags')
            .add_metadata('standard_name', 'air_temperature status_flag')
            .add_metadata('valid_range', None)
            .add_metadata('flag_masks', None)
            .add_metadata('flag_meanings', None))

        self.variables.append(
            VariableSpecification('tasmax_qc', numpy.uint8)
            .add_metadata('long_name', 'Quality control flags')
            .add_metadata('standard_name', 'air_temperature status_flag')
            .add_metadata('valid_range', None)
            .add_metadata('flag_masks', None)
            .add_metadata('flag_meanings', None))

        # UNLIMITED time dimension
        self.dimensions['time'] = None


class SpecificationStationStatus(FileSpecification):
    """Define the specification against which station status files are checked."""

    def __init__(self):

        # Set attributes and empty variables and dimensions
        super(SpecificationStationStatus, self).__init__()

        self.dimensions['detection_time'] = None
        self.dimensions['tasmin_break'] = None
        self.dimensions['tasmax_break'] = None
        self.dimensions['bounds'] = 2

        # Add variables

        self.variables.append(
            VariableSpecification('detection_time', numpy.float32)
            .add_metadata('standard_name', 'time')
            .add_metadata('long_name', 'Start time of period for break detection status report (days)')
            .add_metadata('units', 'days since 1850-01-01T00:00:00Z'))

        self.variables.append(
            VariableSpecification('tasmin_detection_qc', numpy.uint8)
            .add_metadata('long_name', 'Quality control flags for break detection')
            .add_metadata('standard_name', 'air_temperature status_flag')
            .add_metadata('valid_range', [0,2])
            .add_metadata('flag_values', [0,1,2])
            .add_metadata('flag_meanings', 'not_possible possible_but_unreliable reliable'))

        self.variables.append(
            VariableSpecification('tasmax_detection_qc', numpy.uint8)
            .add_metadata('long_name', 'Quality control flags for break detection')
            .add_metadata('standard_name', 'air_temperature status_flag')
            .add_metadata('valid_range', [0,2])
            .add_metadata('flag_values', [0,1,2])
            .add_metadata('flag_meanings', 'not_possible possible_but_unreliable reliable'))

        self.variables.append(
            VariableSpecification('tasmin_break_station', numpy.int32)
            .add_metadata('long_name', 'Index of station at which break was detected'))

        self.variables.append(
            VariableSpecification('tasmax_break_station', numpy.int32)
            .add_metadata('long_name', 'Index of station at which break was detected'))

        self.variables.append(
            VariableSpecification('tasmin_break_amplitude', numpy.float32)
            .add_metadata('long_name', 'Inhomogeneity in minimum surface air temperature')
            .add_metadata('units', 'K'))

        self.variables.append(
            VariableSpecification('tasmax_break_amplitude', numpy.float32)
            .add_metadata('long_name', 'Inhomogeneity in maximum surface air temperature')
            .add_metadata('units', 'K'))

        self.variables.append(
            VariableSpecification('tasmin_break_time_bounds', numpy.float32)
            .add_metadata('long_name', 'Relative time bounds of break evaluation period (days)')
            .add_metadata('units', 'days'))

        self.variables.append(
            VariableSpecification('tasmax_break_time_bounds', numpy.float32)
            .add_metadata('long_name', 'Relative time bounds of break evaluation period (days)')
            .add_metadata('units', 'days'))

        self.variables.append(
            VariableSpecification(
                'tasmin_break_time_affected_bounds', numpy.float32)
            .add_metadata('long_name', 'Relative time bounds of period affected by break (days)')
            .add_metadata('units', 'days since 1850-01-01T00:00:00Z'))

        self.variables.append(
            VariableSpecification(
                'tasmax_break_time_affected_bounds', numpy.float32)
            .add_metadata('long_name', 'Relative time bounds of period affected by break (days)')
            .add_metadata('units', 'days since 1850-01-01T00:00:00Z'))
