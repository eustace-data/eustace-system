"""Store output from file aggregation as NetCDF."""

__version__ = "$Revision: 552 $"
__author__ = "Joel R. Mitchelson"

import time
import numpy
from eustace.outputformats.filebuilder import FileBuilder
from eustace.outputformats import definitions
from eustace.outputformats.outputvariable import OutputVariable
from file_aggregate_grid import FLAG_INVALID

class FileBuilderAggregateField(FileBuilder):
    """Build a NetCDF file for aggregate field."""

    @staticmethod
    def save(filename, fields, source, complevel):
        """Save result as NetCDF file."""

        builder = FileBuilderAggregateField(filename, source, complevel)
        builder.setfields(fields)
        builder.save_and_close()

    def __init__(self, filename, source, complevel):
        """Create with specified attributes."""

        super(FileBuilderAggregateField, self).__init__()

        # make the file
        self.create(filename, title='EUSTACE aggregated data', institution='Met Office',
                    comment='', history='Created ' + time.strftime('%c'), source=source)

        # store compression level
        self.complevel = complevel

    def setfields(self, fields):
        """Add entire dictionary of fields."""

        # global latitude axis
        self.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_LATITUDE,
            variable=definitions.LATITUDE,
            values=fields[definitions.LATITUDE.name])

        # global longitude axis
        self.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_LONGITUDE,
            variable=definitions.LONGITUDE,
            values=fields[definitions.LONGITUDE.name])

        # count of observations used
        self.add_variable(
            OutputVariable(name='ts_number_of_observations',
                           dtype=numpy.int16,
                           fill_value=None,
                           standard_name=None,
                           long_name='Number of observations of surface temperature',
                           units='1'),
            (definitions.DIMENSION_NAME_LATITUDE,
             definitions.DIMENSION_NAME_LONGITUDE),
            fields['ts_number_of_observations'],
            self.complevel)

        # count of observations available (including cloudy)
        self.add_variable(
            OutputVariable(name='total_number_of_observations',
                           dtype=numpy.int16,
                           fill_value=None,
                           standard_name=None,
                           long_name='Total number of satellite observations in grid box',
                           units='1'),
            (definitions.DIMENSION_NAME_LATITUDE,
             definitions.DIMENSION_NAME_LONGITUDE),
            fields['total_number_of_observations'],
            self.complevel)

        # fields
        self.add_temperature_field(name='tsmean', standard_name='surface_temperature',
                                   long_name='Mean surface temperature (K)', values=fields['tsmean'])
        self.add_temperature_field(name='tsmax', standard_name='surface_temperature',
                                   long_name='Maximum surface temperature (K)', values=fields['tsmax'])
        self.add_temperature_field(name='tsmin', standard_name='surface_temperature',
                                   long_name='Minimum surface temperature (K)', values=fields['tsmin'])
        self.add_variance_field(name='tsvariance', standard_name=None,
                                long_name='Surface temperature variance (K2)', values=fields['tsvariance'])
        self.add_uncertainty_field(name='tsmean_unc_ran', standard_name=None,
                                   long_name='Random uncertainty in mean surface temperature (K)', values=fields['tsmean_unc_ran'])
        self.add_uncertainty_field(name='tsmean_unc_loc_atm', standard_name=None,
                                   long_name='Locally correlated uncertainty (atm component) in mean surface temperature (K)', values=fields['tsmean_unc_loc_atm'])
        self.add_uncertainty_field(name='tsmean_unc_loc_sfc', standard_name=None,
                                   long_name='Locally correlated uncertainty (sfc component) in mean surface temperature (K)', values=fields['tsmean_unc_loc_sfc'])
        self.add_uncertainty_field(name='tsmean_unc_sys', standard_name=None,
                                   long_name='Systematic uncertainty in mean surface temperature (K)', values=fields['tsmean_unc_sys'])
        self.add_uncertainty_field(name='tsmean_unc_spl', standard_name=None,
                                   long_name='Sampling uncertainty in mean surface temperature (K)', values=fields['tsmean_unc_spl'])

    def add_field(self, name, standard_name, long_name, units, scale, offset, values):
        """Add field (must have global coverage with respect to latitude and logitude axes)."""

        # Allow a large number of arguments here
        # - could consider having a namespace in future
        # pylint: disable=too-many-arguments

        self.add_variable(
            OutputVariable(
                name=name,
                dtype=numpy.int16,
                fill_value=-32768,
                scale_factor=scale,
                add_offset=offset,
                long_name=long_name,
                standard_name=standard_name,
                units=units),
            (definitions.DIMENSION_NAME_LATITUDE,
             definitions.DIMENSION_NAME_LONGITUDE),
            numpy.ma.masked_array(values, numpy.equal(values, FLAG_INVALID)),
            complevel=self.complevel)

    def add_temperature_field(self, name, standard_name, long_name, values):
        """Add a field with units and precision set for temperature."""

        return self.add_field(name, standard_name, long_name, 'K', 0.005, 273.15, values)

    def add_uncertainty_field(self, name, standard_name, long_name, values):
        """Add a field with units and precision set for uncertainty."""

        self.add_field(name, standard_name, long_name,
                       'K', 0.001, 32.767, values)

    def add_variance_field(self, name, standard_name, long_name, values):
        """Add a field with units and precision set for variance (K2)."""

        self.add_variable(
            OutputVariable(
                name=name,
                dtype=numpy.int32,
                fill_value=-1,
                scale_factor=0.000005,
                add_offset=0.0,
                long_name=long_name,
                standard_name=standard_name,
                units='K2'),
            (definitions.DIMENSION_NAME_LATITUDE,
             definitions.DIMENSION_NAME_LONGITUDE),
            numpy.ma.masked_array(values, numpy.equal(values, FLAG_INVALID)),
            complevel=self.complevel)
