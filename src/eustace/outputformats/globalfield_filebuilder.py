"""Functions to build EUSTACE compliant NetCDF files."""

__version__ = "$Revision: 1300 $"
__author__ = "Joel R. Mitchelson"

import os
import errno
import numpy
from netCDF4 import Dataset
from eustace.outputformats import definitions
from eustace.outputformats.outputvariable import OutputVariable
from eustace.timeutils.epoch import epoch_plus_days
from eumopps.timeutils import datetime_numeric
from filebuilder import FileBuilder
from ensuredirectory import ensuredirectory

class DatasetAttributesGlobalField(object):
    """Hold dataset metadata and build filenames."""

    def __init__(self, dataset, version, mainvariable, institution, comment, history, source):
        """Construct from given parameters."""

        self.dataset = dataset
        self.version = version
        self.mainvariable = mainvariable
        self.institution = institution
        self.comment = comment
        self.history = history
        self.source = source

    def build_pathname_patterns(self):
        """List of pathname components with placeholders for any date fields."""
        
        # Replace YYYY with %Y in subdirectory pattern
        subdirectory = definitions.GLOBAL_FIELD_SUBDIRECTORY_PATTERN.format(YYYY='%Y')

        # And build filename with similar placeholders
        filename = definitions.GLOBAL_FIELD_FILENAME_PATTERN.format(
            variable=self.mainvariable,
            collection=self.dataset.lower(),
            framework=definitions.FRAMEWORK_NAME,
            realisation='0',
            YYYY='%Y', mm='%m', dd='%d')

        return [ subdirectory, filename ]

    def build_pathname(self, basepath, daynumber):
        """ Build whole pathname from basepath and daynumber. """
        
        # Get pattern
        patterns = self.build_pathname_patterns()

        # Convert daynumber to datetime object
        fieldtime = epoch_plus_days(daynumber)

        # initialise pathname to basepath
        pathname = basepath

        # Fill date fields and append
        for pattern in patterns:
            pathname = os.path.join(pathname, datetime_numeric.build_from_pattern(pattern, fieldtime))
        
        return pathname

class FileBuilderGlobalField(FileBuilder):
    """Build a NetCDF file with EUSTACE specification for global fields."""

    def __init__(self, pathname, daynumber, dataset, version, mainvariable, institution, comment, history, source):
        """Create with specified attributes."""

        super(FileBuilderGlobalField, self).__init__()

        # create subdirectory if it doesn't exist
        ensuredirectory(pathname)

        # title to put inside file
        title = definitions.TITLE_PREFIX + dataset

        # make the file
        self.create(pathname, title, institution, comment, history, source)

        # set time in days since the EUSTACE epoch [and set UNLIMITED]
        self.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_TIME,
            variable=definitions.TIME,
            values=numpy.array([daynumber], numpy.float32),
            unlimited=True)

        # also set time bounds array (indicating midnight to midnight)
        self.add_dimension(definitions.DIMENSION_NAME_BOUNDS, 2)

        # set time bounds array
        self.add_variable(definitions.TIME_BOUNDS, (definitions.DIMENSION_NAME_TIME, definitions.DIMENSION_NAME_BOUNDS), numpy.array(
            [[daynumber, daynumber + 1.0]], numpy.float32))

        # global latitude axis
        self.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_LATITUDE,
            variable=definitions.LATITUDE,
            values=numpy.linspace(-90.+definitions.GLOBAL_FIELD_RESOLUTION/2., 90.-definitions.GLOBAL_FIELD_RESOLUTION/2., num=definitions.GLOBAL_FIELD_SHAPE[1]))#numpy.arange(-90.0, 89.7501, 0.25, numpy.float32))

        # global longitude axis
        self.add_dimension_and_variable(
            dimensionname=definitions.DIMENSION_NAME_LONGITUDE,
            variable=definitions.LONGITUDE,
            values=numpy.linspace(-180.+definitions.GLOBAL_FIELD_RESOLUTION/2., 180.-definitions.GLOBAL_FIELD_RESOLUTION/2., num=definitions.GLOBAL_FIELD_SHAPE[2]))#numpy.arange(-180.0, 179.7501, 0.25, numpy.float32))

        # local offset (for local solar time) should correspond to global longitude axis
        self.add_variable(definitions.TIME_OFFSET, (definitions.DIMENSION_NAME_LONGITUDE,),
                          numpy.linspace(-180.+definitions.GLOBAL_FIELD_RESOLUTION/2., 180.-definitions.GLOBAL_FIELD_RESOLUTION/2., num=definitions.GLOBAL_FIELD_SHAPE[2])/360. )

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
