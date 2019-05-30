"""Build NetCDF file using python netCDF4 library directly."""

__version__ = "$Revision: 841 $"
__author__ = "Joel R. Mitchelson"

import os
import numpy
from netCDF4 import Dataset
from eustace.outputformats.definitions import FILE_FORMAT
from eustace.outputformats.definitions import FILE_CONVENTIONS

class FileBuilder(object):
    """Build a NetCDF file."""

    def __init__(self):
        self.nc = None # pylint: disable=invalid-name

    def create(self, pathname, title, institution, comment, history, source):
        """Create with appropriate attributes."""

        self.nc = Dataset(pathname, 'w', FILE_FORMAT)
        self.nc.Conventions = FILE_CONVENTIONS
        self.nc.title = title
        self.nc.institution = institution
        self.nc.history = history
        self.nc.source = source
        self.nc.comment = comment

    def add_dimension(self, dimensionname, size):
        """Add the requested dimension.  If size is None it will be UNLIMITED."""

        # create the dimension
        self.nc.createDimension(dimensionname, size)

    def add_variable(self, variable, dimension_names, values, complevel=0):
        """Add the requested variable, dimensions it depends upon must already exist."""

        # zlib info
        zlib = (complevel > 0)

        # make the NetCDF item
        ncvar = self.nc.createVariable(variable.name, variable.dtype, dimension_names,
                                       fill_value=variable.fill_value, zlib=zlib, complevel=complevel)

        # set auto-mask
        ncvar.set_auto_maskandscale(True)

        # copy all metadata
        for (key, metadata) in variable.get_metadata().iteritems():
            ncvar.__setattr__(key, metadata)

        # assign data
        if len(values.shape) == 1:
            ncvar[0:values.shape[0]] = values
        else:
            ncvar[:] = values

    def add_dimension_and_variable(self, dimensionname, variable, values, unlimited=False):
        """Add the requested dimension and a corresponding variable with axis values."""

        # create requested dimension using values to indicate size (or
        # unlimited if unlimited is requested)
        self.add_dimension(dimensionname, None if unlimited else values.shape[0])

        # populate with values
        self.add_variable(variable, (dimensionname, ), values)

    def save_and_close(self):
        """Close the file and store results."""
        self.nc.close()
