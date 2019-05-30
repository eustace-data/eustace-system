"""Extend iris functionality to provide list of gridded fields on common axes."""

from iris.coords import DimCoord
from iris.cube import Cube
from iris.cube import CubeList
import numpy

INDEXTYPE = numpy.int32
"""Storage type used for indices into the grid."""

class GridFieldAxis(object):
    """Parameters describing one axis of a field."""

    def __init__(self, zeroth, step, count, var_name, circular, units=None, long_name=None, standard_name=None):
        """Create and store the given parameters."""
        self.zeroth = zeroth
        self.step = step
        self.count = count
        self.var_name = var_name
        self.circular = circular
        self.units = units
        self.long_name = long_name
        self.standard_name = standard_name

    def axis_coordinates(self):
        """Create iris DimCoord suitable for building a cube."""
        return DimCoord.from_regular(**self.__dict__)

    def compute_indices(self, values):
        """Indices along this axis at which the given axis values are located."""

        bin_indices = (((values - self.zeroth) / self.step) - 0.5).astype(INDEXTYPE)
        if self.circular:
            wrapped = (bin_indices == INDEXTYPE(self.count)).nonzero()
            bin_indices[wrapped] = 0
        if (bin_indices < 0).any():
            unexpected_indices = numpy.nonzero(bin_indices < 0)
            unexpected_values = values[unexpected_indices]
            message = 'Values out of range: {0} at indices {1} (total {2})'.format(unexpected_values, unexpected_indices, len(unexpected_values))
            raise ValueError(message)
        if (bin_indices >= self.count).any():
            unexpected_indices = numpy.nonzero(bin_indices >= self.count)
            unexpected_values = values[unexpected_indices]
            message = 'Values out of range: {0} at indices {1} (total {2})'.format(unexpected_values, unexpected_indices, len(unexpected_values))
            # message = 'Values out of range: ' + str(unexpected_values) + ' at indices ' + str(unexpected_indices)
            raise ValueError(message)
        return bin_indices

class GridFieldDescriptor:
    """Descriptor for inputs/outputs to retreive."""

    def __init__(self, var_name, dtype=None, standard_name=None, long_name=None, units=None, usemask=True):
        self.var_name = var_name
        self.dtype = dtype
        self.standard_name = standard_name
        self.long_name = long_name
        self.units = units
        self.usemask = usemask
        
class GridFieldList(CubeList):
    """Provide list of data fields over common coordinates by extending iris CubeList."""

    def __init__(self, axes):
        """Create with specified axes."""
        
        # Store
        self.axes = axes

        # And also compute iris axes info
        self.axes_coordinates = [ (axes[index].axis_coordinates(), index) for index in range(len(axes)) ]


    def get_field(self, var_name):
       """Retrieve field corresponding to specified variable name, if exists, None otherwise."""
       
       try:
          return next(field for field in self if field.var_name == var_name)
       except StopIteration:
          return None

    def create_field(self, descriptor, data):
       """Create new field with specified descriptor attributes and data."""

       field = Cube(data,
                    var_name=descriptor.var_name,
                    standard_name=descriptor.standard_name,
                    long_name=descriptor.long_name,
                    units=descriptor.units,
                    dim_coords_and_dims=self.axes_coordinates)
       self.append(field)
       return field

    def get_or_create_field(self, descriptor, initialvalue):
        """Retrieve a field in the output corresponding to given descriptor variable name,
           or if it doesn't exist create and return a new blank one using the descriptor attributes."""

        # retrieve existing field if exists
        field = self.get_field(descriptor.var_name)

        # otherwise build new blank one
        if field is None:
            fielddimensions = [ self.axes[index].count for index in range(len(self.axes)) ]
            initialdata = numpy.empty(fielddimensions, dtype=descriptor.dtype)
            initialdata.fill(initialvalue)
            if descriptor.usemask:
               initialmask = numpy.ones(fielddimensions, dtype=numpy.bool)
               initialcontents = numpy.ma.masked_array(data=initialdata, mask=initialmask)
            else:
               initialcontents = initialdata
            field = self.create_field(descriptor, initialcontents)

        return field
