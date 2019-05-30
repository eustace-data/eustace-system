"""Binning onto rectilinear grid using input coordinates specified by auxiliary variables."""

from gridfieldlist import GridFieldList
import numpy

INDEXTYPE = numpy.int32
"""Storage type used for indices into the grid."""

class GridBinMap(object):
    """Represent mapping between a given input cube object and corresponding grid bins, 
       so that auxiliary variables can be aggregated using the same indices."""

    def __init__(self, outputgrid, inputcube, mapfrom, mapto):
        """Create from given parameters."""
        
        self.outputgrid = outputgrid
        self.inputcube = inputcube
        self.mapfrom = mapfrom
        self.mapto = mapto

    
    def get_input_array(self, descriptor):
        """Retrieve input using the descriptor variable name (descriptor.var_name)."""

        # Get coord (even if auxiliary)
        if self.inputcube.var_name == descriptor.var_name:
            result = self.inputcube.data
        else:
            try:
               auxcoord = next(coord for coord in self.inputcube.coords() if coord.var_name == descriptor.var_name)
               result = auxcoord.points
            except StopIteration:
              raise ValueError('Coord not found: ' + descriptor.var_name)
        
        # Detach mask if exists
        if hasattr(result, 'mask'):
           return result.data
        else:
           return result

class GridBins(GridFieldList):
    """A grid to aggregate onto."""

    def __init__(self, axes):
        """Create with specified 2D axes."""
        
        if (len(axes) != 2) and (len(axes) != 3):
            raise ValueError('Must have either two or three axes.')

        if (len(axes) == 3) and (axes[0].count != 1):
            raise ValueError('When using three axes the first should have one data point only')
        
        super(GridBins, self).__init__(axes)


    def compute_input_map(self, inputcube, alternative_mask=None, input_coordinate_names=None):
        """Compute mapping between a given input cube object and corresponding grid bins, 
           so that auxiliary variables can be aggregated using the same indices.
           Uses the mask associated with the iris cube unless alternative_mask is specified.
           Assumes input coordinate axis names match output grid unless input_coordinate_names specified."""

        # Use alternative mask if specified, otherwise input mask
        mask = alternative_mask if (alternative_mask is not None) else inputcube.data.mask

        # Indices into flat input arrays which contain unmasked data
        mapfrom = numpy.nonzero(numpy.ravel(mask) == False)[0].astype(INDEXTYPE)

        # Build bin indices along each coordinate axis
        axis_indices = [ ]
        for coordinate_index, axis in enumerate(self.axes):
            var_name = input_coordinate_names[coordinate_index] if (input_coordinate_names is not None) else axis.var_name
            input_values_all = inputcube.coord(var_name).points.ravel()
            input_values_unmasked_only = input_values_all[mapfrom]
            indices = axis.compute_indices(input_values_unmasked_only)
            axis_indices.append(indices)

        # Multiply up to make global index
        mapto = (axis_indices[-2] * INDEXTYPE(self.axes[-1].count)) + axis_indices[-1]

        # Create result
        return GridBinMap(self, inputcube, mapfrom, mapto)

    def get_field(self, var_name):
       """Retrieve field corresponding to specified variable name, if exists, None otherwise."""
       
       try:
          return next(field for field in self if field.var_name == var_name)
       except StopIteration:
          return None
