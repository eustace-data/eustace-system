"""Grid axis and two-dimensional grid representations."""

__version__ = "$Revision: 472 $"
__author__ = "Joel R. Mitchelson"

import numpy


class GridAxis(object):
    """Represent one axis of a regular grid by its start, increment, and length (number of increments)."""

    def __init__(self, start, increment, length, wrap=False):
        self.start = start
        self.increment = increment
        self.length = length
        self.wrap = wrap

    def compute_indices(self, values):
        """Determine the integer indices along this axis at which the given axis values are located."""
        bin_indices = ((values - self.start) / self.increment).astype(numpy.int32)
        if self.wrap:
            wrapped = (bin_indices == numpy.int32(self.length)).nonzero()
            bin_indices[wrapped] = 0
        if (bin_indices < 0).any():
            raise ValueError
        if (bin_indices >= self.length).any():
            raise ValueError
        return bin_indices

    def get_points(self):
        """Points along axis (edges of cells)."""
        return (numpy.mgrid[0:self.length] * self.increment) + self.start

    def get_centres(self):
        """Centres of cells along axis (between points)."""
        return (numpy.mgrid[0:self.length] * self.increment) + (self.start + (numpy.float32(0.5) * self.increment))


class GridLatLon(object):
    """Represent a two-dimensional regular grid onto which data can be aggregated."""

    def __init__(self, axis_lat, axis_lon):
        self.axis_lat = axis_lat
        self.axis_lon = axis_lon
        self.dimensions = (axis_lat.length, axis_lon.length)

    def compute_indices(self, coord_lat, coord_lon):
        """Compute the integer indices of grid boxes at which the given coordinates are located."""

        bin_lat = self.axis_lat.compute_indices(coord_lat)
        bin_lon = self.axis_lon.compute_indices(coord_lon)
        bin_indices_flat = (bin_lat * self.axis_lon.length) + bin_lon
        return bin_indices_flat
