"""Output to rectilinear grid."""

import numpy

from eustace.analysis.advanced_standard.elements.element import ObservationStructure

from eustace.outputformats.definitions import GLOBAL_FIELD_OUTPUT_FLAGS

class OutputRectilinearGridStructure(ObservationStructure):

    def __init__(self, time_index, corresponding_datetime, latitudes, longitudes):

        self.time_index_number = time_index
        self.corresponding_datetime = corresponding_datetime
        self.latitudes = latitudes
        self.longitudes = longitudes

    def time_index(self):
        """
        A discretised time index at appropriate resolution for input.
        In EUSTACE fullstace system these correspond to day numbers since
        01/01/1850.
        """
        
        return self.time_index_number

    def time_datetime(self):
        """
        The absolute datetime of these observations
        (needed for seasonal model and factor analysis which 
        must compute a fractional year value).
        """

        return self.corresponding_datetime

    def number_of_observations(self):
        """Return total observations at this time index."""

        return self.latitudes.shape[0] * self.longitudes.shape[0]
        
    def location_polar_coordinates(self):
        """Array of polar coordinates of observation location on unit sphere."""

        lon, lat = numpy.meshgrid(self.longitudes, self.latitudes)
        return numpy.vstack((lat.ravel(), lon.ravel())).T

class RegriddingRectilinearGridStructure(ObservationStructure):

    def __init__(self, time_index, corresponding_datetime, latitudes, longitudes, Nlat, Nlon):

        self.time_index_number = time_index
        self.corresponding_datetime = corresponding_datetime
        self.latitudes = latitudes
        self.longitudes = longitudes
        self.Nlat = Nlat
        self.Nlon = Nlon
	
    def time_index(self):
        """
        A discretised time index at appropriate resolution for input.
        In EUSTACE fullstace system these correspond to day numbers since
        01/01/1850.
        """
        
        return self.time_index_number

    def time_datetime(self):
        """The absolute datetime of these observations
           (needed for seasonal model and factor analysis which 
            must compute a fractional year value)."""

        return self.corresponding_datetime

    def number_of_observations(self):
        """Return total observations at this time index."""

        return self.latitudes.shape[0] * self.longitudes.shape[0]
        
    def rearrange_latitudes(self):
        """Rearrange latitudes for boosting regridding computations"""
        
        latitudes = self.latitudes.reshape(-1, self.Nlat)
        latitudes = numpy.tile(latitudes, (len(self.longitudes)/self.Nlon,self.Nlon))
        reindexer = numpy.argsort(latitudes[:,0])
        
        return latitudes[reindexer]#, longitudes

    def rearrange_longitudes(self):
        """Rearrange longitudes for boosting regridding computations"""
        
        longitudes = self.longitudes.reshape(-1, self.Nlon)
        longitudes = numpy.tile(longitudes, (len(self.latitudes)/self.Nlat,self.Nlat))
        reindexer = numpy.argsort(longitudes[0])
        
        return longitudes[:, reindexer]

    def location_polar_coordinates(self):
        """Array of polar coordinates of observation location on unit sphere."""
        
        lat = self.rearrange_latitudes()
        lon = self.rearrange_longitudes()
        return numpy.vstack((lat.ravel(), lon.ravel())).T

class Regridder(object):
    """Construct grid cell area averages of temperature fields, by using an underlying output grid"""
    
    def __init__(self, outputstructure, grid_points_per_cell, blocking=100):
      
        self.time_index_number = outputstructure.time_index_number
        self.corresponding_datetime = outputstructure.corresponding_datetime
        self.latitudes = outputstructure.latitudes
        self.longitudes = outputstructure.longitudes
        self.Nlat = grid_points_per_cell[0]
        self.Nlon = grid_points_per_cell[1]
        self.blocking = blocking
        self.latitude_spacing = round(self.latitudes[1]-self.latitudes[0], 3)
        self.longitude_spacing = round(self.longitudes[1]-self.longitudes[0],3)
        self.number_of_observations = outputstructure.number_of_observations()
        self.field_flags = GLOBAL_FIELD_OUTPUT_FLAGS
        
    def compute_gridded_expected_value(self, field, component_solution):
        """Performs regridding on input data points, returns regridded values"""      
        
        if self.blocking > len(self.latitudes):
            self.blocking = len(self.latitudes)
        
        rest = len(self.latitudes)%self.blocking
        grid_cell_average = numpy.array([])
        
        index=0
        for index in range(len(self.latitudes)/self.blocking):
            blocked_latitudes = self.latitudes[index*self.blocking:(index+1)*self.blocking]
            grid_cell_average = numpy.concatenate((grid_cell_average, self.compute_blocked_gridded_expected_value(field, blocked_latitudes, self.longitudes, component_solution)))
            
        if rest:
            blocked_latitudes = self.latitudes[(index+1)*self.blocking:]
            grid_cell_average = numpy.concatenate((grid_cell_average, self.compute_blocked_gridded_expected_value(field, blocked_latitudes, self.longitudes, component_solution)))
        
        return grid_cell_average
        
    def compute_blocked_gridded_expected_value(self, field, blocked_latitudes, longitudes, component_solution):
        """Performs regridding on a block of input data points, returns regridded values"""

        intermediate_output_structure = self.create_sub_cell_grids([blocked_latitudes, longitudes])
        latitudes_for_weights = intermediate_output_structure.rearrange_latitudes()
        
        if field == self.field_flags[0]:
            projected_values = component_solution.solution_observation_expected_value(intermediate_output_structure).reshape([-1, self.Nlat*self.Nlon])
        elif field == self.field_flags[1]:
            projected_values = component_solution.solution_observation_expected_uncertainties(intermediate_output_structure).reshape([-1, self.Nlat*self.Nlon])
        elif field == self.field_flags[2]:
            projected_values = component_solution.solution_observation_prior_uncertainties(intermediate_output_structure).reshape([-1, self.Nlat*self.Nlon])
        else:
            message = 'Output field flag {} not available. Available flags are {}'.format(field, self.field_flags)
            raise ValueError(message)
        
        grid_averaged_values = (self.weighting_factors(latitudes_for_weights)*projected_values).sum(axis=1)

        return grid_averaged_values.ravel()

    def create_sub_cell_grids(self, points):
        """Create a finer grid of points for a block  of grid cells"""
        
        Nlat_range = numpy.array(range(1, self.Nlat+1))
        Nlon_range = numpy.array(range(1, self.Nlon+1))
        
        latitudes = numpy.array([latitude - self.latitude_spacing/2.+(2*Nlat_range-1)*self.latitude_spacing/(2*self.Nlat) for latitude in points[0]]).ravel()
        longitudes = numpy.array([longitude - self.longitude_spacing/2.+(2*Nlon_range-1)*self.longitude_spacing/(2*self.Nlon) for longitude in points[1]]).ravel()
        
        return RegriddingRectilinearGridStructure(self.time_index_number, self.corresponding_datetime, latitudes, longitudes, self.Nlat, self.Nlon) 

    def weighting_factors(self, points):
        """Return the weighting factor used for weighting observations within a block of grid cells"""
        
        return self.harmonic_factor(points)/self.normalization_factor(points)
        
    def normalization_factor(self, points):
        """Compute the surface area of a collection of grid cells by taking into account surface curvature"""
        Nlat_tilda = self.harmonic_factor(points).sum(axis=1, keepdims=True)
        return Nlat_tilda

    def harmonic_factor(self, points): 
        """Return the harmonic factor used for computing grid cell area averages for a a block of grid cells"""
        
        return numpy.cos(numpy.radians(points))
