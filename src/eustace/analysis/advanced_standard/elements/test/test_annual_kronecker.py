"""Tests for space-time kronecker product analysis."""

import unittest
import numpy
import scipy.sparse
from datetime import datetime

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.mesh.geometry import cartesian_to_polar2d
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT
from eustace.analysis.advanced_standard.elements.element import ObservationStructure

from eustace.analysis.advanced_standard.elements.spacetimespde import SpaceTimeSPDEHyperparameters

from eustace.analysis.advanced_standard.elements.kronecker_annual import AnnualKroneckerElement
from eustace.analysis.advanced_standard.elements.kronecker_annual import AnnualKroneckerDesign


from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.advanced_standard.stats.spde.lattice import LatticeSPDE, WendlandC4Basis

from eustace.timeutils.decimaltime import datetime_to_decimal_year
from eustace.timeutils import epoch
from eustace.analysis.mesh.geometry import cartesian_to_polar2d

def round_to_quantisation(value, quantisation):
    return (value+quantisation/2) // quantisation * quantisation

def floor_to_quantisation(value, quantisation):
    return value // quantisation * quantisation

def ceil_to_quantisation(value, quantisation):
    return numpy.ceil( value / quantisation ) * quantisation

class TestSpaceTimeKroneckerElement(unittest.TestCase):

    class SimulatedObservationStructure(ObservationStructure):
        
        def __init__(self, polar_locations, obs_datetime):
        
            self.locations = polar_locations
            self.observations = numpy.ones(len(polar_locations))
            self.datetime = obs_datetime
        
        def location_polar_coordinates(self):

            return self.locations

        def time_datetime(self):

            return self.datetime

        def time_index(self):
            return epoch.days_since_epoch(self.datetime)
            
        def number_of_observations(self):
            return len(self.observations)

    def test_init(self):

        # Desired minimum model period
        start_time = datetime_to_decimal_year( datetime(1850, 1, 1) )
        end_time   = datetime_to_decimal_year( datetime(1860, 1, 1) )
        
        # Model time resolution is quarterly
        node_spacing = 0.25
        
        # Model start/end shifted to nearst grid times encompasing the desred start and end times
        model_start = floor_to_quantisation(start_time, node_spacing)
        model_end = ceil_to_quantisation(end_time, node_spacing)
        
        # number of model nodes in period
        n_nodes = int( round( (end_time - start_time) / node_spacing + 1 ) )
        
        element = AnnualKroneckerElement(n_triangulation_divisions=1, alpha=2, starttime=model_start, endtime=model_end, n_nodes=n_nodes, overlap_factor=2.5, H=1.0, wrap_dimensions = None)
        
        numpy.testing.assert_almost_equal( numpy.linspace(start_time, end_time, n_nodes).reshape(-1,1), element.temporal_model.lattice.nodes )
        
    def test_wrap(self):

        # Model time resolution
        node_spacing = 1./12.0
        
        # Model start/end bound a decimal year
        model_start = 1850.0
        model_end = 1851.0
        
        # number of model nodes in as 12 nodes within year
        n_nodes = 12

        element = AnnualKroneckerElement(n_triangulation_divisions=1, alpha=2, starttime=model_start, endtime=model_end, n_nodes=n_nodes, overlap_factor=2.5, H=1.0, wrap_dimensions = [True,])
        
        # Check unit diagonal basis function values when point at beginning of years and a each node in the spatial triangulation
        node_locations_space = cartesian_to_polar2d( element.spatial_model.triangulation.points )
        node_locations_time = element.temporal_model.lattice.nodes
        
        obs = TestSpaceTimeKroneckerElement.SimulatedObservationStructure(node_locations_space, datetime(1849, 1, 1) )
        numpy.testing.assert_almost_equal( numpy.ones(obs.number_of_observations()),  element.element_design(obs).design_matrix().diagonal() )
        obs = TestSpaceTimeKroneckerElement.SimulatedObservationStructure(node_locations_space, datetime(1850, 1, 1) )
        numpy.testing.assert_almost_equal( numpy.ones(obs.number_of_observations()),  element.element_design(obs).design_matrix().diagonal() )
        obs = TestSpaceTimeKroneckerElement.SimulatedObservationStructure(node_locations_space, datetime(1851, 1, 1) )
        numpy.testing.assert_almost_equal( numpy.ones(obs.number_of_observations()),  element.element_design(obs).design_matrix().diagonal() )
        