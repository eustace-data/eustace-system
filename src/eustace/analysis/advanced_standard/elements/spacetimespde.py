"""Base class for elements that use SPDE precision estimates in space and time."""

import numpy
import scipy.sparse

from element import Element
from element import Hyperparameters
from element import ElementDesign
from element import ElementPrior
from element import SPARSEFORMAT

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.advanced_standard.stats.spde.lattice import LatticeSPDE, WendlandC4Basis

class SpaceTimeSPDEHyperparameters(Hyperparameters):

    def __init__(self, space_log_sigma, space_log_rho, time_log_rho):
        """Initialise with specified parameters."""

        self.space_log_sigma = space_log_sigma
        self.space_log_rho = space_log_rho
        self.time_log_rho = time_log_rho
       
    def get_array(self):
        """Retrieve as 1D NumPy array."""
        
        valueslist = [ self.space_log_sigma, self.space_log_rho, self.time_log_rho ]
        return numpy.array(valueslist)
        
    def set_array(self, values):
        """Set from NumPy array."""

        self.space_log_sigma = values[0]
        self.space_log_rho = values[1]
        self.time_log_rho = values[2]
                
class SpaceTimeSPDEElement(Element):
    """Base class for SPDE element.  Inherit from this together with ElementLinear or ElementNonlinear."""

    def __init__(self, n_triangulation_divisions, alpha, starttime, endtime, n_nodes, overlap_factor, H, wrap_dimensions=None):
        """Initialise space and time SPDE definitions."""
        
        super(SpaceTimeSPDEElement, self).__init__()
        self.spatial_model = SphereMeshSPDE(level=n_triangulation_divisions, sparse_format=SPARSEFORMAT)
        self.alpha = alpha
        self.temporal_model = LatticeSPDE.construct(
            dimension_specification = [(starttime, endtime, n_nodes)],
            basis_function = WendlandC4Basis(),
            overlap_factor = overlap_factor,
            wrap_dimensions = wrap_dimensions)
        self.H = H

class SpaceTimeSPDEDesign(ElementDesign):
    """Base class for SPDE design. Inherit from this together with ElementLinearDesign or ElementNonlinearDesign as required."""

    def __init__(self, observationstructure, spatial_model, alpha, temporal_model, H):
        """Initialise space and time SPDE designs."""

        super(SpaceTimeSPDEDesign, self).__init__()

        self.spatial_model = spatial_model
        self.temporal_model = temporal_model

        self.alpha = alpha
        self.H = H

        self.spatial_locations = observationstructure.location_polar_coordinates()
        self.temporal_locations = numpy.array([ [ observationstructure.time_index() ] ])

        self.space_A = spatial_model.build_A(self.spatial_locations)
        self.time_A = temporal_model.build_A(self.temporal_locations)

class SpaceTimeSPDEPrior(ElementPrior):

    def __init__(self, hyperparameters, spatial_model, alpha, temporal_model, H):
        
        super(SpaceTimeSPDEPrior, self).__init__()

        self.hyperparameters = hyperparameters
        self.spatial_model = spatial_model
        self.temporal_model = temporal_model
        self.alpha = alpha
        self.H = H
