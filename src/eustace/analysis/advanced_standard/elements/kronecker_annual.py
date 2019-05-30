"""Space-time kronecker element.  The mesh is just like a 3D grid using two coordinates of space and one of time.

The prior model is a seperable model with independent modeling of covariance in space and in time.

"""

from element import ElementLinearDesign
from element import SPARSEFORMAT
from eustace.timeutils.decimaltime import datetime_to_decimal_year
import scipy.sparse
import numpy

from kronecker import SpaceTimeKroneckerElement

class AnnualKroneckerElement(SpaceTimeKroneckerElement):
    """
    
    Kronecker element with time specified in decimal years.
    
    """
    
    
    def element_design(self, observationstructure):

        return AnnualKroneckerDesign(observationstructure, self.spatial_model, self.alpha, self.temporal_model, self.H)        
        
        
class AnnualKroneckerDesign(ElementLinearDesign):
    """
    
    Replicates functionality of SpaceTimeKroneckerDesign but inherits from AnnualSpaceTimeSPDEDesign.
    
    """
    
    def __init__(self, observationstructure, spatial_model, alpha, temporal_model, H):
        """Initialise space and time SPDE designs."""

        super(AnnualKroneckerDesign, self).__init__()

        self.spatial_model = spatial_model
        self.temporal_model = temporal_model

        self.alpha = alpha
        self.H = H

        self.spatial_locations = observationstructure.location_polar_coordinates()
        
        decimal_year = datetime_to_decimal_year(observationstructure.time_datetime())
        self.temporal_locations = numpy.array([ [ decimal_year ] ])

        self.space_A = spatial_model.build_A(self.spatial_locations)
        self.time_A = temporal_model.build_A(self.temporal_locations)
    
    def design_initial_state(self):

        return numpy.zeros((self.design_number_of_state_parameters(),))

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector (this is space x time)."""

        return self.spatial_model.n_latent_variables() * self.temporal_model.n_latent_variables()
    
    def design_matrix(self):
        
        return scipy.sparse.kron(self.time_A, self.space_A, format=SPARSEFORMAT)
    
    def design_jacobian(self, currentstate):

        return self.design_matrix()

    def design_function(self, currentstate, rectified=False):

        if rectified:
            return numpy.abs(self.design_jacobian(currentstate)).dot(currentstate)
        else:
            return self.design_jacobian(currentstate).dot(currentstate)
            
    