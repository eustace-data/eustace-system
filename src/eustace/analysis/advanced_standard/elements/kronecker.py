"""Space-time kronecker element.  This is just like a 3D grid using two coordinates of space and one of time."""

from element import ElementLinear
from element import ElementLinearDesign
from spacetimespde import SpaceTimeSPDEElement
from spacetimespde import SpaceTimeSPDEDesign
from spacetimespde import SpaceTimeSPDEPrior
from element import SPARSEFORMAT
import scipy.sparse
import numpy

TIME_ALPHA = 1

class SpaceTimeKroneckerElement(ElementLinear, SpaceTimeSPDEElement):
        
    def element_design(self, observationstructure):

        return SpaceTimeKroneckerDesign(observationstructure, self.spatial_model, self.alpha, self.temporal_model, self.H)

    def element_prior(self, hyperparameters):

        return SpaceTimeKroneckerPrior(hyperparameters, self.spatial_model, self.alpha, self.temporal_model, self.H)

class SpaceTimeKroneckerDesign(ElementLinearDesign, SpaceTimeSPDEDesign):

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

class SpaceTimeKroneckerPrior(SpaceTimeSPDEPrior):

    UNIT_SCALE_PARAMETER = numpy.log(1.0)

    def prior_number_of_state_parameters(self):

        return self.spatial_model.n_latent_variables() * self.temporal_model.n_latent_variables()
    
    def prior_precision(self):

        Q_spatial = self.spatial_model.build_Q_stationary(self.hyperparameters.space_log_sigma, self.hyperparameters.space_log_rho, self.alpha)
        Q_temporal = self.temporal_model.build_Q_stationary(SpaceTimeKroneckerPrior.UNIT_SCALE_PARAMETER, self.hyperparameters.time_log_rho, TIME_ALPHA, self.H, SPARSEFORMAT)
        return scipy.sparse.kron(Q_temporal, Q_spatial, format=SPARSEFORMAT)

    def prior_precision_derivative(self, parameter_index):

        if parameter_index < 2:

            dQ_spatial = self.spatial_model.build_dQdp_stationary(self.hyperparameters.space_log_sigma, self.hyperparameters.space_log_rho, self.alpha, parameter_index)
            dQ_temporal = self.temporal_model.build_Q_stationary(SpaceTimeKroneckerPrior.UNIT_SCALE_PARAMETER, self.hyperparameters.time_log_rho, TIME_ALPHA, self.H, SPARSEFORMAT)
            
        else:

            dQ_spatial = self.spatial_model.build_Q_stationary(self.hyperparameters.space_log_sigma, self.hyperparameters.space_log_rho, self.alpha)
            dQ_temporal = self.temporal_model.build_dQdp_stationary(SpaceTimeKroneckerPrior.UNIT_SCALE_PARAMETER, self.hyperparameters.time_log_rho, TIME_ALPHA, self.H, 1, SPARSEFORMAT)

        return scipy.sparse.kron(dQ_temporal, dQ_spatial, format=SPARSEFORMAT)


