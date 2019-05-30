"""Elements for factor analysis."""

import numpy
import scipy.sparse

from element import ElementNonlinear
from element import ElementNonlinearDesign
from element import ElementPrior
from element import SPARSEFORMAT

from spacetimespde import SpaceTimeSPDEHyperparameters
from spacetimespde import SpaceTimeSPDEElement
from spacetimespde import SpaceTimeSPDEDesign
from spacetimespde import SpaceTimeSPDEPrior

from eustace.analysis.advanced_standard.stats.spde.spherical import SphereMeshSPDE
from eustace.analysis.advanced_standard.stats.spde.lattice import LatticeSPDE, WendlandC4Basis

class SpaceTimeFactorElement(ElementNonlinear, SpaceTimeSPDEElement):

    def element_design(self, observationstructure):
        """Build the ElementDesign instance associated with this element's parameters
           and the given observation structure."""

        return SpaceTimeFactorDesign(observationstructure, self.spatial_model, self.alpha, self.temporal_model, self.H)

    def element_prior(self, hyperparameters):

        return SpaceTimeFactorPrior(hyperparameters, self.spatial_model, self.alpha, self.temporal_model, self.H)

class SpaceTimeFactorDesign(ElementNonlinearDesign, SpaceTimeSPDEDesign):

    def design_initial_state(self):

        return numpy.hstack([numpy.zeros((self.spatial_model.n_latent_variables(),)) , numpy.ones((self.temporal_model.n_latent_variables(),))])

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        return self.spatial_model.n_latent_variables() + self.temporal_model.n_latent_variables()
        
    def design_jacobian(self, currentstate):

        # Compute vectors
        v_space = self.space_A.dot(currentstate[:self.space_A.shape[1]])
        v_time = self.time_A.dot(currentstate[self.space_A.shape[1]:])

        # Express as 1 x n matrix (assumes currentstate was a 1D vector object)
        v_space = numpy.expand_dims(v_space, 1)

        # This is a kronecker product but time is always just one element here
        # so write it out explicitly
        jacobian_space = v_time[0] * self.space_A

        # Kronecker product for time part of derivative
        jacobian_time = scipy.sparse.kron(self.time_A, v_space, format=SPARSEFORMAT)

        # Stack 'em up
        return scipy.sparse.hstack([ jacobian_space, jacobian_time ], format=SPARSEFORMAT)

    def design_function(self, currentstate):

        v_space = self.space_A.dot(currentstate[:self.space_A.shape[1]])
        v_time = self.time_A.dot(currentstate[self.space_A.shape[1]:])
        return v_space * v_time

class SpaceTimeFactorPrior(SpaceTimeSPDEPrior):

    def prior_number_of_state_parameters(self):

        return self.spatial_model.n_latent_variables() + self.temporal_model.n_latent_variables()
    
    def prior_precision(self):

        return scipy.sparse.block_diag(
            [ self.spatial_model.build_Q_stationary(self.hyperparameters.space_log_sigma, self.hyperparameters.space_log_rho, self.alpha),
              self.temporal_model.build_Q_stationary(0.0, self.hyperparameters.time_log_rho, self.alpha, self.H, SPARSEFORMAT) ],
            format=SPARSEFORMAT)

    def prior_precision_derivative(self, parameter_index):

        matrix_list = [ ]

        if parameter_index < 2:

            matrix_list = [ 
                self.spatial_model.build_dQdp_stationary(self.hyperparameters.space_log_sigma, self.hyperparameters.space_log_rho, self.alpha, parameter_index),
                scipy.sparse.csc_matrix((self.temporal_model.n_latent_variables(), self.temporal_model.n_latent_variables())) ]
            
        else:

            matrix_list = [
              scipy.sparse.csc_matrix((self.spatial_model.n_latent_variables(), self.spatial_model.n_latent_variables())),
              self.temporal_model.build_dQdp_stationary(0.0, self.hyperparameters.time_log_rho, self.alpha, self.H, 1, SPARSEFORMAT) ]

        return scipy.sparse.block_diag(matrix_list, format=SPARSEFORMAT)
