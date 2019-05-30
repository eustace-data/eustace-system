"""Project output onto grid

TODO Setup for state reading for both spacetime and spatial components

"""

import numpy
import scipy.sparse

from eustace.timeutils.epoch import epoch_plus_days

from eustace.analysis.advanced_standard.components.spatial import SpatialComponentSolutionStorage
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponentSolutionStorage
from eustace.analysis.advanced_standard.components.storage_files_batch import SpatialComponentSolutionStorageIndexed_Files
from eustace.analysis.advanced_standard.components.storage_files_batch import SpaceTimeComponentSolutionStorageBatched_Files
from eustace.analysis.advanced_standard.fileio.output_structure_rectilinear import OutputRectilinearGridStructure

TIME_INDEXED_STORAGE_CLASSES = [SpatialComponentSolutionStorage, SpatialComponentSolutionStorageIndexed_Files]
INDEX_FREE_STORAGE_CLASSES = [SpaceTimeComponentSolutionStorage, SpaceTimeComponentSolutionStorageBatched_Files]


class Projector(object):
    """ Compute linear functions of the state as output = A_projector * state_vector.
    
    where A_projector is a sparse matrix that is retained for application 
    multiple times.
    The projector design_matrix is defined as a product of a linear operation matrix 
    and a the design matrix of a (linear) model component:
    A_projector = A_operation * A_component
    ----
    1 - allocate component
    2 - setup output locations (or biases)
    3 - evaluate design matrix
    4 - save/load design matrix is required
    5 - load model state
    6 - apply the design matrix
    """
    
    def __init__(self, centre_latitudes, centre_longitudes, grid_resolution, time_index, cell_sampling, blocking):
        
        self.outputstructure = None # TODO Allow for an output struture with no spatial location - pure covariate.
        
        self.grid_resolution = grid_resolution
        self.centre_latitudes = centre_latitudes
        self.centre_longitudes = centre_longitudes
        self.time_index = time_index
        self.corresponding_datetime = epoch_plus_days(time_index)
        self.cell_sampling = cell_sampling
        self.blocking = blocking
        
        self.component = None
        self.state_solution = None
        self.state_samples  = None
        self.prior_samples  = None
        
        self.design_matrix = None
    
    def set_component(self, component, keep_design = False):
        """Assign a component to this Projector.
        
        Used to assign a component with the same model structure but different state,
        e.g. for a EUSTACE analysis local component for a different time step.
        Set keep_design to True to allow reuse of a projection for a different component.
        
        """
        print "setting component"
        self.component = component
        self.state_solution = self.read_state_solution(self.component.solutionstorage)
        self.state_samples  = self.read_state_sample(self.component.solutionstorage)
        self.prior_samples  = self.read_prior_sample(self.component.solutionstorage)  # Reading/writing not implemented in solutionstorage classes
        
        print self.component
        print self.state_solution
        print self.state_samples
        print self.prior_samples
        
        if not keep_design:
            self.design_matrix = None
    
    def update_time_index(self, time_index, keep_design = False):
        """Update time time information for the projector
        
        Set keep_design to True to allow reuse of a projection for a different time_index,
        e.g. for a EUSTACE analysis local component for a different time step.
        
        """
    
        self.time_index = time_index
        self.corresponding_datetime = epoch_plus_days(time_index)
    
        if not keep_design:
            self.design_matrix = None
            
        print "New time index:", self.time_index
    
    def read_state_solution(self, solutionstorage):
        """Retrieve a state vector from component storage"""
        
        if any([isinstance(solutionstorage, time_indexed_class) for time_indexed_class in TIME_INDEXED_STORAGE_CLASSES]):
            time_index = self.time_index            
            state_vector = solutionstorage.partial_state_read(time_index)
        elif any([isinstance(solutionstorage, index_free_class) for index_free_class in INDEX_FREE_STORAGE_CLASSES]):
            state_vector = solutionstorage.state_read()
            
        else:
            raise TypeError("Unrecognised component storage type")

        return state_vector
        
    def read_state_sample(self, solutionstorage):
        """Retrieve a state samples from component storage"""
        
        if any([isinstance(solutionstorage, time_indexed_class) for time_indexed_class in TIME_INDEXED_STORAGE_CLASSES]):
            time_index = self.time_index 
            state_vector = solutionstorage.partial_state_sample_read(time_index)            
        elif any([isinstance(solutionstorage, index_free_class) for index_free_class in INDEX_FREE_STORAGE_CLASSES]):
            state_vector = solutionstorage.state_sample_read()
            
        else:
            raise TypeError("Unrecognised component storage type")

        return state_vector
    
    def read_prior_sample(self, solutionstorage):
        """Retrieve a state prior from component storage"""
        
        if any([isinstance(solutionstorage, time_indexed_class) for time_indexed_class in TIME_INDEXED_STORAGE_CLASSES]):
            time_index = self.time_index 
            state_vector = solutionstorage.partial_state_prior_sample_read(time_index)
            
        elif any([isinstance(solutionstorage, index_free_class) for index_free_class in INDEX_FREE_STORAGE_CLASSES]):
            state_vector = solutionstorage.state_prior_sample_read()
            
        else:
            raise TypeError("Unrecognised component storage type")

        return state_vector
    
    
    def project_expected_value(self):
        """Apply design matrix to precomputed component solution"""

        # Get state estimate (might be None before first iteration)
        if self.state_solution is None:
            return None
        else:
            return self.design_matrix.dot(self.state_solution)
        
    def project_sample_values(self, sample_indices = Ellipsis):
        """Project solution samples through design matrix"""

        # Get state estimate (might be None before first iteration)
        if self.state_samples is None:
            return None

        if sample_indices is not None:
            current_state_estimate = self.state_samples[:,sample_indices]
        else:
            current_state_estimate = self.state_samples
        
        return self.design_matrix.dot(current_state_estimate)

    def project_sample_std(self, prior = False):
        """Project solution samples through design matrix and compute the sample standard deviation"""
        
        if prior:
            statesample = self.prior_samples
        else:
            statesample = self.state_samples
        
        if statesample is None:
            return None
              
        sampledepartures = statesample - self.state_solution
        
        n_outputs = self.design_matrix.shape[0]
        output = numpy.zeros((n_outputs,1))
        blocksize = int(n_outputs / self.blocking)
        
        n_blocks = int( numpy.ceil( float(n_outputs) / float(blocksize) ) )
        
        for block in range(n_blocks):
            
            block_start = block*blocksize
            block_end = numpy.min( ( (block+1)*blocksize, n_outputs ) )
            
            output[block_start:block_end,0] = numpy.std( self.design_matrix[block_start:block_end,:].dot(statesample), axis=1 )
            
        return output

    def project_sample_deviation(self, prior = False):
        """Project solution samples through design matrix and compute the sample standard deviation"""
        
        if prior:
            # Assumes zero mean prior
            statesample = self.prior_samples
            sampledepartures = statesample#**2
        else:
            statesample = self.state_samples
            sampledepartures = statesample#(statesample - self.state_solution)**2
            print statesample.shape
            print self.state_solution.shape
            
        if statesample is None:
            return None
        
        n_outputs = self.design_matrix.shape[0]
        output = numpy.zeros((n_outputs,1))
        blocksize = int(n_outputs / self.blocking)
        
        n_blocks = int( numpy.ceil( float(n_outputs) / float(blocksize) ) )
        
        for block in range(n_blocks):
            
            block_start = block*blocksize
            block_end = numpy.min( ( (block+1)*blocksize, n_outputs ) )
            
            output[block_start:block_end,0] = numpy.sqrt( numpy.sum( self.design_matrix[block_start:block_end,:].dot(sampledepartures)**2, axis=1 ) / ( sampledepartures.shape[1] ) )
            
        return output

    def evaluate_design_matrix(self):
        """Produce output design matrices with optional within grid cell sampling"""

        component = self.component

        centre_latitudes = self.centre_latitudes
        centre_longitudes = self.centre_longitudes
        time_index = self.time_index
        corresponding_datetime = self.corresponding_datetime

        latitude_resolution = numpy.float(self.grid_resolution[0])
        longitude_resolution = numpy.float(self.grid_resolution[1])
        
        latitude_delta  = latitude_resolution / numpy.float(self.cell_sampling[0])
        longitude_delta = longitude_resolution / numpy.float(self.cell_sampling[1])
        
        # Loop through within cell sampling to compute cell averaging design matrices
        design_matrix = None
        weight_normalisation_array = None
        for latitude_index in range(self.cell_sampling[0]):
            for longitude_index in range(self.cell_sampling[1]):
                point_latitudes = centre_latitudes - latitude_resolution / 2.0 + (0.5 + latitude_index) * latitude_delta
                point_longitudes = centre_longitudes - longitude_resolution / 2.0 + (0.5 + longitude_index) * longitude_delta
                
                projectionstructure = OutputRectilinearGridStructure(time_index, corresponding_datetime, point_latitudes, point_longitudes)
            
                block_model_matrix = component.storage.element_read().element_design(projectionstructure).design_matrix()
                block_weight_array = self.weight_array(projectionstructure.location_polar_coordinates()[:,0])
                
                if design_matrix is None:
                    design_matrix = scipy.sparse.diags( block_weight_array ).dot(block_model_matrix)
                    weight_normalisation_array = block_weight_array
                else:
                    design_matrix += scipy.sparse.diags( block_weight_array ).dot(block_model_matrix)
                    weight_normalisation_array += block_weight_array
                    
        self.design_matrix = scipy.sparse.diags( 1.0 / weight_normalisation_array ).dot( design_matrix )

    def weight_array(self, point_latitudes):
        """Un-normalised weighting for within gridcell sampling"""
        return numpy.cos( numpy.radians(point_latitudes) )
    
    def save_projection(self, filename):
        """Save a copy of the Projector to disk without saving component_solution or state variables"""
        
        import pickle
        
        output_projector = Projector(self.centre_latitudes, self.centre_longitudes,
                                     self.grid_resolution, self.time_index, self.cell_sampling, self.blocking)
        
        output_projector.design_matrix = self.design_matrix
        output_projector.component = self.component
        
        with open(filename, 'wb') as f:
            pickle.dump( output_projector, f )
    
    @staticmethod
    def load_projection(filename):
        
        import pickle
        
        """Load a save Projector from disk"""
        with open(filename, 'rb') as f:
            projector = pickle.load( f )
            
        return projector

