"""
Geography-based covariate (e.g. altitude, distance from water, vegetation fraction, etc) element.
"""

import scipy
import numpy

from ..stats.spde.lattice import PiecewiseLinearBasis, LatticeMesh
from covariate import CovariatePrior
from element import ElementLinear, ElementLinearDesign, ElementPrior
from element import Hyperparameters
from element import SPARSEFORMAT
from netCDF4 import Dataset
            
class GeographyBasedCovariateFunction(object):
    """Wrapper for mapping observations locations into interpolated physical covariates values."""

    def __init__(self, latitude, longitude, covariate, rescale_factor):
        """Construct build physical covariate design matrix. Get covariate and coordinates values, map observations locations into lineraly interpolated covariate values"""

        self.latitude = latitude
        self.longitude = longitude
        self.covariate = covariate
        self.rescale_factor = rescale_factor

    def compute(self, locations):
        """Compute the altitude design matrix for a given set of locations"""

        resultvector = self.rescale_factor * self.bilinear_interpolation(self.covariate, self.latitude, self.longitude, locations)
        return resultvector.reshape(-1,1)

    @ staticmethod
    def bilinear_interpolation(covariate, latitude, longitude, locations):
        """Apply bilinear interpolation on a regular grid to compute covariate values for a set of locations"""
        
        # generate a 2D grid, with one extra point on each side, enlarge the covariate array to ensure: 
        # a) Periodic Boundary Conditions over longitude: f(L+1,ny) = f(1,ny), f(0,ny) = f(L,ny) ny != (0, l+1)
        # b) Hard Wall Boundary conditions over latitude: f(nx, l+1) = f(nx, 0) = alpha, where alpha is the mean of the values taken from the nx=0,l layer
            # This would allow to correctly compute bilinear interpolation when a locations occur close to the lattice boundaries
        #
        # node nomenclature: node = Nnynx
        #		     N = n(original), e(extra)
        #		     nx = 0,...,L+1
        #		     ny = 0,...,l+1
        #
        #  e00      e01       e02       ...       e0L      e0(L+1)
        #	-----------------------------------------
        #	|         |         |         |         |  
        #  e10	|   n11   |   n12   |   ...   |   n1L   |  e1(L+1)
        #	|         |         |         |         |
        #       -----------------------------------------
        #	|         |         |         |         |  
        #  ...  |   ...   |   ...   |   ...   |   ...   |  ...
        #	|         |         |         |         |        
        #       -----------------------------------------
        #	|         |         |         |         |  
        #  el0	|   nl1   |   nl2   |   ...   |   nlL   |  el(L+1)
        #	|         |         |         |         |        
        #       -----------------------------------------
        #   e(l+1)0 n(l+1)1   n(l+1)2    ...      n(l+1)L   e(l+1)(L+1)
        #
        # here nodes spacing is supposed to be uniform, but for safety we check it

        nodes_spacing = round((latitude.max()-latitude.min())/(latitude.shape[0]-1),2)
        nodes_spacing_check = round((longitude.max()-longitude.min())/(longitude.shape[1]-1),2)
        span = numpy.repeat(nodes_spacing,2)
        
        if nodes_spacing != nodes_spacing_check:
            raise ValueError('Nodes spacing has to be uniform along all the directions!')

        # extend covariate vector over longitude and then latitude
        east_extension, west_extension = covariate[:,-1].reshape(-1,1), covariate[:,0].reshape(-1,1)
        covariate = numpy.hstack((east_extension, covariate))
        covariate = numpy.hstack((covariate, west_extension))
        north_extension = numpy.repeat(numpy.mean(covariate[0,:]),covariate.shape[1])
        south_extension = numpy.repeat(numpy.mean(covariate[-1,:]),covariate.shape[1])

        covariate = numpy.vstack((north_extension, covariate))
        covariate = numpy.vstack((covariate, south_extension))

        dimensions = [(latitude.min()-nodes_spacing, latitude.max()+nodes_spacing, latitude.shape[0]+2),(longitude.min()-nodes_spacing, longitude.max()+nodes_spacing, longitude.shape[1]+2)]

        lattice = LatticeMesh(dimension_specification = dimensions, basis_function = PiecewiseLinearBasis(basis_span = span))

        # to correctly match sparse matrix indexization
        covariate = covariate.ravel(order = 'F')

        basis_function_matrix = lattice.build_A(locations=locations)

        return basis_function_matrix.dot(covariate)

class GeographyBasedElement(ElementLinear):
    """Geography-based model element."""

    NUMBER_OF_STATE_PARAMETERS = 1

    def __init__(self,  filename, latitude_label, longitude_label, covariate_label, rescale_factor):
        """Initialise"""
        self.filename = filename
        self.latitude_label = latitude_label
        self.longitude_label = longitude_label
        self.covariate_label = covariate_label
        self.rescale_factor = rescale_factor

        super(GeographyBasedElement, self).__init__()

    def load(self):
        """Load data from a specific covariate file"""

        covariate_file = Dataset(self.filename)
        self.latitude = covariate_file.variables[self.latitude_label][:]
        self.longitude = covariate_file.variables[self.longitude_label][:]
        self.covariate = covariate_file.variables[self.covariate_label][:]
        covariate_file.close()

    def element_design(self, observationstructure):
        """Build the SeasonalElementDesign instance associated with this element's parameters
           and the given observation structure."""

        return GeographyBasedElementDesign(GeographyBasedCovariateFunction(self.latitude, self.longitude, self.covariate, self.rescale_factor), observationstructure)

    def element_prior(self, hyperparameters):
        """Build prior instance associated with given hyperparameters."""

        return GeographyBasedPrior(hyperparameters, GeographyBasedElement.NUMBER_OF_STATE_PARAMETERS)

class GeographyBasedElementDesign(ElementLinearDesign):

    def __init__(self, function_class, observationstructure):
        """Construct using latitudes and longitudes from given observation structure"""
        
        self.locations = observationstructure.location_polar_coordinates()
        self.function_class = function_class

    def design_number_of_state_parameters(self):
        """Number of parameters expected in state vector."""

        return GeographyBasedElement.NUMBER_OF_STATE_PARAMETERS

    def design_matrix(self):
        """Compute design matrix that maps model parameters to observation space."""

        return scipy.sparse.csr_matrix(self.function_class.compute(self.locations))

class GeographyBasedPrior(CovariatePrior):
    """Prior for geography-based covariates"""

    def __init__(self, hyperparameters, number_of_state_parameters):
        super(GeographyBasedPrior, self).__init__(hyperparameters, number_of_state_parameters)
