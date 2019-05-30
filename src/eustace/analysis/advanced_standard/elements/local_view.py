"""Local model."""

from element import SPARSEFORMAT
from local import LocalElement, LocalHyperparameters, LocalPrior
from eustace.analysis.advanced_standard.stats.spde.spherical_view import SphereMeshViewLocal, SphereMeshViewSuperTriangle, SphereMeshViewGlobal
from eustace.analysis.advanced_standard.elements.element import Hyperparameters
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
import numpy

class LocalSubRegion(LocalElement):    
    """
    Class for evaluating design and precision matrices for SPDE models on a triangulated sphere.    
    """
    
    def __init__(self, n_triangulation_divisions, neighbourhood_level, centre_index_at_level):
        """
        Setup subregion of triangulation of sphere for :class:LocalSubRegion.

        Args:
        
        * n_triangulation_divisions (int):
            Number of subdivisions of icosohedron for triangulation of sphere.
        
        * neighbourhood_level (int):
            Subdivision level at which a neighbourhood is to be defined by super
            triangles.
            
        * centre_index_at_level (int):
            Index to the super triangle at the centre of the neighbourhood.
        
        """
        
        #super(LocalSubRegion, self).__init__()
        self.spde = SphereMeshViewLocal(n_triangulation_divisions, neighbourhood_level, centre_index_at_level, sparse_format=SPARSEFORMAT)

class LocalSuperTriangle(LocalElement):    
    """
    Class for evaluating design and precision matrices for SPDE models on a triangulated sphere.    
    """
    
    def __init__(self, n_triangulation_divisions, super_triangle_level, super_triangle_index):
        """
        Setup subregion of triangulation of sphere for :class:LocalSubRegion.

        Args:
        
        * n_triangulation_divisions (int):
            Number of subdivisions of icosohedron for triangulation of sphere.
        
        * super_triangle_level (int):
            Subdivision level at which a the super_triangle_index exists.
            
        * super_triangle_index (int):
            Index to the super triangle at the specified level.
        
        """
        
        #super(LocalSubRegion, self).__init__()
        self.spde = SphereMeshViewSuperTriangle(n_triangulation_divisions, super_triangle_level, super_triangle_index, sparse_format=SPARSEFORMAT)

class NonStationaryLocal(LocalElement):    
    """
    Class for evaluating design and precision matrices for SPDE models on a triangulated sphere.    
    """
    
    def __init__(self, n_triangulation_divisions):
        """
        Setup global view of triangulation of sphere for :class:LocalSubRegion.
        
        Gives access to merging of local parameterisations.

        Args:
        
        * n_triangulation_divisions (int):
            Number of subdivisions of icosohedron for triangulation of sphere.
        
        * neighbourhood_level (int):
            Subdivision level at which a neighbourhood is to be defined by super
            triangles.
            
        * centre_index_at_level (int):
            Index to the super triangle at the centre of the neighbourhood.
        
        """
        
        #super(NonStationaryLocal, self).__init__()
        self.spde = SphereMeshViewGlobal(n_triangulation_divisions, sparse_format=SPARSEFORMAT)
        
    def element_prior(self, hyperparameters):
        """
        Return LocalPrior instance for local model, using given hyperparameters.
        """
    
        return NonStationaryLocalPrior(hyperparameters, self.spde)

class NonStationaryLocalPrior(LocalPrior):

    def prior_precision(self):
        """
        Builds a precision matrix for the non-stationary local model with expanded parameterisation.
        
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
               
        """
        
        return self.spde.build_Q_expanded(alpha=self.alpha, **self.hyperparameters.get_dictionary())
    
    def prior_precision_derivative(self, parameter_index):
        """
        Builds the derivative of the precision matrix for the non-stationary local model with respect to the precision matrix hyperparmaters.
        
        Kwargs:
        
        * parameter_index:
            An index to the precision matrix parameter that the derivative is to be computed with 
            respect to.
                    
        Returns:
            scipy.sparse.spmatrix of requested sparse matrix format.
                
        """
        
        return self.spde.build_dQdp_expanded(alpha=self.alpha, parameter_index=parameter_index, **self.hyperparameters.get_dictionary())


#class ExtendedLocalHyperparameters(LocalHyperparameters):
    #"""Hyperparameters for local element."""
    
    #@staticmethod
    #def init_from_npy_savefile(filename):
        #"""Load hyperparameter values from npz save file"""
        
        #parameter_dict = numpy.load(filename)
        
        #return ExtendedLocalHyperparameters(**parameter_dict)
        
    #def values_to_npy_savefile(self, filename):
        #"""Write Q_parameters to npz save file"""
        
        #numpy.savez(filename, self.__dict__)


class ExtendedCombinationHyperparameters(CombinationHyperparameters):
    """Hyperparameters for local element."""
    
    def values_from_npy_savefile(self,filename):
        """Load hyperparameter values from npz save file"""
        
        values = numpy.load(filename)
        
        self.set_array(values)
        
    def values_to_npy_savefile(self, filename):
        """Write Q_parameters to npy save file"""
        
        numpy.save(filename, self.get_array())

class ExpandedLocalHyperparameters(LocalHyperparameters):
    """Hyperparameters for local element."""
    
    def get_array(self):        
        """Retrieve as 1D NumPy array."""
        
        valueslist = [ self.log_sigma, self.log_rho ]
        return numpy.concatenate(valueslist)
    
    def set_array(self, values):
        """Set from NumPy array."""

        self.log_sigma = values[0:len(values)/2]
        self.log_rho = values[len(values)/2:]
    
    def values_from_npy_savefile(self, filename):
        """Load hyperparameter values from npy save file"""
        
        values = numpy.load(filename)
        
        self.set_array(values)
        
    def values_to_npy_savefile(self, filename):
        """Write Q_parameters to npy save file"""
        
        numpy.save(filename, self.get_array())

class GeneralHyperparameters(Hyperparameters):
    """Hyperparameters for local element with expanded parameterisation."""
    
    def __init__(self, **Q_parameters):
        """Initialise with specified parameters."""
        
        for key, value in Q_parameters.iteritems():
            setattr(self, key, value)
       
    def get_array(self):
        """Retrieve as 1D NumPy array."""

        valueslist = [ self.__dict__[key] for key in sorted(self.__dict__)]
        return numpy.array(valueslist)
        
    def set_array(self, values):
        """Set from NumPy array."""
        
        startindex = 0
        for key in sorted(self.__dict__):    
            count = len(self.__dict__[key])
            nextstartindex = startindex + count
            setattr(self, key, values[startindex:nextstartindex])
            startindex = nextstartindex
                
    def get_dictionary(self):
        """Get as dictionary of named values."""
        
        return self.__dict__
    
    @staticmethod
    def init_from_npy_savefile(filename):
        """Load Q_parameters values from npz save file"""
        
        parameter_dict = numpy.load(filename)
        
        return ExpandedHyperparameters(**parameter_dict)
        
    def values_to_npy_savefile(self, filename):
        """Write Q_parameters to npy save file"""
        
        numpy.save(filename, self.__dict__)        