"""Connection to ObservationSource that implements a ObservationStructure interface."""

import os.path
import numpy
import scipy.sparse

from eumopps.catalogue.namepatterns import NameCandidate
from eumopps.catalogue.namepatterns import nested_search_patterns
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.fileio.observation_source_surface_effects import insitu_land_covariate_effect
from eustace.analysis.advanced_standard.fileio.observation_source_surface_effects import global_satellite_effect, hemispheric_satellite_effect
from eustace.outputformats.definitions import BIAS_SOURCES_PATTERN_DICTIONARY

from eustace.timeutils.epoch import epoch_plus_days

class ObservationStructureSourceConnector(ObservationStructure):

    def __init__(self, observationsource, observable, corresponding_datetime):
        
        self.location_lookup = observationsource.observation_location_lookup()
        self.observations = observationsource.observations(observable)
        self.valid_indices = numpy.nonzero(~self.observations.mask)
        self.corresponding_datetime = corresponding_datetime
        self.observable_filenames = [os.path.basename(observationsource.filespecs[key].filename) for key in observationsource.filespecs.keys() if observationsource.filespecs[key].filename is not None]


    def time_index(self):
        """
        A discretised time index at appropriate resolution for input.

        In EUSTACE fullstace system these correspond to day numbers since
        01/01/1850.

        In HadCRUT4 example these correspond to month numbers.

        """
        
        # Removed: observation time is not defined to be relative to analysis start time
        return self.observations.time
        
        #return epoch_plus_days(self.corresponding_datetime)

    def time_datetime(self):
        """The absolute datetime of these observations
        (needed for seasonal model and factor analysis which 
        must compute a fractional year value).
        """

        return self.corresponding_datetime

    def number_of_observations(self):
        """Return total observations at this time index."""

        return len(self.valid_indices[0])
        
    def location_polar_coordinates(self):
        """Array of polar coordinates of observation location on unit sphere."""
        locations = self.location_lookup[:, self.observations.location[self.valid_indices] ].T
        
        wrapped_longitudes = numpy.nonzero(locations[:,1] >= numpy.float(180.0))
        locations[wrapped_longitudes,1] -= 360.0
        
        return locations
        
        #return self.location_lookup[:, self.observations.location[self.valid_indices] ].T
                

    def observation_vector(self):
        """Array of observations."""

        return self.observations.measurement[self.valid_indices]

    def observation_precision(self):
        """Observation precision matrix (sparse)."""

        return scipy.sparse.diags(1.0 / (self.observations.uncorrelatederror[self.valid_indices]**2), format='csc') 

    def covariate_effect(self, groupname, **kwargs):
        """
        Retrieve details of bias that affects these observations.
        NumPy matrix of integers with 2 columns.  
        First column gives indices of observations affected.
        Second column is index of corresponding bias parameter affecting it.
        None if no effect defined.
        """
  
        if (groupname in BIAS_SOURCES_PATTERN_DICTIONARY.keys()):
            if self.check_all_candidates(BIAS_SOURCES_PATTERN_DICTIONARY[groupname]):
                if groupname == 'insitu_land' and kwargs:
                    return insitu_land_covariate_effect(self.time_index(), self.observations.location[self.valid_indices], kwargs['breakpoints'])
                elif (groupname == 'surfaceairmodel_land_global') or (groupname == 'surfaceairmodel_ocean_global'):
                    return global_satellite_effect(self.valid_indices[0])
                elif (groupname == 'surfaceairmodel_ice_global'):
                    return hemispheric_satellite_effect(self.valid_indices[0], self.location_polar_coordinates()[:,0])
            else:
                return None
        else:
            raise ValueError('Unrecognized bias effect groupname.')
        
    def check_all_candidates(self, pattern):
      
        candidates = [NameCandidate(name) for name in self.observable_filenames]
        result = nested_search_patterns([pattern], candidates)

        return (len(result) == len(candidates))
