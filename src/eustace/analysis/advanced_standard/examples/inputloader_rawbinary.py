
import os.path

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystemInputLoader
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceBinaryFilenameGenerator
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceSingleDayRawBinary
from eustace.analysis.fileio.observationsource_rawbinary import LocalCorrelationRangeRawBinaryReader
from eustace.analysis.fileio.observationsource_rawbinary import LocationLookupRawBinaryReader
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceBinaryFilenameGenerator
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceSingleDayRawBinary
from eustace.analysis.fileio.observationsource_rawbinary import ObservableFileSpec
#from eustace.analysis.fileio.observationsource_missing import ObservationSourceMissing
from eustace.timeutils.epoch import epoch_plus_days

from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations

from eustace.analysis.advanced_standard.fileio.observation_structure_source_connector import ObservationStructureSourceConnector

class ObservationSourceSingleDayRawBinary_ComputedMean(ObservationSourceSingleDayRawBinary):
    """Extend observation source to provide computed Tmean as (Tmax+Tmin)/2 where Tmean observations are not available directly."""

    def __init__(self, location_lookup_with_id, filespecs, daynumber):

        super(ObservationSourceSingleDayRawBinary_ComputedMean, self).__init__(location_lookup_with_id, filespecs, daynumber)

    def observations(self, observable):
        """Retrieve observations for specified observable quantity.
           Compute mean of max and min if mean is not available directly."""

        if observable in self.filespecs:

            # Where there is a file corresponding directly to observable then just use it
            return super(ObservationSourceSingleDayRawBinary_ComputedMean, self).observations(observable)

        elif (observable == ObservationSource.TMEAN and
              ObservationSource.TMIN in self.filespecs and 
              ObservationSource.TMAX in self.filespecs):

            # Otherwise return the computed mean from min and max
            tmin = super(ObservationSourceSingleDayRawBinary_ComputedMean, self).observations(ObservationSource.TMIN)
            tmax = super(ObservationSourceSingleDayRawBinary_ComputedMean, self).observations(ObservationSource.TMAX)
            return Observations.mean(tmin, tmax)
        else:
            print self.filespecs
            raise RuntimeError('Missing observable input')
        
    def local_correlation_length_scale(self, observable):
        """Get length scale."""

        # Change to TMIN if not available
        if (observable not in self.filespecs) and (observable == ObservationSource.TMEAN):
            observable = ObservationSource.TMIN

        # Lookup
        return super(ObservationSourceSingleDayRawBinary_ComputedMean, self).local_correlation_length_scale(observable)


class AnalysisSystemInputLoaderRawBinary(AnalysisSystemInputLoader):

    def __init__(self, local_correlation_ranges, observable_filenames, fixed_location_lookup=None, mobile_location_lookup_filenames=None, name=None):
        """Construct assuming correlation ranges and any fixed lookup already loaded."""

        self.local_correlation_ranges = local_correlation_ranges
        self.observable_filenames = observable_filenames
        self.fixed_location_lookup = fixed_location_lookup
        self.mobile_location_lookup_filenames = mobile_location_lookup_filenames
        self.name = None
        
    def datetime_at_time_index(self, time_index):
        """Express a time index as days since epoch as a datetime"""
        
        print time_index
        
        return epoch_plus_days(time_index)
    
    def load_observation_structure(self, observable, time_index, log=None):

        # Use fixed location lookup if we have it
        location_lookup = self.fixed_location_lookup

        # Otherwise need to load the mobile one
        if location_lookup is None and self.mobile_location_lookup_filenames is not None:

            # Use daily mobile location (e.g. for ships)
            if time_index in self.mobile_location_lookup_filenames:
                if self.mobile_location_lookup_filenames[time_index] is not None:
                    location_lookup = LocationLookupRawBinaryReader().read(self.mobile_location_lookup_filenames[time_index])
            
        # Observation files for each observable
        filespecs = { 
            observable_name: ObservableFileSpec(filename, self.local_correlation_ranges[observable_name])
            for observable_name, filename in self.observable_filenames[time_index].iteritems()
            #if filename is not None
            }
            
        # Load as observation source            
        observationsource = ObservationSourceSingleDayRawBinary_ComputedMean(location_lookup, filespecs, time_index)

	# Convert to observation stucture instance
        return ObservationStructureSourceConnector(
            observationsource=observationsource,
            observable=observable,
            corresponding_datetime=self.datetime_at_time_index(time_index))

class AnalysisSystemInputLoaderRawBinary_OneDay(AnalysisSystemInputLoaderRawBinary):
    """Load inputs to analysis system from raw binary files."""

    def __init__(self, 
                 local_correlation_ranges_filenames,
                 time_index,
                 observable_filenames,
                 fixed_location_lookup_filename=None,
                 mobile_location_lookup_filenames=None,
                 name=None):
        
        local_correlation_ranges = { observable: LocalCorrelationRangeRawBinaryReader().read(filename) for observable, filename in local_correlation_ranges_filenames.iteritems() }

        if fixed_location_lookup_filename is not None:            
            fixed_location_lookup = LocationLookupRawBinaryReader().read(fixed_location_lookup_filename)
        else:
            fixed_location_lookup = None

        super(AnalysisSystemInputLoaderRawBinary_OneDay, self).__init__(
            local_correlation_ranges,
            { time_index: observable_filenames },
            fixed_location_lookup,
            { time_index: mobile_location_lookup_filenames.values()[0] } if mobile_location_lookup_filenames else None)

class AnalysisSystemInputLoaderRawBinary_Sources(AnalysisSystemInputLoaderRawBinary):
    """Load inputs to analysis system from raw binary files."""

    # Pre-defined location lookup names for fixed sources
    LOCATIONLOOKUP_SATELLITE = 'locationlookup_satellite.bin'
    LOCATIONLOOKUP_INSITU_LAND = 'locationlookup_insitu_land.bin'
    FIXEDLOCATIONLOOKUP = {
        'surfaceairmodel_land': LOCATIONLOOKUP_SATELLITE,
        'surfaceairmodel_ocean': LOCATIONLOOKUP_SATELLITE,
        'surfaceairmodel_ice': LOCATIONLOOKUP_SATELLITE,
        'insitu_land': LOCATIONLOOKUP_INSITU_LAND
    }
    OBSERVABLES = {
        'surfaceairmodel_ice': [ ObservationSource.TMEAN ],
        'surfaceairmodel_land': [ ObservationSource.TMIN, ObservationSource.TMAX ],
        'surfaceairmodel_ocean': [ ObservationSource.TMEAN ],
        'insitu_land': [ ObservationSource.TMIN, ObservationSource.TMAX ],
        'insitu_ocean': [ ObservationSource.TMEAN ]
    }

    def __init__(self, basepath, source, time_indices):

        # Cache source
        self.source = source

        # Cache observable
        observables = AnalysisSystemInputLoaderRawBinary_Sources.OBSERVABLES[source]

        # Object to make filenames of interest        
        namebuilder = ObservationSourceBinaryFilenameGenerator(source, basepath)

        # Retrieve local correlation range values
        local_correlation_ranges = { observable: LocalCorrelationRangeRawBinaryReader().read(namebuilder.filename_local_correlation_ranges(observable)) for observable in observables }

        # Name for location lookup
        if source in AnalysisSystemInputLoaderRawBinary_Sources.FIXEDLOCATIONLOOKUP:

            # if in known fixed sources use fixed name
            location_lookup_filename = os.path.join(basepath, AnalysisSystemInputLoaderRawBinary_Sources.FIXEDLOCATIONLOOKUP[source])

            # read fixed lookup
            fixed_location_lookup = LocationLookupRawBinaryReader().read(location_lookup_filename)

            # no mobile lookup in this case
            mobile_location_lookup_filenames = None

        else:

            # No fixed location in this case
            fixed_location_lookup = None

            # New mobile lookup
            mobile_location_lookup_filenames = { 
                time_index: namebuilder.filename_mobile_locations(observables[0], time_index) for time_index in time_indices 
                }
            
        # Individual daily files
        observable_filenames = {
            time_index: { observable_name: namebuilder.filename_observations(observable_name, time_index)
                          for observable_name in observables }
            for time_index in time_indices
            }

        # Build
        super(AnalysisSystemInputLoaderRawBinary_Sources, self).__init__(
            local_correlation_ranges,
            observable_filenames,
            fixed_location_lookup,
            mobile_location_lookup_filenames)


    def load_observation_structure(self, observable, time_index, log=None):
    
        if log is not None:
            log.write('AnalysisSystemInputLoaderRawBinary.load_observation_structure: {0}\n'.format(self.source))

        return super(AnalysisSystemInputLoaderRawBinary_Sources, self).load_observation_structure(observable, time_index, log)

