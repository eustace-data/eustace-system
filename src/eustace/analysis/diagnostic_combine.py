"""EUMOPPS module to run combination of gridded data."""

import numpy
from eustace.timeutils.epoch import days_since_epoch
from fileio.observationsource_rawbinary import LocationLookupRawBinaryReader
from fileio.observationsource_rawbinary import LocalCorrelationRangeRawBinaryReader
from fileio.observationsource_rawbinary import ObservationSourceSingleDayRawBinary
from fileio.observationsource_rawbinary import ObservableFileSpec
from diagnosticgrid import DiagnosticGridBins
from diagnosticgrid import DiagnosticGridObservationConnector
from eustace.outputformats.globalfield import GlobalFieldAxes2D
from eustace.outputformats.globalfield import GlobalFieldAxes2DWithTime
from eustace.outputformats.iriscubesaver import NetCDFSaverEUSTACE
from eustace.outputformats.ensuredirectory import ensuredirectory

class DiagnosticInput(object):
    """Class to hold description of each input source."""

    def __init__(self, sourcename, observable, correlationrangesfilename, locationfilename, observationfilename):

        self.sourcename = sourcename
        self.observable = observable
        self.correlationrangesfilename = correlationrangesfilename
        self.locationfilename = locationfilename
        self.observationfilename = observationfilename

def run(outputfilename, inputs, processdate):
    """EUMOPPS run commands."""

    # Need day number for inclusion in output
    daynumber = numpy.int64(days_since_epoch(processdate))

    # Grid onto this
    axes = GlobalFieldAxes2DWithTime(daynumber).aslist()
    outputgrid = DiagnosticGridBins(axes)

    # Cache location lookup as we might refer to same one multiple times
    locationlookup_cache = { }

    # One set of results per available source [land|sea|ice|lakes|in-situ land|in-situ ocean] per available observable [Tmean|Tmax|Tmin]
    for inputindex, descriptor in enumerate(inputs):
        
        # Cached loading of location lookup
        try:

            # Attempt to find in cache in case it's already been loaded for other sources/observables
            locationlookup = locationlookup_cache[descriptor.locationfilename]

        except KeyError:

            # Not found - load it for the first time this operation
            locationlookup = LocationLookupRawBinaryReader().read(descriptor.locationfilename)
            locationlookup_cache[descriptor.locationfilename] = locationlookup
        
        # Read correlation ranges for this item
        ranges = LocalCorrelationRangeRawBinaryReader().read(descriptor.correlationrangesfilename)

        # Observation files for each observable
        filespecs = { descriptor.observable: ObservableFileSpec(descriptor.observationfilename, ranges) }
            
        # Load as observation source            
        filesource = ObservationSourceSingleDayRawBinary(locationlookup, filespecs, daynumber)

        # Show stats
        print_stats(filesource, descriptor.sourcename, descriptor.observable)

        # Connector for gridding this
        connector = DiagnosticGridObservationConnector(axes, filesource)

        # Grid each observable
        dailydata = connector.get_day(descriptor.observable, daynumber)
        outputgrid.create_fields_from_sparse_observations(descriptor.sourcename, dailydata)
        outputgrid.compute_weighted_mean(descriptor.sourcename, descriptor.observable)

    # Store result
    ensuredirectory(outputfilename)
    saver = NetCDFSaverEUSTACE(outputfilename)
    saver.write_cubes(outputgrid)

def print_stats(source, sourcename, observable):
    """Write basic stats to standard output"""

    obs = source.observations(observable).measurement
    print 'min({0}:{1}): {2}'.format(sourcename, observable, min(obs))
    print 'max({0}:{1}): {2}'.format(sourcename, observable, max(obs))
