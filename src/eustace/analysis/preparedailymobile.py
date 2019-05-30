"""EUMOPPS module to output raw binary files for subsequent processing using mobile location data (e.g. for ships)."""

import uuid
import argparse
import numpy
from eustace.timeutils.epoch import days_since_epoch
from fileio.observationsource_rawbinary import ObservationRawBinaryWriter
from fileio.observationsource_rawbinary import LocationLookupRawBinaryWriter
from fileio.observationsource_rawbinary import LocationLookupWithID
from eustace.preprocess.fileio.insitu_ocean import ObservationSourceInsituOceanHadNMAT2Default

SOURCE_INSITU_OCEAN = 'insitu_ocean'
SOURCECLASS = { 
    SOURCE_INSITU_OCEAN: ObservationSourceInsituOceanHadNMAT2Default
}

def run_day(outputobservationsfilename, outputlocationfilename, inputfilename, processdate, sourcename, observable):
    """EUMOPPS run commands."""

    # Get the info
    source = SOURCECLASS[sourcename](inputfilename)

    # Compute day number (since EUSTACE epoch)
    daynumber =  numpy.int64(days_since_epoch(processdate))
    
    # Retrieve daily data (and daily locations)
    dailydata = source.observations(observable)
    dailylocations = source.observation_location_lookup()

    # Location lookup structure with new unique ID
    dailylookup = LocationLookupWithID(uuid.uuid1(), dailylocations)

    # Store daily data
    ObservationRawBinaryWriter().write_day(outputobservationsfilename, dailylookup.uuid, dailydata, daynumber)

    # Store corresponding location lookup
    LocationLookupRawBinaryWriter().write(outputlocationfilename, dailylookup)
