"""EUMOPPS module to generate location lookup tables where a single table for each surface is useful."""

import uuid
from fileio.observationsource_rawbinary import LocationLookupRawBinaryWriter
from fileio.observationsource_rawbinary import LocationLookupWithID
import preparedailyfixed
from observationsource import ObservationSource

def run(outputfilename, inputfilename, sourcename):
    """EUMOPPS run commands."""

    print 'Loading: ', inputfilename

    # load observation source (slightly different constructor for in-situ land)
    if sourcename == preparedailyfixed.SOURCE_INSITU_LAND:
        source = preparedailyfixed.SOURCECLASS[sourcename](inputfilename)
    else:
        source = preparedailyfixed.SOURCECLASS[sourcename]([ ObservationSource.TMIN ] , inputfilename)

    # store locations and newly-generated UUID
    lookup_with_uuid = LocationLookupWithID(uuid.uuid1(), source.observation_location_lookup())
    LocationLookupRawBinaryWriter().write(outputfilename, lookup_with_uuid)
