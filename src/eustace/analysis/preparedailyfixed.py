"""EUMOPPS module to output raw binary files for subsequent processing."""

import uuid
import argparse
import numpy
from eustace.timeutils.epoch import days_since_epoch
from fileio.observationsource_rawbinary import ObservationRawBinaryWriter
from fileio.observationsource_rawbinary import LocationLookupRawBinaryReader
from eustace.surfaceairmodel.fileio.consistent import ObservationSourceSatstace
from eustace.surfaceairmodel.fileio.lake import ObservationSourceLakeReading
from eustace.preprocess.fileio.insitu_land import ObservationSourceInsituLand
from eustace.preprocess.fileio.insitu_ocean import ObservationSourceInsituOceanHadNMAT2Default

SOURCE_SATELLITE_LAND = 'surfaceairmodel_land'
SOURCE_SATELLITE_OCEAN = 'surfaceairmodel_ocean'
SOURCE_SATELLITE_ICE = 'surfaceairmodel_ice'
SOURCE_SATELLITE_LAKE = 'surfaceairmodel_lake'
SOURCE_INSITU_LAND = 'insitu_land'
SOURCECLASS = { 
    SOURCE_SATELLITE_LAND: ObservationSourceSatstace,
    SOURCE_SATELLITE_OCEAN: ObservationSourceSatstace,
    SOURCE_SATELLITE_ICE: ObservationSourceSatstace,
    SOURCE_SATELLITE_LAKE: ObservationSourceLakeReading,
    SOURCE_INSITU_LAND: ObservationSourceInsituLand,
    }

def run_day(outputfilename, locationfilename, dailyfilename, processdate, sourcename, observable, hold_out_list=None, qc_mask=False, altitude_adjustment = None):
    """EUMOPPS run commands."""

    # First item is always the fixed location info
    reference_location = LocationLookupRawBinaryReader().read(locationfilename)

    # Compute day number (since EUSTACE epoch)
    daynumber =  numpy.int64(days_since_epoch(processdate))

    # Get the info
    if sourcename == SOURCE_INSITU_LAND:
        # special case for land to load only data relevant to given day number
        if hold_out_list == None:
            source = SOURCECLASS[sourcename](dailyfilename, daynumber, altitude_adjustment)
        else:
            source = SOURCECLASS[sourcename](dailyfilename, daynumber, hold_out_list, altitude_adjustment)
        
        # Load it
        dailydata = source.observations(observable, qc_mask)
    else:
        source = SOURCECLASS[sourcename]([ observable ] , dailyfilename)

        # Load it
        dailydata = source.observations(observable)

    # Verify location
    mylocation = source.observation_location_lookup()
    location_difference = (reference_location.lookuptable[:,dailydata.location] - mylocation[:,dailydata.location])
    if max(abs(location_difference.ravel())) > 0.0001:
        raise ValueError('location mismatch:\nexpected:\n{0}\ngot:\n{1}'
                         .format(str(reference_location.lookuptable[:,dailydata.location]), str(mylocation[:,dailydata.location])))

    # Store daily data
    ObservationRawBinaryWriter().write_day(outputfilename, reference_location.uuid, dailydata, daynumber)
    