"""Direct inspection of a specific input using the raw binary files."""

import argparse
import json
import os
from datetime import datetime
from eustace.analysis.fileio.observationsource_rawbinary import LocalCorrelationRangeRawBinaryReader
from eustace.analysis.fileio.observationsource_rawbinary import LocationLookupRawBinaryReader
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceBinaryFilenameGenerator
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceSingleDayRawBinary
from eustace.analysis.fileio.observationsource_rawbinary import ObservableFileSpec
from eustace.timeutils.epoch import days_since_epoch

# Pre-defined location lookup names for fixed sources
LOCATIONLOOKUP_SATELLITE = 'locationlookup_satellite.bin'
LOCATIONLOOKUP_INSITU_LAND = 'locationlookup_insitu_land.bin'
FIXEDLOCATIONLOOKUP = {
    'surfaceairmodel_land': LOCATIONLOOKUP_SATELLITE,
    'surfaceairmodel_lake': LOCATIONLOOKUP_SATELLITE,
    'surfaceairmodel_ice': LOCATIONLOOKUP_SATELLITE,
    'insitu_land': LOCATIONLOOKUP_INSITU_LAND
}

def main():
    """EUSTACE inspectinput program
    
       Run like:
       
       python -m eustace.analysis.inspectinput --path /work/scratch/joel/diagnostic2 --source insitu_ocean --observable Tmean --date 19000101

    """
    

    # Parse command line arguments
    parser = argparse.ArgumentParser('inspectinput')
    parser.add_argument('--path', required=True, help='basepath')
    parser.add_argument('--source', required=True,
                        help='name of source [surfaceairmodel_land|surfaceairmodel_ocean|surfaceairmodel_ice|insitu_land|insitu_ocean]')
    parser.add_argument('--observable', required=True, help='name of observable [Tmin|Tmax|Tmean]')
    parser.add_argument('--date', required=True, help='date to inspect in form YYYYmmdd')
    args = parser.parse_args()

    # Convert datetime to daynumber since 01/01/1850
    daynumber = days_since_epoch(datetime.strptime(args.date, ObservationSourceBinaryFilenameGenerator.DATEFORMAT))

    # Object to make filenames of interest
    namebuilder = ObservationSourceBinaryFilenameGenerator(args.source, args.path)

    # Retrieve local correlation range values
    local_correlation_ranges = LocalCorrelationRangeRawBinaryReader().read(namebuilder.filename_local_correlation_ranges(args.observable))

    # Name for location lookup
    if args.source in FIXEDLOCATIONLOOKUP:
        # if in known fixed sources use fixed name
        location_lookup_filename = os.path.join(args.path, FIXEDLOCATIONLOOKUP[args.source])
    else:
        # otherwise use per-day location files (e.g. for ships)
        location_lookup_filename = namebuilder.filename_mobile_locations(args.observable, daynumber)

    # Read location lookup
    location_lookup = LocationLookupRawBinaryReader().read(location_lookup_filename)
    
    # Observation files for each observable
    # -- in this case we are looking at just one observable
    filespecs = { args.observable: ObservableFileSpec(namebuilder.filename_observations(args.observable, daynumber), local_correlation_ranges) }

    # Observation source
    observationsource = ObservationSourceSingleDayRawBinary(location_lookup, filespecs, daynumber)

    # Get observations for this day
    observations = observationsource.observations(args.observable)

    # Make a list of dictionaries we can print out
    results = [ ]
    for obs_index in range(observations.number_of_observations()):

        # Get measurement, uncorrelated uncertainty, and locally correlated components
        measurement = observations.measurement[obs_index]
        uncorrelatederror = observations.uncorrelatederror[obs_index]
        locallycorrelatederror = [ c[obs_index] for c in observations.locallycorrelatederror ]

        # Location ID number
        location_id = observations.location[obs_index]

        # Lookup for latitude and longitude
        coords = observationsource.observation_location_lookup()[:,location_id]

        # Append measurement and uncorrelated uncertainty
        results.append( { 'm': measurement, 'u':  uncorrelatederror, 'c': locallycorrelatederror, 'lat': coords[0], 'lon': coords[1] } )

    # Show results
    print json.dumps(results)


if __name__ == '__main__':
    main()
