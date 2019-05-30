""" Write example files of raw binary observations.  Useful for testing from C++ side.
    Run like: python -m eustace.analysis.fileio.test.rawbinary_write_examples
    Outputs to current directory:
        example_rawbinary_localcorrelationrange.bin
        example_rawbinary_locationlookup.bin
        example_rawbinary_observations_57499.bin
        example_rawbinary_observations_57500.bin
        example_rawbinary_locationpointers.bin
"""

import uuid
from ...test.test_diagnosticgrid import SimulatedObservationSource
from ..observationsource_rawbinary import LocalCorrelationRangeRawBinaryWriter
from ..observationsource_rawbinary import LocationLookupWithID
from ..observationsource_rawbinary import LocationLookupRawBinaryWriter
from ..observationsource_rawbinary import ObservationRawBinaryWriter
from ..locationpointers import LocationPointersRawBinaryWriter

def main():
    """Entry point.  Write files."""

    # Will store this number in location lookup and observation files
    location_id = uuid.UUID('4996797e-f0dc-4495-b91d-4345189fcb71')

    # Make simulated source data
    source = SimulatedObservationSource()

    # Write correlation range file (two ranges)
    localcorrelationrange = source.local_correlation_length_scale('pretend')
    LocalCorrelationRangeRawBinaryWriter().write('example_rawbinary_localcorrelationrange.bin', localcorrelationrange)

    # Write location ID lookup file
    location_lookup = LocationLookupWithID(location_id, source.observation_location_lookup())
    LocationLookupRawBinaryWriter().write('example_rawbinary_locationlookup.bin', location_lookup)

    # Write observations file
    obs = source.observations('pretend')
    ObservationRawBinaryWriter().write_day('example_rawbinary_observations_57499.bin', location_id, obs, 57499)
    ObservationRawBinaryWriter().write_day('example_rawbinary_observations_57500.bin', location_id, obs, 57500)

    # Create location pointers file using the above data
    pointerwriter = LocationPointersRawBinaryWriter('example_rawbinary_locationlookup.bin', 'example_rawbinary_localcorrelationrange.bin')
    pointerwriter.append_day('example_rawbinary_observations_57499.bin', 57499)
    pointerwriter.append_day('example_rawbinary_observations_57500.bin', 57500)
    pointerwriter.write('example_rawbinary_locationpointers.bin')

if __name__ == "__main__":
    main()

