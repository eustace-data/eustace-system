"""EUMOPPS module to output raw binary files for subsequent processing."""

from fileio.observationsource_rawbinary import LocalCorrelationRangeRawBinaryWriter
import preparedailyfixed
import preparedailymobile

def run(outputfilename, inputfilename, sourcename, observable):
    """EUMOPPS run commands."""

    # Build a source class list that comes from daily and mobile modules (this module works with either)
    sourceclasses = { }
    sourceclasses.update(preparedailyfixed.SOURCECLASS)
    sourceclasses.update(preparedailymobile.SOURCECLASS)

    # Construct source reader (slightly different constructors for in-situ)
    if (sourcename == preparedailyfixed.SOURCE_INSITU_LAND) or (sourcename == preparedailymobile.SOURCE_INSITU_OCEAN):
        source = sourceclasses[sourcename](inputfilename)
    else:
        source = sourceclasses[sourcename]([ observable ], inputfilename)

    # Extract ranges
    local_correlation_ranges = source.local_correlation_length_scale(observable)

    # Write ranges
    rangewriter = LocalCorrelationRangeRawBinaryWriter()
    rangewriter.write(outputfilename, local_correlation_ranges)
