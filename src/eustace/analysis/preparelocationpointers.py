"""EUMOPPS module to output raw binary files for subsequent processing."""

import numpy
from eustace.timeutils.epoch import days_since_epoch
from fileio.locationpointers import LocationPointersRawBinaryWriter
from datetime import datetime

def extractdaynumber(filename):
    """This is a hack to use filename to get daynumber, assuming filename ends with YYmmmdd.bin
       Ideally EUMOPPS would provide us with this but it doesn't at present."""

    # check file extension
    if filename[-4:] != '.bin':
        raise ValueError('Filename \"{0}\" expected to end with .bin but does not'.format(filename))

    # extract date string
    datestring = filename[-12:-4]

    # Convert to datetime object
    t = datetime.strptime(datestring, '%Y%m%d')

    # Convert to daynumber
    return numpy.int64( days_since_epoch(t) )


def run(outputfilename, locationfilename, correlationrangesfilename, dailyfilenames):
    """EUMOPPS run commands."""

    # Create writer object
    writer = LocationPointersRawBinaryWriter(locationfilename, correlationrangesfilename)

    # Compile pointer info in memory
    for input_index, inputfilename in enumerate(dailyfilenames):

        # Print progress
        if input_index % 100 == 0:
            pattern = '{0}: input {1} / {2}'
            message = pattern.format(__name__, input_index, len(dailyfilenames))
            print message

        # Append this info
        writer.append_day(inputfilename, extractdaynumber(inputfilename))

    # Store result
    writer.write(outputfilename)
    
