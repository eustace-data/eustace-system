"""Compute checksum of file on disk."""

import subprocess
import numpy
import os


class Checksum(object):
    """Class to compute checksum for a file on disk."""

    def __init__(self, pathname):
        """Compute checksum and file size for file with given pathname.
           Results are stored as class members checksum (string containing checksum) and size (64-bit integer containing size)."""
        output = subprocess.Popen(['cksum', pathname], stdout=subprocess.PIPE, stderr=open(os.devnull, 'w')).communicate()
        checksum_string, size_string = output[0].split()[0:2]
        self.checksum = checksum_string
        self.size = numpy.int64(size_string)
