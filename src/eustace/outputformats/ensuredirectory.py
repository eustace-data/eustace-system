"""Ensure that parent directory for output files exists."""

import os
import errno

def ensuredirectory(pathname):
    """Attempt to make parent directory of the file pathname if it does not already exist."""

    # get subdirectory name
    subdirectory = os.path.dirname(pathname)

    # If none or empty we're using CWD, which doesn't need creating
    if subdirectory:

        # Use a try-except structure rather than attempting to check if it exists
        # because other processes might create it inbetween
        try:
            
            os.makedirs(subdirectory)

        except OSError as exception:

            # If it exists that's fine
            if exception.errno == errno.EEXIST and os.path.isdir(subdirectory):
                pass
            else:
                raise
