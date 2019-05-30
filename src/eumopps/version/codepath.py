"""
codepath module
---------------
Get the path from which a python module is loaded.
"""

import os


def get_module_path(pythonmodule):
    """Retrieve the full path from which the given python module is loaded."""
    return os.path.dirname(os.path.abspath(pythonmodule.__file__))
