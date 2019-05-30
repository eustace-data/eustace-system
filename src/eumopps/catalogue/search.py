"""Search for all available files."""

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

import os


def find_all_files(basepath):
    """Make a list of all files on the specified path (returns path names relative to the base path)."""
    entries = []
    for dirpath, dirnames, dircontents in os.walk(basepath):
        for filename in dircontents:
            pathname = os.path.relpath(os.path.join(dirpath, filename), basepath)
            entries.append(pathname)
    return entries
