"""Describe the layout of files in a data set and search for all available files."""

__version__ = "$Revision: 1019 $"
__author__ = "Joel R. Mitchelson"

import namepatterns
import os
import tarfile
import sys


class DataStorage(object):
    """Base class for describing the layout of files in a data set and searching for matching files."""

    def __init__(self):
        pass

    def match_files(self, catalogue_subset, datasetpath, candidates):
        """Search for matching files on the specified path (base class does nothing)."""
        pass


class DataStorageFiles(DataStorage):
    """Describe a data set in which a file is the smallest unit of storage."""

    CLASSID = 'files'
    """Identifier for kind of storage.  Used in descriptors and catalogue files."""

    def __init__(self, patterns):

        super(DataStorageFiles, self).__init__()
        self.patterns = patterns

    def match_files(self, catalogue_subset, datasetpath, candidates):
        """Search for matching files on the specified path."""

        catalogue_subset.append_file_results(namepatterns.nested_search_patterns(self.patterns, candidates))


class DataStorageArchive(DataStorage):
    """Base class for data set comprising archive files (such as zip), each containing multiple individual files."""

    CLASSID = 'archive'
    """Identifier for kind of storage.  Used in descriptors and catalogue files."""

    def __init__(self, archivepathname, internalpathname):

        super(DataStorageArchive, self).__init__()
        self.pathname = archivepathname
        self.internalpathname = internalpathname

    def match_files(self, catalogue_subset, datasetpath, candidates):
        """Search for matching files on the specified path and inside the archive files"""

        # match those corresponding to archive name pattern
        archive_results = namepatterns.nested_search_patterns([self.pathname], candidates)

        # search within archives
        for archive_entry in archive_results:

            # make full path to archive
            archive_pathname = os.path.join(datasetpath, archive_entry.name)
            sys.stderr.write(archive_pathname + '\n')

            # search archive to get internal file names
            internals = [namepatterns.NameCandidate(filename) for filename in self.get_archive_internals(archive_pathname)]
            # print internals

            # search filenames
            internal_results = namepatterns.nested_search_patterns([self.internalpathname], internals)

            # these ones were not touched
            internals_unused = [candidate.get_name() for candidate in internals if not candidate.is_matched()]

            # store results
            catalogue_subset.append_archive_contents(archive_entry.name, internal_results)
            catalogue_subset.append_archive_unused(archive_entry.name, internals_unused)

    def get_archive_internals(self, archive_pathname):
        """Get list of filenames in an archive - must be overidden in derived class for different archive types."""
        pass


class DataStorageTar(DataStorageArchive):
    """Data set comprising tar archives, each containing multiple files."""

    CLASSID = 'tar'
    """Identifier for kind of storage.  Used in descriptors and catalogue files."""

    def __init__(self, archivepathname, internalpathname):
        super(DataStorageTar, self).__init__(archivepathname, internalpathname)

    def get_archive_internals(self, archive_pathname):
        """Get list of filenames in .tar archive."""
        archive = tarfile.open(archive_pathname, self.tarmode())
        internals = [internal_file.name for internal_file in archive]
        archive.close()
        return internals

    def tarmode(self):
        """Mode used internally to open TAR archive."""
        return 'r'


class DataStorageTarGZip(DataStorageTar):
    """Data set comprising tar/gzip archives, each containing multiple files."""

    CLASSID = 'tar-gzip'
    """Identifier for kind of storage.  Used in descriptors and catalogue files."""

    def __init__(self, archivepathname, internalpathname):
        super(DataStorageTarGZip, self).__init__(archivepathname, internalpathname)

    def tarmode(self):
        """Mode used internally to open TAR archive (as tar-gzip in this case)."""
        return 'r:gz'

LAYOUTCLASSES = [DataStorageFiles, DataStorageTar, DataStorageTarGZip]
"""All layout classes."""

LAYOUTCLASSLOOKUP = {cls.CLASSID: cls for cls in LAYOUTCLASSES}


def layout_builder(name, parameters):
    """Build appropriate storage layout class given the parameters dictionary."""

    return LAYOUTCLASSLOOKUP[name](**parameters)
