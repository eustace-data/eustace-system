"""
svn module
----------
Retrieve and parse global SVN repository revision number.
"""

import subprocess
import codepath

allow_unversioned_code = False  # pylint: disable=invalid-name
"""Global variable to bypass consistent SVN state."""

ID_PATTERN = 'R{:06d}'
SVN_STATE_NOT_CLEAN = 'Local copy of SVN has modifications or is in a mixed state.'
PROVENANCE_BYPASS_WARNING = 'UNVERSIONED'


class RepositoryStateException(Exception):
    """Error thrown when subversion repository is not in a clean state."""

    def __str__(self):
        """Returns a fixed error message."""
        return SVN_STATE_NOT_CLEAN


def set_allow_unversioned_code(flag):
    """Set global provenance bypass flag to permit mixed repository states."""
    globals()['allow_unversioned_code'] = flag


def get_revision_string(path_to_local_copy):
    """Get the revision string for the system (might have several numbers in it if in mixed state)."""
    proc = subprocess.Popen('svnversion', cwd=path_to_local_copy, stdout=subprocess.PIPE)
    return proc.stdout.read()


def parse_revision_string(revisionstring):
    """Return integer if s is a single latest revision, raise error otherwise."""
    try:
        return int(revisionstring)
    except Exception:
        raise RepositoryStateException()


def get_revision_id_for_number(number):
    """Convert revision number into ID."""
    return ID_PATTERN.format(number)


def get_revision_id(path_to_local_copy):
    """Get revision ID for the specified local copy or raise error if source is in a mixed state."""

    try:

        # Attempt to get revision ID
        return get_revision_id_for_number(parse_revision_string(get_revision_string(path_to_local_copy)))

    except RepositoryStateException as errorinstance:

        if allow_unversioned_code:

            # Swallow the exception if instructed to do so
            # and return default version indicator
            return PROVENANCE_BYPASS_WARNING

        else:

            # Repository was not in a clean state and no bypass
            # was requested, so pass repository exception on to caller
            raise errorinstance


def get_revision_id_for_module(pythonmodule):
    """Get revision ID for local copy in which specified python module is contained, or raise exception if source is in a mixed state.
       Optionally disable the exception if noprovenance is set True. In this case defaultmessage is returned instead of a version."""

    return get_revision_id(codepath.get_module_path(pythonmodule))
