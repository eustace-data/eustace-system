"""Return base toolchain used by eumopps."""

from version.svn import get_revision_id_for_module
import eumopps


def versionlist():
    """List of versions used in toolchain."""

    return ['eumopps-' + get_revision_id_for_module(eumopps)]
