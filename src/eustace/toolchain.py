"""Module used by eumopps to generate toolchain versioning information."""

import os
import subprocess
from eumopps.version.svn import get_revision_id_for_module
import eustace

RECORD_JASMIN_TOOLCHAIN=True

def versionlist():
    """Get list of version identifiers for tools used."""

    tools = [ ]

    if RECORD_JASMIN_TOOLCHAIN:
        tools.append( commandline(['rpm', '-q', 'jasmin-common-vm']) )

    tools.append('eustace-' + get_revision_id_for_module(eustace))

    return tools

def commandline(arglist):
    """Run specified command line."""

    return subprocess.Popen(
        arglist,
        stdout=subprocess.PIPE,
        stderr=open(os.devnull, 'w')).communicate()[0].strip()
