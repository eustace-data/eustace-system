#!/bin/bash
#
# Call EUMOPPS method to retrieve version of eustace system
#

SCRIPT_REVISION="\$Revision: 0 $"

python2.7 -c "from eumopps.version.svn import get_revision_id_for_module; import eustace; print get_revision_id_for_module(eustace)"
