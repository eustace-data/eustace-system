#!/bin/bash
#
# Paths for preparing catalogue of commands to run
#

#--------------------------------------------------------
# Configuration

# Spec of inputs
export INPUTCATALOGUE_JSON=inputs.json
export INPUTCATALOGUE_ICE_JSON=ice.json
export INPUTCATALOGUE_LAND_JSON=land.json
export INPUTCATALOGUE_OCEAN_JSON=ocean.json
# Catalogue to produce
export OUTPUTDIRECTORY=/gws/nopw/j04/eustace/data/internal/D2.2/R001335/20190111/
export INPUTCATALOGUE_NC=${OUTPUTDIRECTORY}/catalogue.nc

# Enable/disable enforcement of version control (enabled when blank)
#export VERSIONCONTROL=
#export VERSIONCONTROL=--allow_unversioned_code
#export CHECKSUM=--nochecksum