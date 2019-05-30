#!/bin/bash
#
# runlsf_locationpointers.sh
#

# parameters
LSF_FILENAME='lsf_locationpointers.sh'  # script to run

# Check command line parameters exist
if [ "$#" -lt 2 ];  then
    echo "Usage: runlsf_locationpointers.sh <cataloguefilename> <modulename> [--allow_unversioned_code]"
    exit 1
fi

# Catalogue
export EUMOPPS_CATALOGUE=$1

# Module name from command line
export EUMOPPS_MODULENAME=$2

# Options
export EUMOPPS_OPTIONS=$3

# find full path of LSF script
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
LSF_PATHNAME=${SCRIPT_PATH}'/'${LSF_FILENAME}

# Command count
export EUMOPPS_COMMANDCOUNT=`python -m eumopps.catalogue.commandcount ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME}`

# Batches
export BATCHARRAY="${EUMOPPS_MODULENAME}[1-$EUMOPPS_COMMANDCOUNT]"

# Tell the user what's going to happen
echo "Running: $LSF_PATHNAME for EUMOPPS module $EUMOPPS_MODULENAME with $EUMOPPS_COMMANDCOUNT commands"

# Run
bsub -M 64000000 -J ${BATCHARRAY} < ${LSF_PATHNAME}
