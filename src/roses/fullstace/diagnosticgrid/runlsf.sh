#!/bin/bash
#
# runlsf.sh
#

# parameters
LSF_FILENAME=lsf_commandrun.sh  # script to run [lsf_commandrun.sh|lsf_commandrun_solver.sh]

# Catalogue from command line
# The readlink command gives the absolute path so
# that jobs can use this even if run from other directories
export EUMOPPS_CATALOGUE=`readlink -f $1`

# Module name from command line
export EUMOPPS_MODULENAME=$2

# Optional wait dependency to pass to bsub
# - see bsub documentation for more info
#   should be a string like "done(<job id>)"
LSF_WAIT_OPTIONS=$3

# find full path of LSF script
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
LSF_PATHNAME=${SCRIPT_PATH}'/'${LSF_FILENAME}

# Batch size to use
export EUMOPPS_BATCHSIZE=366

# Get batch count for this module name
export EUMOPPS_BATCHCOUNT=`python -m eumopps.catalogue.commandcount ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} --batchsize ${EUMOPPS_BATCHSIZE}`

# Tell the user what's going to happen
echo "Running: \"$LSF_PATHNAME\" for EUMOPPS module \"$EUMOPPS_MODULENAME\" with $EUMOPPS_BATCHCOUNT batches of $EUMOPPS_BATCHSIZE"
echo "Dependencies: \"${LSF_WAIT_OPTIONS}\""

# Run batches
for EUMOPPS_BATCHNUMBER in `seq 0 $((EUMOPPS_BATCHCOUNT - 1))` ;
do
    
    # Define this batch
    BATCHSTART=$((1 + EUMOPPS_BATCHNUMBER*EUMOPPS_BATCHSIZE))
    BATCHEND=$((BATCHSTART - 1 + EUMOPPS_BATCHSIZE))
    BATCHARRAY="${EUMOPPS_MODULENAME}_batch${EUMOPPS_BATCHNUMBER}[$BATCHSTART-$BATCHEND]"
    export EUMOPPS_BATCHNUMBER
    echo ${BATCHARRAY}

    # Run it
    bsub -w "${LSF_WAIT_OPTIONS}" -J ${BATCHARRAY} < ${LSF_PATHNAME}

done
