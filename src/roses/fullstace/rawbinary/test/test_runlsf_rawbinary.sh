#!/bin/bash
#
# test_runlsf_rawbinary.sh
#

# this script submits the jobs to Cycl
# batch size is set - for test purposes - reduced to batches of 10 


# Parameters
LSF_FILENAME='test_lsf_rawbinary.sh'  # script to run

# Batch size to use
export EUMOPPS_BATCHSIZE=10

# Check command line parameters exist
if [ "$#" -lt 2 ];  then
    echo "Usage: test_runlsf_rawbinary.sh <cataloguefilename> <modulename> [--allow_unversioned_code]"
    exit 1
fi

# Catalogue
export EUMOPPS_CATALOGUE=$1

# Module name from command line
export EUMOPPS_MODULENAME=$2

# Options to pass to EUMOPPS
export EUMOPPS_OPTIONS=$3

echo "Eumopps options: $EUMOPPS_OPTIONS"
echo "module: $EUMOPPS_MODULENAME"

# find full path of LSF script
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
LSF_PATHNAME=${SCRIPT_PATH}'/'${LSF_FILENAME}

# Get batch count for this module name
export EUMOPPS_BATCHCOUNT=`python -m eumopps.catalogue.commandcount ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} --batchsize ${EUMOPPS_BATCHSIZE}`

# Tell the user what's going to happen
echo "Running: $LSF_PATHNAME for EUMOPPS module $EUMOPPS_MODULENAME with $EUMOPPS_BATCHCOUNT batches of $EUMOPPS_BATCHSIZE"

echo "$EUMOPPS_CATALOGUE"
echo "$EUMOPPS_OPTIONS"

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
    bsub -r -J ${BATCHARRAY} < ${LSF_PATHNAME}

done
