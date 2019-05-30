#!/bin/bash
#
# runlsf_rawbinary.sh
#

# parameters
LSF_FILENAME='lsf_commandrun.sh'  # script to run

# find full path of LSF script
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
LSF_PATHNAME=${SCRIPT_PATH}'/'${LSF_FILENAME}

# Catalogue
export EUMOPPS_CATALOGUE="/work/scratch/joel/example_eustace/catalogue.nc"

# Module name from command line
export EUMOPPS_MODULENAME=$1

# Batch size to use
export EUMOPPS_BATCHSIZE=29

# Get batch count for this module name
export EUMOPPS_BATCHCOUNT=`python -m eumopps.catalogue.commandcount ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} --batchsize ${EUMOPPS_BATCHSIZE}`

# Tell the user what's going to happen
echo "Running: $LSF_PATHNAME for EUMOPPS module $EUMOPPS_MODULENAME with $EUMOPPS_BATCHCOUNT batches of $EUMOPPS_BATCHSIZE"

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
    bsub -J ${BATCHARRAY} < ${LSF_PATHNAME}

done
