#!/bin/bash
#
# runlsf_solver.sh
#

# Catalogue from command line
# The readlink command gives the absolute path so
# that jobs can use this even if run from other directories
export EUMOPPS_CATALOGUE=`readlink -f $1`

# Module name from command line
export EUMOPPS_MODULENAME="$2"

# Optional wait dependency to pass to bsub
# - see bsub documentation for more info
#   should be a string like "done(<job id>)"
LSF_WAIT_OPTIONS="$3"

# parameters
LSF_FILENAME=lsf_commandrun_solver.sh  # script to run [lsf_commandrun.sh|lsf_commandrun_solver.sh]

# find full path of LSF script
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
LSF_PATHNAME=${SCRIPT_PATH}'/'${LSF_FILENAME}

# Batch size to use
export EUMOPPS_BATCHSIZE=1

# Get batch count for this module name
export EUMOPPS_BATCHCOUNT=`python -m eumopps.catalogue.commandcount ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} --batchsize ${EUMOPPS_BATCHSIZE}`

# Tell the user what's going to happen
echo "Running: \"$LSF_PATHNAME\" for EUMOPPS module \"$EUMOPPS_MODULENAME\" with $EUMOPPS_BATCHCOUNT batches of $EUMOPPS_BATCHSIZE"

# Build string of job ids indicating finished condition
export LSF_FINISHED_CONDITION=""

# Run batches
for EUMOPPS_BATCHNUMBER in `seq 0 $((EUMOPPS_BATCHCOUNT - 1))` ;
do
    
    # Define this batch
    BATCHSTART=$((1 + EUMOPPS_BATCHNUMBER*EUMOPPS_BATCHSIZE))
    BATCHEND=$((BATCHSTART - 1 + EUMOPPS_BATCHSIZE))
    BATCHARRAY="${EUMOPPS_MODULENAME}_batch${EUMOPPS_BATCHNUMBER}[$BATCHSTART-$BATCHEND]"
    export EUMOPPS_BATCHNUMBER
    echo ${BATCHARRAY}

    # Run it and retrieve job information
    export JOBSTRING=`bsub -w "${LSF_WAIT_OPTIONS}" -M 64000000 -J ${BATCHARRAY} < ${LSF_PATHNAME}`

    # Parse job string using the sed utility
    # Job strings output from bsub look like:
    #
    #     "Job <1239736> is submitted to queue <short-serial>"
    #
    # And this will parse it so that it's just the string representing the number, e.g. "1239736"
    #
    export JOBID=`echo $JOBSTRING | sed -r "s/Job <([0-9]*)>.*/\1/"`

    # Build finished condition for whole batch, for example
    # "done(1239736) && done(1239737) && done(1239738)"
    if [ -n "${LSF_FINISHED_CONDITION}" ]; then
        export LSF_FINISHED_CONDITION="${LSF_FINISHED_CONDITION} && done(${JOBID})"
    else
        export LSF_FINISHED_CONDITION="done(${JOBID})"
    fi

done

# Output the finished condition which the next process can check
echo "Finished Condition: ${LSF_FINISHED_CONDITION}"
