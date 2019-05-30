#!/bin/bash
#
# submit_cell_percentile_calculation.sh
#
# Submit jobs to derive percentiles for each grid cell to later be used in location checks
#

# parameters
LSF_FILENAME=cell_percentile_job.sh  # script to run

export START_YEAR=1981
export END_YEAR=2010

export ITERATION=4
export ANALYSIS_DIRECTORY="/work/scratch/cmorice/advanced_standard.rawbinary8.1850-2015.R1400.20190329/"
export ANALYSIS_VERSION="R1400"

#export ITERATION=9
#export ANALYSIS_DIRECTORY="/work/scratch/cmorice/advanced_standard/"
#export ANALYSIS_VERSION="R1413"

export INPUT_DIRECTORY=$ANALYSIS_DIRECTORY
export LOGDIR="/work/scratch/cmorice/lsf/cell_percentile_calculation/${ANALYSIS_VERSION}/"
export OUTPUT_DIRECTORY="/work/scratch/cmorice/masking/${ANALYSIS_VERSION}/cell_percentiles"


# Optional wait dependency to pass to bsub
# - see bsub documentation for more info
#   should be a string like "done(<job id>)"
LSF_WAIT_OPTIONS="$3"


# find full path of LSF script
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
LSF_PATHNAME=${SCRIPT_PATH}'/'${LSF_FILENAME}

export N_OPERATIONS=`python -m operation_count ${START_YEAR} ${END_YEAR}`

# Tell the user what's going to happen
echo "Running cell percentile calculation with 12 operations"
echo "Dependencies: \"${LSF_WAIT_OPTIONS}\""

# Build string of job ids indicating finished condition
export LSF_FINISHED_CONDITION=""

# Define this batch
BATCHSTART=1
BATCHEND=12

BATCHARRAY="cell_percentile_calc[$BATCHSTART-$BATCHEND]"

echo ${BATCHARRAY}

# Run it and retrieve job information
export JOBSTRING=`bsub -r -w "${LSF_WAIT_OPTIONS}" -J ${BATCHARRAY} < ${LSF_PATHNAME}`

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

# Output the finished condition which the next process can check
echo "Finished Condition: ${LSF_FINISHED_CONDITION}"

