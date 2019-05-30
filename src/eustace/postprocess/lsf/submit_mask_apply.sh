#!/bin/bash
#
# submit_mask_apply.sh
#
# Submit jobs to apply bit masks to analysis output
#

# parameters
LSF_FILENAME=masking_job_apply.sh  # script to run

export START_YEAR="$1"
export END_YEAR="$2"

export REFERENCE_TIME_STRING="${START_YEAR}-01-01"

export ITERATION=4
export ANALYSIS_DIRECTORY="/work/scratch/cmorice/advanced_standard.rawbinary8.1850-2015.R1400.20190329/"
export ANALYSIS_VERSION="R1400"
export OUTPUT_DIRECTORY="/gws/nopw/j04/eustace_vol2/data/internal/analysis/release_candidate_1.0/"

#export ITERATION=9
#export ANALYSIS_DIRECTORY="/work/scratch/cmorice/advanced_standard/"
#export ANALYSIS_VERSION="R1413"
#export OUTPUT_DIRECTORY="/gws/nopw/j04/eustace_vol2/data/internal/analysis/release_candidate_1.1/"

# threshold settings to load
export CLIMATOLOGY_THRESHOLD=0.7
export LARGESCALE_THRESHOLD=0.7
export CONSTRAINT_THRESHOLD=0.6

# windowing settings to load
export WINDOW_RANGE=10
export COUNT_THRESHOLD=5

export OUTLIER_WINDOW_RANGE=10
export OUTLIER_COUNT_THRESHOLD=10



export LOGDIR="/work/scratch/cmorice/lsf/mask_apply/${ANALYSIS_VERSION}/"
export FLAG_DIRECTORY="/work/scratch/cmorice/masking/${ANALYSIS_VERSION}/outlier_window_W${OUTLIER_WINDOW_RANGE}_C${OUTLIER_COUNT_THRESHOLD}_threshold_C${CLIMATOLOGY_THRESHOLD}_L${LARGESCALE_THRESHOLD}_T${CONSTRAINT_THRESHOLD}_window_W${WINDOW_RANGE}_C${COUNT_THRESHOLD}_calendar_extended/"
export FLAG_DIRECTORY2="/work/scratch/cmorice/masking/${ANALYSIS_VERSION}/location_threshold/"
#export OUTPUT_DIRECTORY="/work/scratch/cmorice/masking/${ANALYSIS_VERSION}/masked_outlier_window_W${OUTLIER_WINDOW_RANGE}_C${OUTLIER_COUNT_THRESHOLD}_threshold_C${CLIMATOLOGY_THRESHOLD}_L${LARGESCALE_THRESHOLD}_T${CONSTRAINT_THRESHOLD}_window_W${WINDOW_RANGE}_C${COUNT_THRESHOLD}_calendar/"

#export OUTPUT_DIRECTORY="/gws/nopw/j04/eustace_vol2/data/internal/analysis/${ANALYSIS_VERSION}/fullflag/"
#export OUTPUT_DIRECTORY="/gws/nopw/j04/eustace_vol2/data/internal/analysis/${ANALYSIS_VERSION}/nolargescaleflag/"


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
echo "Running: \"$LSF_PATHNAME\" for \"$START_YEAR\" to \"$END_YEAR\" with \"$N_OPERATIONS\" operations"
echo "Dependencies: \"${LSF_WAIT_OPTIONS}\""

# Build string of job ids indicating finished condition
export LSF_FINISHED_CONDITION=""

# Define this batch
BATCHSTART=1
BATCHEND=$N_OPERATIONS

BATCHARRAY="apply_masking[$BATCHSTART-$BATCHEND]"

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

