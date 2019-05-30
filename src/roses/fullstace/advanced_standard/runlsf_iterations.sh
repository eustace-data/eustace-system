#!/bin/bash
#
# runlsf_iterations.sh
#
# Multiple runs of runlsf.sh and runlsf_solver.sh using LSF wait features
# so that operation dependencies are satisfied
#

# Catalogue from command line
# The readlink command gives the absolute path so
# that jobs can use this even if run from other directories
EUMOPPS_CATALOGUE="$1"

# Iteration parameters
FIRST_ITERATION=0
LAST_ITERATION=2

# find full path of LSF scripts
# assuming in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
RUNLSF_PATHNAME="${SCRIPT_PATH}/runlsf.sh"
RUNLSF_SOLVER_PATHNAME="${SCRIPT_PATH}/runlsf_solver.sh"

# Call child scripts
# Each script should export its done condition as ${LSF_FINISHED_CONDITION}
# which is then fed as the wait condition for the next task

# Loop over iterations of system
for ITERATION in `seq ${FIRST_ITERATION} ${LAST_ITERATION}`; do

  source ${RUNLSF_PATHNAME} "${EUMOPPS_CATALOGUE}" "climatology_input_${ITERATION}" "${LSF_FINISHED_CONDITION}"
  source ${RUNLSF_SOLVER_PATHNAME} "${EUMOPPS_CATALOGUE}" "climatology_solve_${ITERATION}" "${LSF_FINISHED_CONDITION}"
  source ${RUNLSF_PATHNAME} "${EUMOPPS_CATALOGUE}" "large_scale_input_${ITERATION}" "${LSF_FINISHED_CONDITION}"
  source ${RUNLSF_SOLVER_PATHNAME} "${EUMOPPS_CATALOGUE}" "large_scale_solve_${ITERATION}" "${LSF_FINISHED_CONDITION}"
  source ${RUNLSF_PATHNAME} "${EUMOPPS_CATALOGUE}" "local_input_and_solve_${ITERATION}" "${LSF_FINISHED_CONDITION}"

done

# Final grid output
source ${RUNLSF_PATHNAME} "${EUMOPPS_CATALOGUE}" "output_grid" "${LSF_FINISHED_CONDITION}"
