#!/bin/bash
#
#BSUB -c 06:00               # cpu time (hrs:mins)
#BSUB -n 1                  # number of tasks per array element
#BSUB -W 04:00               # wall time
##BSUB -q eustace
#BSUB -q short-serial
#BSUB -R rusage[mem=4000]
#BSUB -M 4000
#BSUB -R "select[hname!='host157.jc.rl.ac.uk']"   # exclude known problem hosts

OPERATION_INDEX=$((LSB_JOBINDEX - 1))

mkdir -p ${LOGDIR}
LOGFILE=${LOGDIR}/${OPERATION_INDEX}.out
ERRFILE=${LOGDIR}/${OPERATION_INDEX}.err

python -m eustace.postprocess.merge_flag_files --reference_time_string $REFERENCE_TIME_STRING --operation_index $OPERATION_INDEX --input_directory_1 $INPUT_DIRECTORY_1 --input_directory_2 $INPUT_DIRECTORY_2 --output_directory $OUTPUT_DIRECTORY --iteration $ITERATION 2> $ERRFILE > $LOGFILE
