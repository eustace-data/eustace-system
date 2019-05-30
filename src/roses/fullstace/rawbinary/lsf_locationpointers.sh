#!/bin/bash
#
#BSUB -o %J.out
#BSUB -c 02:00               # cpu time (hrs:mins)
#BSUB -W 02:00               # wall clock time (hrs:mins)
#BSUB -n 1                   # number of tasks per array element
#BSUB -q short-serial
#BSUB -R rusage[mem=62500]

EUMOPPS_COMMANDINDEX=$((LSB_JOBINDEX - 1))
OUTDIR=/work/scratch/${USER}/lsf/${EUMOPPS_MODULENAME}
mkdir -p ${OUTDIR}
OUTFILE=${OUTDIR}/${EUMOPPS_COMMANDINDEX}.out
ERRFILE=${OUTDIR}/${EUMOPPS_COMMANDINDEX}.err
python -m eumopps.catalogue.commandrun ${EUMOPPS_OPTIONS} ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} $EUMOPPS_COMMANDINDEX 2> $ERRFILE > $OUTFILE
