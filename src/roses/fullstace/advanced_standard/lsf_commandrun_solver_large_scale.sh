#!/bin/bash
#
#BSUB -c 12:00               # cpu time (hrs:mins)
#BSUB -n 16                  # number of tasks per array element
#BSUB -W 06:00               # wall time
#BSUB -q eustace
##BSUB -q par-single
#BSUB -R rusage[mem=64000]
#BSUB -M 64000
#BSUB -R "select[hname!='host157.jc.rl.ac.uk']"   # exclude known problem hosts

export MKL_NUM_THREADS=16
module load intel/cce/15.0.090
module load intel/mkl/11.3.1.150
EUMOPPS_COMMANDINDEX=$((LSB_JOBINDEX - 1))
EUMOPPS_SUBINDEX=$((EUMOPPS_COMMANDINDEX % EUMOPPS_BATCHSIZE))
OUTDIR=/work/scratch/$USER/lsf/${EUMOPPS_MODULENAME}/$EUMOPPS_BATCHNUMBER
mkdir -p ${OUTDIR}
OUTFILE=${OUTDIR}/${EUMOPPS_SUBINDEX}.out
ERRFILE=${OUTDIR}/${EUMOPPS_SUBINDEX}.err
python -m eumopps.catalogue.commandrun ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} $EUMOPPS_SUBINDEX --use_batch_size $EUMOPPS_BATCHSIZE --batch $EUMOPPS_BATCHNUMBER --allow_unversioned_code 2> $ERRFILE > $OUTFILE
