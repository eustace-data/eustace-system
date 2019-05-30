#!/bin/bash
#
#BSUB -c 00:15 #02:00               # cpu time (hrs:mins)
#BSUB -W 00:10 #02:00               # Wall time
#BSUB -n 1                   # number of tasks per array element
##BSUB -q eustace # short-serial
#BSUB -q short-serial
#BSUB -R "select[hname!='host157.jc.rl.ac.uk' && hname!='host291.jc.rl.ac.uk']"   # exclude known problem hosts

module load intel/cce/15.0.090
module load intel/mkl/11.3.1.150
EUMOPPS_COMMANDINDEX=$((LSB_JOBINDEX - 1))
EUMOPPS_SUBINDEX=$((EUMOPPS_COMMANDINDEX % EUMOPPS_BATCHSIZE))
OUTDIR=/work/scratch/$USER/lsf/${EUMOPPS_MODULENAME}/$EUMOPPS_BATCHNUMBER
mkdir -p ${OUTDIR}
OUTFILE=${OUTDIR}/${EUMOPPS_SUBINDEX}.out
ERRFILE=${OUTDIR}/${EUMOPPS_SUBINDEX}.err

python -m eumopps.catalogue.commandrun ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} $EUMOPPS_SUBINDEX --use_batch_size $EUMOPPS_BATCHSIZE --batch $EUMOPPS_BATCHNUMBER --allow_unversioned_code 2> $ERRFILE > $OUTFILE
