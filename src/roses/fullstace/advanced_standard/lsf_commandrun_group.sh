#!/bin/bash
#
#BSUB -c 06:00               # cpu time (hrs:mins)
#BSUB -W 06:00               # Wall time
#BSUB -n 1                   # number of tasks per array element
#BSUB -q short-serial
#BSUB -R "select[hname!='host157.jc.rl.ac.uk' && hname!='host291.jc.rl.ac.uk']"   # exclude known problem hosts

module load intel/cce/15.0.090
module load intel/mkl/11.3.1.150

EUMOPPS_COMMANDCOUNT=`python -m eumopps.catalogue.commandcount ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME}`
EUMOPPS_JOBINDEX=$((LSB_JOBINDEX - 1))                  # Indexes jobs from zero across all batches
EUMOPPS_SUBINDEX=$((EUMOPPS_JOBINDEX % EUMOPPS_NJOBS))  # Indexes jobs within this batch

SUPEROUTDIR=/work/scratch/$USER/lsf/${EUMOPPS_MODULENAME}/$EUMOPPS_BATCHNUMBER
mkdir -p ${SUPEROUTDIR}

OUTDIR=/work/scratch/$USER/lsf/${EUMOPPS_MODULENAME}/$EUMOPPS_BATCHNUMBER/$EUMOPPS_JOBINDEX
mkdir -p ${OUTDIR}

for EUMOPPS_TASKINDEX in `seq 0 $((EUMOPPS_NTASKS - 1))` ;
do

    EUMOPPS_INNERINDEX=$((EUMOPPS_SUBINDEX * EUMOPPS_NTASKS + EUMOPPS_TASKINDEX))              # Task index relative to the start of the batch
    EUMOPPS_COMMANDINDEX=$((EUMOPPS_BATCHNUMBER * EUMOPPS_BATCHSIZE + EUMOPPS_INNERINDEX))     # Task index relative to the start of batch zero
    
    # Break if we have exceeded the end of the total nuber of commands
    if [ $EUMOPPS_COMMANDINDEX -ge $EUMOPPS_COMMANDCOUNT ]; then
      break
    fi

    OUTFILE=${OUTDIR}/${EUMOPPS_TASKINDEX}.out
    ERRFILE=${OUTDIR}/${EUMOPPS_TASKINDEX}.err

    python -m eumopps.catalogue.commandrun ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} $EUMOPPS_INNERINDEX --use_batch_size $EUMOPPS_BATCHSIZE --batch $EUMOPPS_BATCHNUMBER --allow_unversioned_code 2> $ERRFILE > $OUTFILE

done
