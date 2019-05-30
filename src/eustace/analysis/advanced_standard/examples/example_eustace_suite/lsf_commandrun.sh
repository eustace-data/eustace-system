#!/bin/bash
#
#BSUB -c 00:10               # cpu time (hrs:mins)
#BSUB -n 1                   # number of tasks per array element
#BSUB -q short-serial

module load intel/cce/15.0.090
module load intel/mkl/11.3.1.150
source /gws/nopw/j04/eustace/software/python-gpl/bin/activate
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gws/nopw/j04/eustace/software/SuiteSparse/lib
export LD_PRELOAD=$MKL_LIBS/libmkl_core.so:$MKL_LIBS/libmkl_sequential.so
export PYTHONPATH=$HOME/code/deveustace/system/src:$HOME/code/deveustace/research/metoffice-mockup/src:$HOME/code/deveustace/research/metoffice-analysis/eustace
export PATH=$HOME/code/deveustace/system/bin:$PATH:$HOME/.local/bin

EUMOPPS_COMMANDINDEX=$((LSB_JOBINDEX - 1))
EUMOPPS_SUBINDEX=$((EUMOPPS_COMMANDINDEX % EUMOPPS_BATCHSIZE))

OUTDIR=/work/scratch/joel/lsf/${EUMOPPS_MODULENAME}/$EUMOPPS_BATCHNUMBER
mkdir -p ${OUTDIR}
OUTFILE=${OUTDIR}/${EUMOPPS_SUBINDEX}.out
ERRFILE=${OUTDIR}/${EUMOPPS_SUBINDEX}.err
python -m eumopps.catalogue.commandrun ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} $EUMOPPS_SUBINDEX --use_batch_size $EUMOPPS_BATCHSIZE --batch $EUMOPPS_BATCHNUMBER --allow_unversioned_code 2> $ERRFILE > $OUTFILE
