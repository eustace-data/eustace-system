#!/bin/bash

# Load modules relevant for running the system
toLoad=(
    "idl"
    "intel/cce/15.0.090"
    "intel/mkl/11.3.1.150"
)

for i in "${toLoad[@]}"
do
   module load $i
done

