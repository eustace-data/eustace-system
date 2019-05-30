#!/bin/bash
#
# buildall_locationpointers.sh
#
# Submit all jobs to build location pointers from raw binary files
#

# Command to call for each operation
COMMAND_FILENAME=runlsf_locationpointers.sh

# Names of modules to run
declare -a EUMOPPS_MODULENAMES=(
    "preparelocationpointers_insitu_land_Tmin"
    "preparelocationpointers_insitu_land_Tmax"
    "preparelocationpointers_surfaceairmodel_ice_Tmin"
    "preparelocationpointers_surfaceairmodel_ice_Tmax"
    "preparelocationpointers_surfaceairmodel_ice_Tmean"
    "preparelocationpointers_surfaceairmodel_land_Tmin"
    "preparelocationpointers_surfaceairmodel_land_Tmax"
    "preparelocationpointers_surfaceairmodel_ocean_Tmean"    
)

# Check catalogue argument
if [ "$#" -lt 1 ]
  then
    echo "Usage: buildall_locationpointers.sh <cataloguefilename> [--allow_unversioned_code]"
    exit 1
fi

# Catalogue
EUMOPPS_CATALOGUE=$1

# Options
EUMOPPS_OPTIONS=$2

# find full path of command to call
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
COMMAND=${SCRIPT_PATH}'/'${COMMAND_FILENAME}

# Run operations
for EUMOPPS_MODULENAME in "${EUMOPPS_MODULENAMES[@]}"
do
   source ${COMMAND} ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} ${EUMOPPS_OPTIONS}
done
