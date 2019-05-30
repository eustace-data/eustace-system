#!/bin/bash
#
# test_buildall_rawbinary.sh
#
# Submit all jobs to build raw binary files
#

# Command to call for each operation
COMMAND_FILENAME="test_runlsf_rawbinary.sh"

# Names of modules to run
declare -a EUMOPPS_MODULENAMES=(
"preparedailyfixed_insitu_land_Tmin"
"preparedailyfixed_insitu_land_Tmax"
#"preparedailymobile_insitu_ocean_Tmean"
#"preparedailyfixed_surfaceairmodel_ice_Tmin"
#"preparedailyfixed_surfaceairmodel_ice_Tmax"
#"preparedailyfixed_surfaceairmodel_ice_Tmean"
#"preparedailyfixed_surfaceairmodel_land_Tmin"
#"preparedailyfixed_surfaceairmodel_land_Tmax"
#"preparedailyfixed_surfaceairmodel_ocean_Tmean"   
)

# Check catalogue argument
if [ "$#" -lt 1 ]
  then
    echo "Usage: test_buildall_lookups_land.sh <cataloguefilename> [--allow_unversioned_code]"
    exit 1
fi

# Catalogue
EUMOPPS_CATALOGUE=$1

# Options
EUMOPPS_OPTIONS=$2

echo " --- "
echo "Eumopps options: $EUMOPPS_OPTIONS"
echo " --- "

echo " --- "
echo "Building raw binaries - buildall_rawbinary"
echo " --- "

# find full path of command to call
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
COMMAND=${SCRIPT_PATH}'/'${COMMAND_FILENAME}

# Run operations
for EUMOPPS_MODULENAME in "${EUMOPPS_MODULENAMES[@]}"
do
   echo "$SCRIPT_PATH"
   echo "$COMMAND"
   echo "$EUMOPPS_OPTIONS"
   source ${COMMAND} ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} ${EUMOPPS_OPTIONS} |& tee testoutput/${EUMOPPS_MODULENAME}.txt
done
