#!/bin/bash
#
# test_buildall_lookups_land.sh
#
# Build lookup files for correlation ranges and fixed locations
#

# -- resticted to just land --

# Names of modules to run
declare -a EUMOPPS_MODULENAMES=(
"preparecorrelationranges_insitu_land_Tmin"
"preparecorrelationranges_insitu_land_Tmax"
#"preparecorrelationranges_insitu_ocean_Tmean"
#"preparecorrelationranges_surfaceairmodel_ice_Tmin"
#"preparecorrelationranges_surfaceairmodel_ice_Tmax"
#"preparecorrelationranges_surfaceairmodel_ice_Tmean"
#"preparecorrelationranges_surfaceairmodel_land_Tmin"
#"preparecorrelationranges_surfaceairmodel_land_Tmax"
#"preparecorrelationranges_surfaceairmodel_ocean_Tmean"
"preparelocationlookup_insitu_land"
#"preparelocationlookup_satellite"
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
echo "Eumopps Options"
echo " --- "

echo "Eumopps options: $EUMOPPS_OPTIONS"

echo " --- "
echo "Building lookups"
echo " --- "

# Run all modules
for EUMOPPS_MODULENAME in "${EUMOPPS_MODULENAMES[@]}"
do
  python -m eumopps.catalogue.commandrun ${EUMOPPS_OPTIONS} ${EUMOPPS_CATALOGUE} ${EUMOPPS_MODULENAME} 0
done
