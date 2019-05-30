#!/bin/bash
#
# create_rawbinary_files.sh
#
# Run the full rawbinary creation process  
# Use input specified in .json files
# Default is to use full range of input data


DIRECTORY="/work/scratch/$USER/rawbinary"

if [ -d "$DIRECTORY" ]; then
  # Control will enter here if $DIRECTORY exists.
  rm -rf "$DIRECTORY"
  echo "Previous run deleted"
fi

mkdir -p "$DIRECTORY"
echo "Directory created in $DIRECTORY"

# run eumopps build - 

# in production mode this will only permit use for versioned code
# for testing, flag can be modified to allow unversioned code
# for testing, nochecksum flag will speed up process 

python -m eumopps.catalogue.build \
--pathdefault $DIRECTORY \
# --nochecksum \
# --allow_unversioned_code \
rawbinary_inputs.json \
$DIRECTORY/catalogue.nc

echo " --- "
echo "Catalogue file created"
echo " --- "

python -m eumopps.catalogue.build \
--pathdefault $DIRECTORY \
# --nochecksum \
# --allow_unversioned_code \
--update \
rawbinary_generation.json \
$DIRECTORY/catalogue.nc

echo " --- "
echo "Updated the catalogue"
echo " --- "

source buildall_lookups.sh \
# --allow_unversioned_code \
$DIRECTORY/catalogue.nc 

echo " --- "
echo "Built lookups"
echo " --- "

source buildall_rawbinary.sh \
# --allow_unversioned_code \
$DIRECTORY/catalogue.nc 

echo " --- "
echo "Jobs sent to Cylc for building rawbinary files"
echo " --- "

