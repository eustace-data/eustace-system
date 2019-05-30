#!/bin/bash
#
# test_create_rawbinary_files_land_noqc_ho.sh
#
# Run the full rawbinary creation process using limited set of input data as 
# specified in test .json files
# Options are set to speed up the tests = allow unversioned code, no checksum
# Data are only used for 2006 for the main multi year datasets
# Analysis is filtered by dataset for speed
# Index 0 = land, etc
# this script relates to modified copies of the main rawbinary generation files
# eg prefix test, and suffix  _land
 
# this script does the following:
# creates a new test directory for output 
# builds the catalogue - but with dataset filter of - insitu_land
# updates the catalogue using json files specifying only the insitu land operations 
# creates the lookup binaries
# creates the annual / daily rawbinaries
#  

DIRECTORY="/work/scratch/$USER/rawbinary_test_land_noqc_ho"
echo "$DIRECTORY"

if [ -d "$DIRECTORY" ]; then
  # Control will enter here if $DIRECTORY exists.
  rm -rf "$DIRECTORY"
  echo "Previous test run deleted"
fi

mkdir -p "$DIRECTORY"
echo "Test directory created in $DIRECTORY"

# run eumopps build - with test based parameters
# include options to allow unversioned code - for testing 

python -m eumopps.catalogue.build \
--dataset_filter "insitu_land" \
--nochecksum \
--allow_unversioned_code \
--pathdefault $DIRECTORY \
test_rawbinary_inputs.json \
$DIRECTORY/catalogue.nc

echo " --- "
echo "Catalogue file created"
echo " --- "

python -m eumopps.catalogue.build \
--nochecksum \
--allow_unversioned_code \
--pathdefault $DIRECTORY \
--update \
test_rawbinary_generation_land_noqc_with_holdout.json \
$DIRECTORY/catalogue.nc

echo " --- "
echo "Updated the catalogue"
echo " --- "

source test_buildall_lookups_land.sh \
$DIRECTORY/catalogue.nc \
--allow_unversioned_code \
|& tee testoutput/test_buildall_lookups_out_noqc_ho.txt

echo " --- "
echo "Built lookups"
echo " --- "

source test_buildall_rawbinary_land.sh \
$DIRECTORY/catalogue.nc \
--allow_unversioned_code \
|& tee testoutput/test_buildall_rawbinary_out_noqc_ho.txt

echo " --- "
echo "Jobs sent to Cylc for building rawbinary files"
echo " --- "

