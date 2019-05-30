#!/bin/bash
#
# Build input catalogue
#

#--------------------------------------------------------
# Configuration

source environment.sh

#--------------------------------------------------------
# Ensure output directory exists

mkdir -p $OUTPUTDIRECTORY

#--------------------------------------------------------
# Make catalogue of inputs and operations

python2.7 -m eumopps.catalogue.build $INPUTCATALOGUE_ICE_JSON $INPUTCATALOGUE_NC $VERSIONCONTROL $CHECKSUM --pathdefault $OUTPUTDIRECTORY

python2.7 -m eumopps.catalogue.build $INPUTCATALOGUE_LAND_JSON $INPUTCATALOGUE_NC $VERSIONCONTROL $CHECKSUM --pathdefault $OUTPUTDIRECTORY --update

python2.7 -m eumopps.catalogue.build $INPUTCATALOGUE_OCEAN_JSON $INPUTCATALOGUE_NC $VERSIONCONTROL $CHECKSUM --pathdefault $OUTPUTDIRECTORY --update
