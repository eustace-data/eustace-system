#!/bin/bash
#
# Check style of python code in eustace and eumopps packages against EUSTACE coding standards.
#

# retrieve CODE_PATH string from EUSTACE system configuration
SYSTEM_PATH=`python2.7 -c "import eustaceconfig;print eustaceconfig.SYSTEM_PATH"`
SOURCE_PATH=${SYSTEM_PATH}/src
RCFILE=${SYSTEM_PATH}/data/style/pylint.rc

# do check 
python2.7 -m pylint --rcfile ${RCFILE} eumopps eustace
