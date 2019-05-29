#!/bin/bash

# retrieve CODE_PATH string from EUSTACE system configuration
CODE_PATH=`python2.7 -c "import eustaceconfig;print eustaceconfig.CODE_PATH"`

# discover and run tests on code path
python2.7 -m unittest discover -t $CODE_PATH -s $CODE_PATH/eumopps
