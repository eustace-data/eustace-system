#!/bin/bash
#
# Run station data format verification
#

python2.7 -m eustace.outputformats.verify.verify_stations "$@"
