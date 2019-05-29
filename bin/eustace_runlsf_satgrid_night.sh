#!/bin/bash
#
# eustace_runlsf_satgrid_night.sh
#
# Launch satgrid on parallel computing cluster.
# Options set for processing nighttime measurements.
#
#

# version control
SCRIPT_REVISION="\$Revision: 1336 $"
SCRIPT_NAME=`basename ${BASH_SOURCE[0]}`
SYSTEM_REVISION=`eustace_print_system_revision_id.sh`

# must have non-empty system revision
if [ -z $SYSTEM_REVISION ]; then
    exit 1
fi

# parameters
LSF_FILENAME='lsf_satgrid.sh'  # script to run
export STARTDATE="20100101"    # first date to process
export DAYS=366                # total days (including start date)
export RESOLUTION=0.25         # grid resolution (degrees)
export NLON=1440               # number of grid points for longitude
export NLAT=720                # number of grid points for latitude
export DAYNIGHT="night"        # day/night flag to include in filename
export QC_MASK_OBS=1           # QC bits that determine all observations of interest
export QC_MASK_VALID=7         # QC bits that determine valid observations
export QC_FILTER_OBS=1         # QC bit values at all observations of interest
export QC_FILTER_VALID=1       # QC bit values at valid observations

# use these optional parameters to override input path and filename patterns
# export BASEPATH="/gws/nopw/j04/eustace/data/incoming/MODIS"
# export PATTERN_LST="%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc"
# export PATTERN_AUX="%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc"

# version information to store in output
export OUTPUT_SOURCE=${SYSTEM_REVISION}' '${SCRIPT_NAME}' '${SCRIPT_REVISION}

# find full path of LSF script
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
LSF_PATHNAME=${SCRIPT_PATH}'/'${LSF_FILENAME}

# run
bsub -q short-serial -J "satgrid[1-$DAYS]" < ${LSF_PATHNAME}
