#!/bin/bash
#
# satgrid script to submit to parallel computing cluster.
#
# $Revision$ 
#
# NOTE: this outputs to the current working directory.
# Usually run from another script so that environment variables can be set.
#
# Run like: bsub -J "satgrid[1-366]" < lsf_satgrid.sh
#
# Requires environment variables:
#   OUTPUT_SOURCE
#   STARTDATE
#   DAYNIGHT
#   QC_MASK_OBS
#   QC_FILTER_OBS
#   QC_MASK_VALID
#   QC_FILTER_VALID
#   RESOLUTION
#   NLON
#   NLAT
#
# Optional environment variables (to override MODIS Aqua defaults):
#   BASEPATH
#   PATTERN_LST
#   PATTERN_AUX
#
#
#BSUB -o satgrid.%J.o
#BSUB -e satgrid.%J.e
#BSUB -c 00:10               # cpu time (hrs:mins)
#BSUB -n 1                   # number of tasks per array element
#

# defaults for filename patterns if not defined
: ${BASEPATH:="/gws/nopw/j04/eustace/data/incoming/MODIS"}
: ${PATTERN_LST:="%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_LST_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc"}
: ${PATTERN_AUX:="%Y/%m/%d/GT_MYG_2P/GT_SSD-L2-MYGSV_AUX_2-%Y%m%d_%H%M00-CUOL-0.01X0.01-V2.0.nc"}

# offset by LSF job array index
TIMESTRING=`date '+%C%y%m%d' -d "$STARTDATE + $LSB_JOBINDEX days - 1 days"`

# run command 
python2.7 -m eustace.satgrid_iris -source "${OUTPUT_SOURCE}" -qc_mask_obs=${QC_MASK_OBS} -qc_filter_obs=${QC_FILTER_OBS} -qc_mask_valid=${QC_MASK_VALID} -qc_filter_valid=${QC_FILTER_VALID} -xn ${NLON} -yn ${NLAT} -xs ${RESOLUTION} -ys ${RESOLUTION} -path ${BASEPATH} -pattern_lst ${PATTERN_LST}  -pattern_aux ${PATTERN_AUX} -o "satgrid.${DAYNIGHT}.${TIMESTRING}.nc" "${TIMESTRING}"
