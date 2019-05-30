"""Define field names for MODIS input data."""

import warnings
from satellite import SatelliteFieldNames

# Ignore Iris warnings about MODIS format
warnings.filterwarnings(action='ignore', message='Missing CF-netCDF auxiliary coordinate variable')
warnings.filterwarnings(action='ignore', message='Missing CF-netCDF label variable')
warnings.filterwarnings(action='ignore', message='NetCDF default loading behaviour currently does not expose variables which define reference surfaces')

class SatelliteFieldNamesMODIS(SatelliteFieldNames):
    """Define standard field names in MODIS data."""

    def __init__(self):

        super(SatelliteFieldNamesMODIS, self).__init__(
            primary='LST',
            qc='QC',
            coordinate_fields=['latitude','longitude'],
            uncertainty_fields=[ 'LST_unc_ran', 'LST_unc_loc_atm', 'LST_unc_loc_sfc' ],
            uncertainty_scalars=[ 'LST_unc_sys', ],
            # uncertainty_scalars=[ ],
            output_primary='ts',
            output_uncertainty=[ 'tsmean_unc_ran', 'tsmean_unc_loc_atm', 'tsmean_unc_loc_sfc' ],
            output_uncertainty_binned_scalars=['tsmean_unc_sys'])
            # output_uncertainty_binned_scalars=[ ])
