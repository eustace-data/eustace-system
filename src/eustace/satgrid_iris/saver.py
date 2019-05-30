"""Saving of output from satellite gridding."""

from eustace.outputformats.iriscubesaver import VariableStorage
from eustace.outputformats.iriscubesaver import NetCDFAttributesEUSTACE
from eustace.outputformats.iriscubesaver import NetCDFSaverEUSTACE
import numpy

# Declare storage for different units
STORAGE_K = VariableStorage(numpy.int16, 0.005, 273.15, numpy.int16(-32768))
STORAGE_RELATIVE_K = VariableStorage(numpy.int16, 0.001, 32.767, numpy.int16(-32768))
STORAGE_RELATIVE_K2 = VariableStorage(numpy.int32, 5.0E-6, 0.0, numpy.int32(-1))

# set global storage types
NetCDFSaverEUSTACE.set_storage('tsmean', STORAGE_K)
NetCDFSaverEUSTACE.set_storage('tsmin', STORAGE_K)
NetCDFSaverEUSTACE.set_storage('tsmax', STORAGE_K)
NetCDFSaverEUSTACE.set_storage('tsvariance', STORAGE_RELATIVE_K2)
NetCDFSaverEUSTACE.set_storage('tsmean_unc_ran', STORAGE_RELATIVE_K)
NetCDFSaverEUSTACE.set_storage('tsmean_unc_loc_atm', STORAGE_RELATIVE_K)
NetCDFSaverEUSTACE.set_storage('tsmean_unc_loc_sfc', STORAGE_RELATIVE_K)
NetCDFSaverEUSTACE.set_storage('tsmean_unc_spl', STORAGE_RELATIVE_K)
NetCDFSaverEUSTACE.set_storage('tsmean_unc_sys', STORAGE_RELATIVE_K)

def write_cubes(filename, source, cubelist):
    """Write cubelist to specified filename."""

    # apply attributes so that they appear in final NetCDF
    NetCDFAttributesEUSTACE(title='EUSTACE aggregated data', institution='Met Office', source=source, comment='').apply(cubelist)

    # use Iris to write NetCDF, with appropriate modifications for EUSTACE compliance
    NetCDFSaverEUSTACE(filename).write_cubes(cubelist)
