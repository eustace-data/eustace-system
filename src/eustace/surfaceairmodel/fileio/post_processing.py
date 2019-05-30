"""Collection of post-processing functions for any surface"""

from eustace.outputformats.definitions import TASMIN, TASMAX, TASMINUNCERTAINTY, TASMAXUNCERTAINTY
import numpy
import netCDF4

def apply_land_post_processing(results, surface_dependent_variables):
    """Post processing procedure applied on Satstace output of LSAT. Filter out pixels where Tmin >Tmax"""

    tasmindata = results[TASMIN.name][:]
    tasmaxdata = results[TASMAX.name][:]
    update_mask = numpy.logical_and(numpy.logical_and((tasmindata.data!=tasmindata.fill_value) ,(tasmaxdata.data!=tasmaxdata.fill_value)), tasmindata.data > tasmaxdata.data)

    results['problematic_pixels_count'] = update_mask.sum()
    for variable in [TASMIN, TASMINUNCERTAINTY, TASMAX, TASMAXUNCERTAINTY]+surface_dependent_variables:
      results[variable.name][:].mask = numpy.logical_or(results[variable.name][:].mask, update_mask)
      results[variable.name][:].data[:] = results[variable.name][:].filled()

