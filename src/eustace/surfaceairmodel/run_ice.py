"""
Ice model analysis
------------------
"""

__version__ = "$Revision: 1060 $"

from fileio.ice import IceSurfaceTemperatureQualityControlNetCDF
from fileio.consistent import ConsistentModelOutputNetCDF
from fileio.consistent import DatasetAttributesConsistentModelOutput
from model_ice import ModelIce
from eumopps.catalogue.operation import RunModule
from eustace.timeutils.epoch import EPOCH
from eustace.timeutils.epoch import days_since_epoch
from eustace.timeutils.epoch import epoch_plus_days
from eumopps.timeutils import datetime_numeric
import os.path
import argparse

def run_day(catalogue_id, institution, output_main, output_ancillary, input_model_sh, input_model_nh, input_satellite_sh, input_satellite_nh,modulename=None):
    """Run for the given set of daily files."""

    # Default module name
    if modulename is None:
        modulename = __name__

    # check output files do not already exist
    if os.path.isfile(output_main):
       raise ValueError('ERROR: file already exists: {0}'.format(output_main))
    if os.path.isfile(output_ancillary):
       raise ValueError('ERROR: file already exists: {0}'.format(output_ancillary))

    model = ModelIce({ ModelIce.SOUTHERN_HEMISPHERE: input_model_sh, ModelIce.NORTHERN_HEMISPHERE: input_model_nh })
    southern_hemisphere = IceSurfaceTemperatureQualityControlNetCDF(input_satellite_sh)
    northern_hemisphere = IceSurfaceTemperatureQualityControlNetCDF(input_satellite_nh)
    results = model.process_global(southern_hemisphere, northern_hemisphere)
    attributes = DatasetAttributesConsistentModelOutput(modulename, institution, catalogue_id)
    writer = ConsistentModelOutputNetCDF(attributes, ModelIce.ANCILLARY_VARIABLES)
    writer.write_primary(output_main, results)
    writer.write_ancillary(output_ancillary, results)

def main():
    """Utility to run from command line."""

    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--institution', required=True)
    parser.add_argument('--output_main', required=True)
    parser.add_argument('--output_ancillary')
    parser.add_argument('--input_model_nh', required=True)
    parser.add_argument('--input_model_sh', required=True)
    parser.add_argument('--input_satellite_nh', required=True)
    parser.add_argument('--input_satellite_sh', required=True)
    args = parser.parse_args()

    # run it
    run_day(modulename='eustace.surfaceairmodel.run_ice', catalogue_id='(NO CATALOGUE IN USE)', **args.__dict__).run()

if __name__ == '__main__':
    main()
