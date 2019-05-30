"""
Ocean model analysis
--------------------
"""

__version__ = "$Revision: 1060 $"

from fileio.consistent import ConsistentModelOutputNetCDF
from fileio.consistent import DatasetAttributesConsistentModelOutput
from fileio.ocean import RemapNetCDFOceanWP1
from model_ocean import ModelOcean
from eustace.timeutils.epoch import EPOCH
from eustace.timeutils.epoch import days_since_epoch
from eustace.timeutils.epoch import epoch_plus_days
from eumopps.timeutils import datetime_numeric
import os.path
import argparse
import tempfile

def run_day(catalogue_id, institution, output_main, output_ancillary, input_model_parameters, input_model_uncertainty, input_model_stdev, input_satellite_ocean):
    """Run for the given set of daily files."""

    # Make model     
    model = ModelOcean(input_model_parameters, input_model_uncertainty, input_model_stdev)

    # Make temp files for intermediate results
    temporary_file = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.run_ocean', suffix='.nc')

    # Build intermediate results
    model.run_idl([temporary_file.name] , [input_satellite_ocean])

    # A remap class
    remapper = RemapNetCDFOceanWP1()

    # Attributes and writer for output
    attributes = DatasetAttributesConsistentModelOutput(__name__, institution, catalogue_id)
    writer = ConsistentModelOutputNetCDF(attributes, RemapNetCDFOceanWP1.ANCILLARY_VARIABLES)
        
    # Apply remapper and read fields
    fields = remapper.read_fields(temporary_file.name)

    # Write
    writer.write_primary(output_main, fields)
    writer.write_ancillary(output_ancillary, fields)

def return_daily_inputs(catalogue_id, institution, output_main, output_ancillary, input_model_parameters, input_model_uncertainty, input_model_stdev, input_satellite_ocean):
    """Collect single-day relevant information into a nested-keys dictionary"""

    details = {'catalogue_id':catalogue_id, 'institution':institution}
    input_parameters = {'input_model_parameters':input_model_parameters,
			'input_model_uncertainty':input_model_uncertainty, 
			'input_model_stdev':input_model_stdev,}
    input_files = {'input_satellite_ocean':input_satellite_ocean}

    output_files = {'output_main':output_main,'output_ancillary':output_ancillary}

    return {'details':details, 'input_parameters':input_parameters, 'input_files':input_files, 'output_files':output_files}

def run_multiple_days(list_of_daily_inputs):
    """Run for the given set of daily files."""
           
    outputs_main = []
    outputs_ancillary = []
    inputs_satellite_ocean = []
     
    # These information should be the same for every day
    catalogue_id = list_of_daily_inputs[0]['details']['catalogue_id']
    institution = list_of_daily_inputs[0]['details']['institution']
    input_model_parameters = list_of_daily_inputs[0]['input_parameters']['input_model_parameters']
    input_model_uncertainty = list_of_daily_inputs[0]['input_parameters']['input_model_uncertainty']
    input_model_stdev = list_of_daily_inputs[0]['input_parameters']['input_model_stdev']

    for dictionary in list_of_daily_inputs:
	inputs_satellite_ocean.append(dictionary['input_files']['input_satellite_ocean'])
	outputs_main.append(dictionary['output_files']['output_main'])	
	outputs_ancillary.append(dictionary['output_files']['output_ancillary'])
     
    # Make model     
    model = ModelOcean(input_model_parameters, input_model_uncertainty, input_model_stdev)

    # Make temp files for intermediate results
    tempfiles = [ tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.run_ocean.', suffix='.nc') for index in range(len(outputs_main)) ]

    # Build intermediate results
    model.run_idl([ t.name for t in tempfiles ], inputs_satellite_ocean)

    # A remap class
    remapper = RemapNetCDFOceanWP1()

    # Attributes and writer for output
    attributes = DatasetAttributesConsistentModelOutput(__name__, institution, catalogue_id)
    writer = ConsistentModelOutputNetCDF(attributes, RemapNetCDFOceanWP1.ANCILLARY_VARIABLES)

    # Read and remap
    for index, output_main in enumerate(outputs_main):
        
        # Apply remapper and read fields
        fields = remapper.read_fields(tempfiles[index].name)

        # Write
        writer.write_primary(output_main, fields)
        output_ancillary = outputs_ancillary[index]
        writer.write_ancillary(output_ancillary, fields)

def main():
    """Utility to run run_day from command line."""

    # parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument('--institution', required=True)
    parser.add_argument('--output_main', required=True)
    parser.add_argument('--output_ancillary')
    parser.add_argument('--input_model_parameters', required=True)
    parser.add_argument('--input_model_uncertainty', required=True)
    parser.add_argument('--input_model_stdev', required=True)
    parser.add_argument('--input_satellite_ocean', required=True,)

    args = parser.parse_args()

    run_day(catalogue_id='(NO CATALOGUE IN USE)',institution = args.institution,
	      output_main=args.output_main,output_ancillary=args.output_ancillary,
	      input_model_paramters=args.input_model_parameters,input_model_uncertainty=args.input_model_uncertainty,input_model_stdev=args.input_model_stdev,
	      input_satellite_land_day=args.inpu_day,input_satellite_land_night=args.inputnight,
	      input_satellite_ocean=args.input_satellite_ocean)

if __name__ == '__main__':
    main()
