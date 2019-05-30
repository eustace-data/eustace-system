"""
Land model analysis
-------------------
"""

__version__ = "$Revision: 1060 $"

from fileio.land import RemapNetCDFLandWP1
from fileio.consistent import ConsistentModelOutputNetCDF
from fileio.consistent import DatasetAttributesConsistentModelOutput
from model_land import ModelLand
from eustace.timeutils.epoch import EPOCH
from eustace.timeutils.epoch import days_since_epoch
from eustace.timeutils.epoch import epoch_plus_days
from eumopps.timeutils import datetime_numeric
import os.path
import argparse
import tempfile

def run_day(catalogue_id,institution,
	    output_main,output_ancillary,
	    input_model_1,input_model_2,input_model_3,
	    input_satellite_land_day,input_satellite_land_night,
	    input_satellite_land_snow,
	    input_satellite_land_mask_north,input_satellite_land_mask_south,
	    input_satellite_land_DEM,
	    input_satellite_land_fvc,
	    frac_sat_obs_threshold, 
	    sampling_threshold):
    """Run a single day analysis"""

    # Make model
    input_model_parameters = [input_model_1,input_model_2,input_model_3]
    model = ModelLand(input_model_parameters, frac_sat_obs_threshold, sampling_threshold)

    # Make temp file for intermediate results
    temporary_file = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.run_land.', suffix='.nc')

    input_mask = [input_satellite_land_mask_north,input_satellite_land_mask_south]
    inputs_day = [input_satellite_land_day,input_satellite_land_night]
    inputs_satellite = [[input_satellite_land_fvc, input_satellite_land_DEM, input_mask, input_satellite_land_snow, inputs_day]]

    
    # Build intermediate results
    model.run_idl([[tempfile.gettempdir()+'/' ,os.path.basename(temporary_file.name)]], inputs_satellite)
    
    # A remap class
    remapper = RemapNetCDFLandWP1()

    # Attributes and writer for output
    attributes = DatasetAttributesConsistentModelOutput(__name__, institution, catalogue_id)
    writer = ConsistentModelOutputNetCDF(attributes, RemapNetCDFLandWP1.VARIABLES_ANCILLARY)

    # Apply remapper and read fields
    fields = remapper.read_fields(temporary_file.name)

    # Filter out pixels where Tmin > Tmax
    writer.post_process(fields,'LAND')
    writer.dataset_attributes.comment = 'Percentage of (removed) Tmin > Tmax pixels: '+str(fields['problematic_pixels_count'])+'%'

    # Write
    writer.write_primary(output_main, fields)
    writer.write_ancillary(output_ancillary, fields)

def return_daily_inputs(catalogue_id,institution,
			output_main,output_ancillary,
			input_model_1,input_model_2,input_model_3,
			input_satellite_land_day,input_satellite_land_night,
			input_satellite_land_snow,
			input_satellite_land_mask_north,input_satellite_land_mask_south,
			input_satellite_land_DEM,
			input_satellite_land_fvc,
			frac_sat_obs_threshold, 
			sampling_threshold):
    """Collect single-day relevant information into a nested-keys dictionary"""

    details = {'catalogue_id':catalogue_id, 'institution':institution}
    input_parameters = {'input_model_parameters':[input_model_1,input_model_2,input_model_3],
			'frac_sat_obs_threshold':frac_sat_obs_threshold, 
			'sampling_threshold':sampling_threshold,}
    input_files = {'input_fvc':input_satellite_land_fvc, 
		   'input_DEM':input_satellite_land_DEM, 
		   'input_mask':[input_satellite_land_mask_north,input_satellite_land_mask_south],
		   'inputs_day':[input_satellite_land_day,input_satellite_land_night],
		   'input_snow':input_satellite_land_snow}
    output_files = {'output_main':output_main,'output_ancillary':output_ancillary}

    return {'details':details, 'input_parameters':input_parameters, 'input_files':input_files, 'output_files':output_files}

def run_multiple_days(list_of_daily_inputs):
    """Run for the given set of daily files."""

    outputs_main = []
    outputs_ancillary = []
    inputs_day = []
    inputs_snow = []
    
    # These information should be the same for every day
    catalogue_id = list_of_daily_inputs[0]['details']['catalogue_id']
    institution = list_of_daily_inputs[0]['details']['institution']
    input_model_parameters = list_of_daily_inputs[0]['input_parameters']['input_model_parameters']
    frac_sat_obs_threshold = list_of_daily_inputs[0]['input_parameters']['frac_sat_obs_threshold']
    sampling_threshold = list_of_daily_inputs[0]['input_parameters']['sampling_threshold']
    input_fvc = list_of_daily_inputs[0]['input_files']['input_fvc']
    input_DEM = list_of_daily_inputs[0]['input_files']['input_DEM']
    input_mask = list_of_daily_inputs[0]['input_files']['input_mask']

    for dictionary in list_of_daily_inputs:
	inputs_day.append(dictionary['input_files']['inputs_day'])
	if dictionary['input_files']['input_snow']!=None:
	    inputs_snow.append(dictionary['input_files']['input_snow'])
	else:
	    inputs_snow.append('MISSING')
	outputs_main.append(dictionary['output_files']['output_main'])	
	outputs_ancillary.append(dictionary['output_files']['output_ancillary'])

    # Make model
    model = ModelLand(input_model_parameters, frac_sat_obs_threshold, sampling_threshold)

    # Make temp files for intermediate results
    tempfiles = [ tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.run_land.', suffix='.nc') for index in range(len(outputs_main)) ]

    inputs_satellite = [[input_fvc, input_DEM, input_mask, inputs_snow[index], inputs_day[index]] for index in range(len(outputs_main))]
    
    # Build intermediate results
    model.run_idl([ [tempfile.gettempdir()+'/' ,os.path.basename(t.name)] for t in tempfiles ], inputs_satellite)
    
    # A remap class
    remapper = RemapNetCDFLandWP1()

    # Attributes and writer for output
    attributes = DatasetAttributesConsistentModelOutput(__name__, institution, catalogue_id)
    writer = ConsistentModelOutputNetCDF(attributes, RemapNetCDFLandWP1.VARIABLES_ANCILLARY)

    # Read and remap
    for index, output_main in enumerate(outputs_main):
        
        # Apply remapper and read fields
        fields = remapper.read_fields(tempfiles[index].name)
	
	# Filter out pixels where Tmin > Tmax
	writer.post_process(fields,'LAND')
	writer.dataset_attributes.comment = 'Number of (removed) Tmin > Tmax pixels: '+str(fields['problematic_pixels_count'])

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
    parser.add_argument('--input_model1_parameters', required=True)
    parser.add_argument('--input_model2_parameters', required=True)
    parser.add_argument('--input_model3_parameters', required=True)
    parser.add_argument('--inputfvc', required=True, help='directory of fvc dataset')
    parser.add_argument('--inputDEM', required=True, help='name of DEM datafile')
    parser.add_argument('--inputmasknord', required=True, help='nord mask datafile')
    parser.add_argument('--inputmasksouth', required=True, help='south mask datafile')
    parser.add_argument('--inputsnow', required=True, help='name of snow data file')
    parser.add_argument('--inputday', required=True, help='name of day input satgrid file')
    parser.add_argument('--inputnight', required=True, help='name of night input satgrid file')
    parser.add_argument('--frac_sat_obs_threshold', required=True,help='threshold for satellite observation')
    parser.add_argument('--sampling_threshold', required=True,help='sampling threshold')

    args = parser.parse_args()

    run_day(catalogue_id='(NO CATALOGUE IN USE)',institution = args.institution,
	      output_main=args.output_main,output_ancillary=args.output_ancillary,
	      input_model_1=args.input_model1_parameters,input_model_2=args.input_model2_parameters,input_model_3=args.input_model3_parameters,
	      input_satellite_land_day=args.inpu_day,input_satellite_land_night=args.inputnight,
	      input_satellite_land_snow=args.inputsnow,
	      input_satellite_land_mask_north=args.inputmasknord,input_satellite_land_mask_south=args.inputmasksouth,
	      input_satellite_land_DEM=args.inputDEM,
	      input_satellite_land_fvc=args.inputfvc,
	      frac_sat_obs_threshold=args.frac_sat_obs_threshold, 
	      sampling_threshold=args.sampling_threshold)

if __name__ == '__main__':
    main()
