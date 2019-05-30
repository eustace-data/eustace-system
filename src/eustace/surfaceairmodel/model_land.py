"""
Surface-air land model
----------------------
Python module for running land module (calls through to IDL)
"""

from datetime import datetime
import json
import os
import re
import subprocess
import tempfile

class ModelLand(object):
    """The wrapper for IDL core code."""
    def __init__(self, filename_model_clim_parameters, frac_sat_obs_threshold, sampling_threshold):
        """Initialise with locations of model files."""

        # Store parameters
        self.filename_model_clim_parameters = filename_model_clim_parameters
        self.frac_sat_obs_threshold = frac_sat_obs_threshold
        self.sampling_threshold = sampling_threshold

    @staticmethod
    def read_input_file_time(filename):
        """Helper to get time from input_file file."""

        pattern = re.search('(.[0-9]{4}[0-9]{2}[0-9]{2}.nc)',filename).group(0)
        date_string = re.search('([0-9]{4}[0-9]{2}[0-9]{2})',pattern).group(0)
        timevariable = datetime.strptime(date_string,'%Y%m%d')
        return timevariable

    def build_config_dictionary(self, filenames_output, filenames_input):
        """Make configuration dictionary (suitable for JSON-encoding) using model filenames and the given daily input/output files."""

        # Dictionary for JSON-encoding
        jsondict = { 'model': { 'filename_model_clim_parameters': self.filename_model_clim_parameters,
                                'frac_sat_obs_threshold': self.frac_sat_obs_threshold,
                                'sampling_threshold': self.sampling_threshold},
                     'data': [ ] }


        # Build data part of dictonary for pairs of filenames
        for index, filename_output in enumerate(filenames_output):
           
            # Also need information about input filenames and directories for fvc,dem,snow,lst
            directory_input_fvc = filenames_input[index][0]
            filename_input_dem = filenames_input[index][1]
            filename_input_mask = filenames_input[index][2]
            filename_input_snow = filenames_input[index][3]
            filename_input_lst = filenames_input[index][4]
            
            directory_output = filename_output[0]
	    filename_output_id = filename_output[1]
	   
            # file time
            time = ModelLand.read_input_file_time(filename_input_lst[0])
            timecheck = ModelLand.read_input_file_time(filename_input_lst[1])
            
            if(time != timecheck):
	      raise ValueError("LST day/night data files have mismatching dates!")

            # dictionary for this pair of filenames
            dailydata = { 'date': { 'year': time.year, 'month': time.month, 'day': time.day },
                          "dir_input_fvc": directory_input_fvc,
			  "filename_input_dem": filename_input_dem,
                          "filename_input_mask": filename_input_mask,
			  "filename_input_snow": filename_input_snow,
			  "filename_input_lst": filename_input_lst,
			  "filename_output_id": filename_output_id,
			  "dir_output": directory_output}
            # append
            jsondict['data'].append(dailydata)

        # Return this dictionary
        return jsondict
        
    def write_config_file(self, filename_config, filenames_output, filenames_input):
	"""Make configuration dictionary (suitable for JSON-encoding) using model filenames and the given daily input/output files."""

	# Make config dictionary
	config_dictionary = self.build_config_dictionary(filenames_output, filenames_input)

	# Write to file as JSON
	configfile = open(filename_config, 'w')
	json.dump(config_dictionary, configfile)
	configfile.close()

    @staticmethod
    def build_idl_commandline(filename_config):
        """Create command line for running IDL script against given config file."""

        # IDL module filename
        idlmodule = 'model_land_idl.pro'

        # Build command line
        return [ 'idl',
                 '-e', '.run ' + idlmodule,
                 '-args',
                 filename_config ]

    def run_idl(self, filenames_output, filenames_input):
        """Run on given input file."""

        # Directory of this module
        codepath = os.path.dirname(os.path.abspath(__file__))

        # A temporary file to store config
        tempconfig = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.model_land.', suffix='.json')

        # Write config
        self.write_config_file(tempconfig.name, filenames_output, filenames_input)
            
        # Command line
        commandline = ModelLand.build_idl_commandline(tempconfig.name)

        try:
            # Run in this directory
            process = subprocess.Popen(commandline, cwd=codepath)

        except OSError:
        
            raise OSError('ERROR launching IDL subprocess. Has IDL not been installed or module not loaded?')

        # Wait until complete
        result = process.wait()
        if result != 0:
            raise ValueError('ERROR running IDL script')
