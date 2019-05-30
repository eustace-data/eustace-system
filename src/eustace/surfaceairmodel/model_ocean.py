"""
Surface-air ocean model
-----------------------
Python module for running land module (calls through to IDL)
"""

import netCDF4
import os
import subprocess
import json
import tempfile

class ModelOcean(object):
    """The wrapper for IDL core code."""
    def __init__(self, filename_model_model_clim_parameters, filename_model_clim_parameters_uncertainty, filename_model_STDEV_clim_parameters):
        """Initialise with locations of model files."""

        # Store parameters
        self.filename_model_model_clim_parameters = filename_model_model_clim_parameters
        self.filename_model_clim_parameters_uncertainty = filename_model_clim_parameters_uncertainty
        self.filename_model_STDEV_clim_parameters = filename_model_STDEV_clim_parameters

    @staticmethod
    def read_netcdf_time(filename):
        """Helper to get time from netcdf file."""

        dataset = netCDF4.Dataset(filename, 'r')
        timevariable = dataset.variables['time']
        return netCDF4.num2date(timevariable[0], units=timevariable.units)

    def build_config_dictionary(self, filenames_output, filenames_input):
        """Make configuration dictionary (suitable for JSON-encoding) using model filenames and the given daily input/output files."""

        # Dictionary for JSON-encoding
        jsondict = { 'model': { 'filename_parameters': self.filename_model_model_clim_parameters,
                                'filename_uncertainty': self.filename_model_clim_parameters_uncertainty,
                                'filename_stdev': self.filename_model_STDEV_clim_parameters },
                     'data': [ ] }

        # Build data part of dictonary for pairs of filenames
        for index, filename_output in enumerate(filenames_output):
           
            # Also need input filename
            filename_input = filenames_input[index]

            # file time
            time = ModelOcean.read_netcdf_time(filename_input)

            # dictionary for this pair of filenames
            dailydata = { 'date': { 'year': time.year, 'month': time.month, 'day': time.day },
                          'filename_output': filename_output,
                          'filename_input': filename_input }

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
        idlmodule = 'model_ocean_idl.pro'

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
        tempconfig = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.model_ocean.', suffix='.json')

        # Write config
        self.write_config_file(tempconfig.name, filenames_output, filenames_input)
            
        # Command line
        commandline = ModelOcean.build_idl_commandline(tempconfig.name)

        try:
            # Run in this directory
            process = subprocess.Popen(commandline, cwd=codepath)

        except OSError:
        
            raise OSError('ERROR launching IDL subprocess. Has IDL not been installed or module not loaded?')

        # Wait until complete
        result = process.wait()
        if result != 0:
            raise ValueError('ERROR running IDL script')
