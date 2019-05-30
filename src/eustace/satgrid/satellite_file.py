"""Filenames and loading of satellite data files."""

__version__ = "$Revision: 472 $"
__author__ = "Joel R. Mitchelson"

import os
import numpy
from netCDF4 import Dataset


class SatelliteFilename(object):
    """Represent LST and AUX file names in MODIS satellite data."""

    def __init__(self, filename_lst, filename_aux):
        """Construct from known names."""
        self.filename_lst = filename_lst
        self.filename_aux = filename_aux

    def check_readable(self):
        """Verify that all file components exist and have at least read access."""

        return SatelliteFilename.check_component_readable(self.filename_lst) and \
            SatelliteFilename.check_component_readable(self.filename_aux)

    @staticmethod
    def from_pattern_and_time(path, pattern_lst, pattern_aux, time):
        """Build names from path, patterns, and time, and return a SatelliteFilename object."""
        return SatelliteFilename(
            filename_lst=SatelliteFilename.build_component_filename(path, pattern_lst, time),
            filename_aux=SatelliteFilename.build_component_filename(path, pattern_aux, time))

    @staticmethod
    def build_component_filename(path, pattern, time):
        """Build one filename from path and pattern."""
        return os.path.join(path, time.strftime(pattern))

    @staticmethod
    def check_component_readable(component):
        """Helper to check that one file exists."""
        return os.path.isfile(component) and os.access(component, os.R_OK)


class SatelliteFile(object):
    """NetCDF file of satellite data together with corresponding file of auxillary data."""

    VARIABLES_CORE = ['LST', 'QC', 'lat', 'lon']
    VARIABLES_AUX = ['LST_unc_ran', 'LST_unc_loc_atm', 'LST_unc_loc_sfc', 'LST_unc_sys']
    VARIABLES = VARIABLES_CORE + VARIABLES_AUX

    def __init__(self, filename):
        """Open file specified by filename."""
        self.handle_lst = Dataset(filename.filename_lst, mode='r')
        self.handle_aux = Dataset(filename.filename_aux, mode='r')
        # self.handle.variables['LST'].set_auto_maskandscale(False)

    def get_variable(self, variablename, indices=None):
        """Get the specified variable from the appropriate component file. Name should be a member of VariableName."""
        if variablename not in SatelliteFile.VARIABLES:
            raise ValueError
        handle = self.handle_aux if variablename in SatelliteFile.VARIABLES_AUX else self.handle_lst
        return handle.variables[variablename][:] # pylint: disable=unsubscriptable-object

    def get_indices_from_qc(self, qc_mask, qc_filter):
        """Get relevant indices of data by filtering according to QC flags."""
        return SatelliteFile.find_qc_flag_indices(self.get_variable(SatelliteFile.VARIABLES[1]), qc_mask, qc_filter)

    def get_variable_at_indices(self, variablename, indices):
        """Retrieve the specified variable at specified indices."""
        values = self.get_variable(variablename)
        return values.flat[indices]

    @staticmethod
    def find_qc_flag_indices(qc, qc_mask, qc_filter):
        """Find integer patterns which are equal to qc_filter after masking with qc_mask."""

        # qc is a good name as it is a known abbreviation
        # pylint: disable=invalid-name

        masked = numpy.bitwise_and(qc, qc_mask)
        filtered = (masked == qc_filter)
        return numpy.flatnonzero(filtered)
