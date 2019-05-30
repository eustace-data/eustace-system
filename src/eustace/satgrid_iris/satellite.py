"""Filenames and loading of satellite data files."""

__version__ = "$Revision: 978 $"
__author__ = "Joel R. Mitchelson"

import os
import numpy
import sys
from datetime import datetime
from iris import load_cube
from iris import load
from iris import Constraint
from iris.coords import AuxCoord
from iris.cube import Cube

ERRORVALUE=-32768.


class SatelliteFilenames(list):
    """List of filename strings corresponding to a satelllite input source."""

    def check_readable(self):
        """Verify that all file components exist and have at least read access."""
        return all(SatelliteFilenames.check_component_readable(filename) for filename in self)

    @staticmethod
    def from_pattern_and_time(path, patterns, time):
        """Build names from path, patterns, and time, and return a SatelliteFilename object."""
        filenames = SatelliteFilenames()
        for pattern in patterns:
            filenames.append(SatelliteFilenames.build_component_filename(path, pattern, time))
        return filenames

    @staticmethod
    def build_component_filename(path, pattern, time):
        """Build one filename from path and pattern."""
        return os.path.join(path, time.strftime(pattern))

    @staticmethod
    def check_component_readable(filename):
        """Helper to check that one file exists."""
        return os.path.isfile(filename) and os.access(filename, os.R_OK)

class SatelliteInputList(object):
    """Represent a list of input sources built from filename pattern specification.

       Multiple filenames are needed for each input source to allow for uncertainty fields
       containted in auxiliary data files. Hence this is a 2D array in which rows are input sources 
       and columns are the individual filenames corresponding to each source."""

    TIME_MODE_SINGLE = 'single'
    TIME_MODE_DAY = 'day'

    def __init__(self, args):

        # initialise to empty
        times = []

        if len(args.date) == 12:
            # requested a specific time
            times.append(datetime.strptime(args.date, '%Y%m%d%H%M'))
            time_mode = SatelliteInputList.TIME_MODE_SINGLE
        else:
            # requested all times on a given day
            time_mode = SatelliteInputList.TIME_MODE_DAY
            selected_date = datetime.strptime(args.date, '%Y%m%d').date()
            for hour in range(0, 24):
                for minute in range(0, 60, 5):
                    times.append(datetime(selected_date.year, selected_date.month, selected_date.day, hour, minute))

        # store mode used
        self.time_mode = time_mode

        # build files from times
        self.filenames = [SatelliteFilenames.from_pattern_and_time(args.path, [ args.pattern_lst, args.pattern_aux ], t) for t in times]

class SatelliteFieldNames(object):
    """Names describing fields in satellite file."""

    def __init__(self, primary, qc, coordinate_fields, uncertainty_fields, uncertainty_scalars, output_primary, output_uncertainty, output_uncertainty_binned_scalars):
        self.primary = primary
        self.qc = qc
        self.coordinate_fields = coordinate_fields
        self.uncertainty_fields = uncertainty_fields
        self.uncertainty_scalars = uncertainty_scalars
        self.output_primary = output_primary
        self.output_uncertainty = output_uncertainty
        self.output_uncertainty_binned_scalars = output_uncertainty_binned_scalars

class SatelliteCollection(object):
    """Load satellite data into Iris cube and allow filtering based on QC flags."""

    @staticmethod
    def load(fields, filenames):
        """Load from specified filenames using the given field names instance."""

        # Constraint to indicate the main variable
        constraint_primary = Constraint(cube_func=lambda cube: cube.var_name==fields.primary)

        # Constraint to indicate aux variables representing 2D fields
        # (note coordinate fields sometimes loaded automatically by Iris, but included here anyway to be sure)
        fieldnames_aux = [ fields.qc ] + fields.coordinate_fields + fields.uncertainty_fields
        constraint_aux_fields = Constraint(cube_func=lambda cube: cube.var_name in fieldnames_aux)

        # Constraint to indicate aux variables representing scalar value
        constraint_aux_scalars = Constraint(cube_func=lambda cube: cube.var_name in fields.uncertainty_scalars)

        # Load cubes
        cube_primary = load_cube(filenames[0], constraint=constraint_primary)
        cubes_aux_field = load(filenames, constraints=[constraint_aux_fields])
        cubes_aux_scalar = load(filenames, constraints=[constraint_aux_scalars])
	
	# Append to primary cube as aux variables and check it covers known data
        for cube in cubes_aux_field:
            if hasattr(cube.data, 'mask'):
                missing_aux_pixels = numpy.sum(numpy.logical_and(numpy.logical_not(cube_primary.data.mask), cube.data.mask))
                if missing_aux_pixels:		    
                    original_pixel_count = numpy.sum(numpy.logical_not(cube_primary.data.mask))
                    message = 'Auxiliary field {0} does not cover primary field {1} with {2}/{3} pixels affected'.format(
                        cube.var_name, cube_primary.var_name, missing_aux_pixels, original_pixel_count)
                    # raise ValueError(message) # original behaviour was to generate an error here
                    # but instead restrict coverage to items with all fields valid
                    sys.stderr.write('WARNING: ' + message + ': removing those pixels\n')
                    cube_primary.data.mask |= cube.data.mask
            cube_primary.add_aux_coord(AuxCoord(points=cube.data, var_name=cube.var_name), (0,1))
	    
        # Append scalar values to primary cube (unless masked)
        for cube in cubes_aux_scalar:
            if hasattr(cube.data, 'mask'):
                if cube.data.mask[0]:
                    # Aux scalar is unexpectedly masked
                    pixels_affected = numpy.sum(numpy.logical_not(cube_primary.data.mask))
                    if pixels_affected:
                        # If there were supposed to be any valid measurement pixels, raise an error here
                        message = 'Auxiliary scalar {0} is masked out with {1} pixels affected'.format(cube.var_name, pixels_affected)
                        raise ValueError(message)
                    else:
                        # If there were no valid pixels anyway, just set to zero and report no error
                        cube.data.data.fill(0)
                        cube.data.mask.fill(False)
            cube_primary.add_aux_coord(AuxCoord(points=cube.data, var_name=cube.var_name))

        # Check required coordinates have been loaded
        # (may have happened automatically by Iris)
        for coordname in fields.coordinate_fields:
            try:
                auxcoord = next(coord for coord in cube_primary.coords() if coord.name() == coordname)
            except StopIteration:
                raise ValueError('Coordinate not found: ' + coordname)
            if hasattr(auxcoord.points.data, 'mask'):
                raise ValueError('Coordinate field {0} has masked elements'.format(coordname))

        # Create result
        return SatelliteCollection(fields, cube_primary)

    def __init__(self, fields, cube):
        """Create with specified field names and Iris cube."""

        self.fields = fields
        self.cube = cube
	self.safety_flag = False
	self.safety_mask = numpy.array([False])

    def apply_filter(self, filtermask):
        """Modify validity mask to exclude any observations where filtermask is True."""

        # Logical OR with existing validity mask, and store back to cube
        if hasattr(self.cube.data, 'mask'):
            self.cube.data.mask = numpy.logical_or(self.cube.data.mask, filtermask)
        else:
            self.cube.data = numpy.ma.masked_array(data=self.cube.data, mask=filtermask, dtype=self.cube.data.dtype)
    
    def generate_safety_mask(self):
	""" Check for invalid latitude-longitude values, create a safety mask to exclude them """
		
	for name in ['latitude', 'longitude']:
	    coord = self.cube.coord(name)
	    if (coord.points==ERRORVALUE).any():
	      self.safety_flag = True
	      coordinate_mask = coord.points==ERRORVALUE
	      message='Invalid {0} values: detected  {1}/{2} wrong values'.format(coord.name(),coordinate_mask.sum(),coordinate_mask.size)
	      sys.stderr.write('WARNING: ' + message + ': removing those pixels\n')
	      self.safety_mask = numpy.logical_or(self.safety_mask,coordinate_mask)
    
    def get_filter_from_qc_flags(self, qc_mask, qc_filter):
        """Examine QC flags and return array with non-zeros where (qc AND qc_mask) is NOT EQUAL to qc_filter."""

        # Read QC flags
        qc = self.cube.coords(self.fields.qc)[0].points
	
        # Apply mask
        masked = numpy.bitwise_and(qc, qc_mask)

	if self.safety_flag:
	  return numpy.logical_or((masked != qc_filter),self.safety_mask)
	else:
	  # Result is non-zero whenever flags don't match the filter
	  return (masked != qc_filter)

        
