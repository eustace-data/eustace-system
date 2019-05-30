"""Read the intermediate formats from WP1 land and ocean with the purpose of remapping to different variable names
   for output in consistent formats."""

import netCDF4
import numpy
from eustace.timeutils.epoch import days_since_epoch
from consistent import ConsistentModelOutputNetCDF
from eustace.outputformats import definitions

class RemapNetCDFSpecification(object):
    """Abstract base class for speciying a remap operation."""

    def __init__(self):
        """Base constructor does nothing."""
        pass

    def operate(self, results, netcdf):
        """Operate on netcdf dataset object and populate results dictionary."""

        raise NotImplementedError

class RemapNetCDFSpecDayNumber(RemapNetCDFSpecification):
    """Read time field and populate as day number since EUSTACE epoch."""

    def operate(self, results, netcdf):
        """Operate on netcdf dataset object and populate results dictionary."""

        timevariable = netcdf.variables['time']
        timevalue = netCDF4.num2date(timevariable[0], units=timevariable.units)
        results[ConsistentModelOutputNetCDF.FIELDNAME_DAYNUMBER] = int(days_since_epoch(timevalue))

class RemapNetCDFSpecCopy(RemapNetCDFSpecification):
    """Copy output to input with optional offset."""

    def __init__(self, outputname, inputname, offset=None):
        """Copy output to input with optional offset."""

        self.outputname = outputname
        self.inputname = inputname
        self.offset = offset

    def operate(self, results, netcdf):
        """Operate on netcdf dataset object and populate results dictionary."""

        # Retrieve field using input name
        field = netcdf.variables[self.inputname][:]

        # Apply offset if given (e.g. used for ocean to map Celsius -> Kelvin)
        if self.offset is not None:
            field += self.offset

        # Remove 1D items (because WP1 land format puts time at end rather than start of dimensions)
        field = numpy.squeeze(field)

        # Set using output name
        results[self.outputname] = field

class RemapNetCDFSpecTotalUncertainty(RemapNetCDFSpecification):
    """Read a number of uncertainty components and output total as square root of sum of squares."""

    def __init__(self, outputname, inputcomponents):
        """Read a number of uncertainty components and output total as square root of sum of squares."""

        self.outputname = outputname
        self.inputcomponents = inputcomponents

    def operate(self, results, netcdf):
        """Operate on netcdf dataset object and populate results dictionary."""

        # Build sum squares field
        sum_squares = numpy.ma.zeros(shape=definitions.GLOBAL_FIELD_SHAPE[1:])

        # Iterate over input components
        for componentname in self.inputcomponents:

            # Retrieve
            field = netcdf.variables[componentname][:]

            # Remove 1D items (because WP1 land format puts time at end rather than start of dimensions)
            field = numpy.squeeze(field)

            # Perform sum and propagate masks
            sum_squares += numpy.square(field)

        # Compute result
        result = numpy.sqrt( sum_squares )

        # Set sqrt of sum squares using output name
        results[self.outputname] = result


class RemapNetCDF(object):
    """Read NetCDF fields and place in dictionary with alternate names."""

    def __init__(self, mapping):
        """Construct with specified mapping (list of RemapNetCDFSpecification fields)."""

        self.mapping = mapping

    def read_fields(self, filename):
        """Read from filename and remap the fields."""

        # Empty results dictionary to populate
        results = { }

        # Open input NetCDF
        netcdf = netCDF4.Dataset(filename, 'r')
        # netcdf.set_auto_maskandscale(True)

        # Apply mapping and populate results
        for spec in self.mapping:
            
            spec.operate(results, netcdf)

        # Return populated results
        return results
