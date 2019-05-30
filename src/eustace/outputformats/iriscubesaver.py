"""Save Iris cubes as NetCDF with appropriate adaptations to meet EUSTACE specifications."""

from iris.fileformats import netcdf
import netCDF4
import time
from definitions import FILE_FORMAT
from definitions import FILE_CONVENTIONS

class VariableStorage(object):
    """Definition of storage required for specific variables."""

    def __init__(self, datatype, scale_factor, add_offset, fill_value):

        self.datatype = datatype
        self.scale_factor = scale_factor
        self.add_offset = add_offset
        self.fill_value = fill_value

class DatasetOverride(netCDF4.Dataset):
    """Override of netCDF4 dataset to intercept Iris choice of datatype."""

    def __init__(self, filename, **kwargs):

        super(DatasetOverride, self).__init__(filename, **kwargs)
        # self.storage_override = storage_override

    def createVariable(self, varname, datatype, dimensions=(), zlib=False, complevel=4, shuffle=True, fletcher32=False, contiguous=False, chunksizes=None, endian='native', least_significant_digit=None, fill_value=None):

        # Automatically on compression whenever non-zero complevel requested
        zlib = (complevel > 0)

        # Override storage information for known variable types
        override = NetCDFSaverEUSTACE.storage_override[varname] if varname in NetCDFSaverEUSTACE.storage_override else None
        if override:
            datatype = override.datatype
            fill_value = override.fill_value

        result = super(DatasetOverride, self).createVariable(varname, datatype, dimensions, zlib, complevel, shuffle, fletcher32, contiguous, chunksizes, endian, least_significant_digit, fill_value)

        if override:
            result.scale_factor = override.scale_factor
            result.add_offset = override.add_offset
            result.set_auto_maskandscale(True)

        return result

class NetCDFAttributesEUSTACE(object):
    """Place a standard set of attributes onto an iris cube list prior to saving."""

    def __init__(self, title, institution, comment, source):
        """Construct with specified attributes to place in cube."""

        self.attributes = {
            'title': title,
            'institution': institution,
            'history': 'Created ' + time.strftime('%c'),
            'comment': comment,
            'source': source
            }
        
    def apply(self, cubelist):
        for cube in cubelist:
            for key, value in self.attributes.iteritems():
                cube.attributes[key] = value
        
class NetCDFSaverEUSTACE(netcdf.Saver):
    """Adaptation of Iris NetCDF saving to meet EUSTACE specifications."""

    storage_override = { }
    """Global variable for storage types."""

    @staticmethod
    def set_storage(var_name, storage):
        NetCDFSaverEUSTACE.storage_override[var_name] = storage
    
    def __init__(self, filename):
        """Instantiate saving for specified filename."""

        # Call base class constructor
        super(NetCDFSaverEUSTACE, self).__init__(filename, FILE_FORMAT)

        # Replace default dataset loading with our overridden class
        self._dataset.close()
        self._dataset = DatasetOverride(filename, mode='w', format=FILE_FORMAT)


    def write_cubes(self, cubes):
        """Override writing code to impose required CF conventions string."""

        for cube in cubes:
            self.write(cube, unlimited_dimensions=[])

        self.update_global_attributes(Conventions=FILE_CONVENTIONS)

