"""Input/Output of catalogue in NetCDF format."""

from eumopps.netcdfobjects import netcdfobjects
from eumopps.catalogue.catalogue import Catalogue
from reader import CatalogueReader
from writer import CatalogueWriter
from eumopps.catalogue.datareference import reference_integers_to_list
from netCDF4 import chartostring
import os.path

ATTRIBUTE_IDENTIFIER = 'identifier'
DATASETS_LIST_NAME = 'datasets'
OPERATIONS_LIST_NAME = 'operations'

class CatalogueReaderNetCDF(CatalogueReader):
    """Define class to encapsulate this reading operation."""

    def __init__(self):
        super(CatalogueReaderNetCDF, self).__init__()

    def load(self, pathname):
        """Load datasets where pathname corresponds to a NetCDF file."""

        ncfile = netcdfobjects.nc_open_for_load(pathname)
        catalogue = netcdfobjects.load_object(ncfile)
        ncfile.close()
        return catalogue

    def load_identifier(self, pathname):
        """Load unique identifier from catalogue."""

        ncfile = netcdfobjects.nc_open_for_load(pathname)
        result = ncfile.getncattr(ATTRIBUTE_IDENTIFIER)
        ncfile.close()
        return result

    def findgroups(self, ncfile, listname):
        """Quickly find the specified groupkey and value in the requested list of objects."""

        list_count = ncfile.dimensions[listname].size
        list_items = [ netcdfobjects.LISTITEM_OBJECT_SUFFIX.format(listname, index) for index in range(list_count) ]
        return [ ncfile.groups[itemname] for itemname in list_items ]

    def findoperation(self, ncfile, modulename):
        """Specified operation group object or None if no such group exists."""

        # All operations groups
        groups = self.findgroups(ncfile, OPERATIONS_LIST_NAME)

        # Filter on module name specified int the 'runmodule' member variable
        discovered = [ ]
        for group in groups:

            # Search by name, class, function name
            # We deliberately avoid elseif here in order to
            # detect ambiguous names
            if 'name' in group.ncattrs() and (group.getncattr('name') == modulename):
                discovered.append(group)
            runmodule = group['runmodule']
            if ('python_class' in runmodule.ncattrs()) and (runmodule.python_class == modulename):
                discovered.append(group)
            if ('python_function' in runmodule.ncattrs()) and (runmodule.python_function == modulename):
                discovered.append(group)

        # Return operation if found
        if len(discovered) == 1:
            return discovered[0]

        # Should be unique
        if len(discovered) > 1:
            raise ValueError('ambigious module name \"{0}\"'.format(modulename))

    def operationcount(self, pathname, modulename, batchsize=1):
        """Quickly establish number of records in NetCDF dataset, or number of batches if batchsize is > 1."""

        # Load file and get op group
        ncfile = netcdfobjects.nc_open_for_load(pathname)
        group = self.findoperation(ncfile, modulename)

        # Retrieve count value from operation
        record_count = group.getncattr('count') if group else 0

        # Done with file
        ncfile.close()

        # Express in batch numbers if requested
        batch_count = record_count / batchsize
        if (record_count % batchsize) != 0:
            batch_count = batch_count + 1

        return batch_count

    def load_operation(self, pathname, modulename):
        """Load the class associated with given module name."""

        ncfile = netcdfobjects.nc_open_for_load(pathname)
        group = self.findoperation(ncfile, modulename)
        result = netcdfobjects.load_object(group)
        ncfile.close()
        return result

    def load_data_pathnames(self, pathname, reflist):
        """Read path names based on list of data references."""

        ncfile = netcdfobjects.nc_open_for_load(pathname)
        names = []
        for ref in reflist:
            dataset = ncfile.groups[netcdfobjects.LISTITEM_OBJECT_SUFFIX.format('datasets', ref.dataset)]
            subset = dataset.groups[netcdfobjects.LISTITEM_OBJECT_SUFFIX.format('subsets', ref.subset)]
            matches = subset.groups['matches']
            
            # Construct pathname if the entry index is valid, otherwise no match is available in the catalogue
            if ref.entry>=0:
                name = str(chartostring(matches.variables['name'][ref.entry, :]))
                path = dataset.getncattr('path')
                names.append(os.path.join(path, name))
            else:
                names.append(None)
            #except:
            #    raise ValueError('ERROR: {0} references data which does not exist in the catalogue ({1}, {2}, {3})'
            #                     .format(__name__, ref.dataset, ref.subset, ref.entry))
        ncfile.close()
        return names

class CatalogueWriterNetCDF(CatalogueWriter):
    """Write catalogue using NetCDF format."""

    def __init__(self):

        super(CatalogueWriterNetCDF, self).__init__()

    def save(self, pathname, catalogue):
        """Save specified catalogue."""

        ncfile = netcdfobjects.nc_open_for_save(pathname)
        netcdfobjects.save_complex_python_object(ncfile, item=catalogue, membername=None)
        ncfile.close()
    
