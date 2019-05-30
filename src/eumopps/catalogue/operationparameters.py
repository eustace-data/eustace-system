"""Lists of parameters for operations stored on catalogue."""

from eumopps.catalogue.fileio.formatnetcdf import CatalogueReaderNetCDF
from datareference import reference_integers_to_list

class OperationParameter(object):
    """Base class for all operation parameters.  These may be global or per-operation."""

    def resolve_single_operation(self, cataloguepathname, operationindex):

        raise NotImplementedError

class OperationParameterList(OperationParameter):
    """Base class stored on final catalogues which lists parameters for multiple operations."""

    def __init__(self, operation_parameters):
        """Initialise with given list."""

        self.operation_parameters = operation_parameters

    def restrict_to_subset(self, restrict_to_operations):
        """Restrict to specified list of operation indices. 
           Update this object and return reference to ourselves."""

        self.operation_parameters = [ self.operation_parameters[index] for index in restrict_to_operations ]
        return self

    def resolve_single_operation(self, cataloguepathname, operationindex):
        """Extract parameters corresponding to specified operation index."""

        return self.operation_parameters[operationindex]

    def count(self):
        """Total operations represented."""
        
        return len(self.operation_parameters)
      
class OperationFileReference(OperationParameterList):
    """Intermediate class holding dataset references grouped by step index."""
    
    def __init__(self, operation_parameters):
        """Store reference."""
        
        super(OperationFileReference, self).__init__(operation_parameters)

    def resolve_single_operation(self, cataloguepathname, operationindex):

        reflist = reference_integers_to_list(self.operation_parameters[operationindex])
        pathnames = CatalogueReaderNetCDF().load_data_pathnames(cataloguepathname, reflist)
        if pathnames:
            return pathnames[0]

class OperationFileListReference(OperationParameterList):
    """Intermediate class holding dataset references grouped by step index."""
    
    def __init__(self, operation_parameters):
        """Store reference."""
        
        super(OperationFileListReference, self).__init__(operation_parameters)

    def resolve_single_operation(self, cataloguepathname, operationindex):

        reflist = reference_integers_to_list(self.operation_parameters[operationindex])
        return CatalogueReaderNetCDF().load_data_pathnames(cataloguepathname, reflist)

class OperationTimeStamp(OperationParameter):
    """Indicates that a timestamp of the operation will be output in place of this field.
       This class is checked for directly when expanding operations."""

    def __init__(self):
        pass

class OperationCatalogueID(OperationParameter):
    """Indicates that current catalogue ID will be output in place of this field."""

    def resolve_single_operation(self, cataloguepathname, operationindex):

        return CatalogueReaderNetCDF().load_identifier(cataloguepathname)
