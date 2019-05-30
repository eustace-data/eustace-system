"""Describe a sequence of operations that create a new dataset."""

from eumopps import toolchain as systemtoolchain
import placeholder
from datareference import DataReference
from datareference import reference_list_to_integers
from operationparameters import OperationParameter
from operationparameters import OperationParameterList
from operationparameters import OperationTimeStamp
from operationparameters import OperationCatalogueID
from dataset import CatalogueFileEntry
import os.path
import importlib

FUNCTIONID = 'python_function'

class Operation(object):

    def __init__(self, runmodule, step, name=None, toolchain=None, newdatasets=None, count=None):
        """Initialise with module name to run, operation step class instance, and module parameters."""

        # We support runmodules which are classes or dictionaries containing a python method
        if not hasattr(runmodule, 'run'):
            if not (isinstance(runmodule, dict) and (FUNCTIONID in runmodule)):
                raise ValueError(
                    'runmodule must either be a class with run method or a dictionary containing a {0} key'.format(FUNCTIONID))

        self.runmodule = runmodule
        self.step = step
        self.name = name
        self.newdatasets = newdatasets
        self.count = count
        self.toolchain = toolchain


    def run(self):
        """Run the operation. This calls the run method of the runmodule or python_method if runmodule is a dictionary."""

        if hasattr(self.runmodule, 'run'):

            return runmodule.run()

        else:

            [ modulename, functionname ] = self.runmodule[FUNCTIONID].rsplit('.', 1)
            run_module = importlib.import_module(modulename)
            run_function = getattr(run_module, functionname)
            function_parameters = { key: value for key, value in self.runmodule.iteritems() if key != FUNCTIONID }
            
            print "Running:", repr(run_function)+' ('+repr(run_function.__module__)+')'
            print "With:"
            print function_parameters.keys()
            
            return run_function(**function_parameters)

    def current_runtime_toolchain(self):
        """Evaluate the current runtime toolchain of the runmodule."""

        # get module name in use
        if hasattr(self.runmodule, 'run'):
            modulename = self.runmodule.__module__
        else:
            modulename = self.runmodule[FUNCTIONID]

        # always use top-level parent package name
        basename = modulename.split('.')[0]

        # look for toolchain info in top-level
        processingtoolchain = importlib.import_module(basename + '.toolchain')

        # version list for entire toolchain (system and processing module)
        return systemtoolchain.versionlist() + processingtoolchain.versionlist()

    def resolve_single_operation(self, cataloguepathname, operationindex, timestamp):
        """After operation references are resolved, we can filter down to just one operation."""

        self.update_runmodule_parameters(OperationResolveSingleOperation(self.count, cataloguepathname, operationindex, timestamp))

    def resolve_operation_references(self, catalogue, pathdefault=None):
        """Update the catalogue with details of expected outputs, and return a modified parameters object
           in which placeholder.Input and placeholder.Output instances are replaced with lists of DataReferences
           on the catalogue. If newdatasets is non-empty, these will be added to the catalogue."""

        # Set our toolchain to be the current module toolchain
        self.toolchain = self.current_runtime_toolchain()

        # Make new empty datasets if requested and remove these from requestlist
        if self.newdatasets:

            # Set default path if none there
            for dataset in self.newdatasets:
                if not dataset.path:
                    if pathdefault:
                        dataset.path = pathdefault
                    else:
                        raise ValueError('missing path for dataset \"{0}\" and no pathdefault option specified'.format(dataset.name))

            # Append to catalogue
            catalogue.datasets.extend(self.newdatasets)

            # Remove from our list
            self.newdatasets = None

        # Store count value
        self.count = self.step.count()

        # Some inputs may request a skip of operations where no input is present
        # In this case the request_skip will be computed
        request_skip = [ False for operation_index in range(self.count) ]

        # Expand any input references 
        #
        # This will:
        #
        #     * Replace any placeholder.Input instance with a list of lists of DataReference objects
        #     * Populate the request_skip list if any inputs request a skip of the processing step
        #     * Raise an exception if data is missing when missing data is disallowed by the reference object
        #
        self.update_runmodule_parameters(OperationResolveInputReferences(self.step, request_skip, catalogue))

        # ...if we get this far it means no exception was raised due to missing inputs
        # so we can build outputs...

        # Will restrict to indices of available inputs (the ones we're not skipping)
        restrict_to_operations = [ operation_index for operation_index, skip in enumerate(request_skip) if not skip ]

        # Reduce count value
        self.count = len(restrict_to_operations)

        # Build outputs across the restricted set
        self.update_runmodule_parameters(OperationResolveOutputReferences(self.step, catalogue, restrict_to_operations))

        # Return updated runmodule
        return self.runmodule

    def update_runmodule_parameters(self, update):
        """Update all class members according to process() method of given update object.
           This will recurse into all member objects, lists, and dictionaries.
           The update object should be derived from OperationParameterUpdate."""

        self.runmodule = Operation.updatewalk(update, self.runmodule)

    @staticmethod
    def updatewalk(update, template_object):
        """Helper to update the template object using the process method of given update object."""

        # Attempt process
        processed = update.process(template_object)

        # Check result class 
        if isinstance(processed, OperationParameterUpdate.NotProcessed):

            # The update process did nothing - so just continue the walk over the tree
            if hasattr(template_object, '__dict__'):

                # Recurse into class members and resolve
                template_object.__dict__.update( Operation.updatewalk(update, template_object.__dict__) )

                # Then return the object
                return template_object

            elif isinstance(template_object, dict):
                    
                # Recurse into dictionaries
                return { key: Operation.updatewalk(update, element) for key, element in template_object.iteritems() }
                    
            elif isinstance(template_object, list):
                    
                # Iterate over lists
                return [ Operation.updatewalk(update, element) for element in template_object ]

            else:

                # Everything else treated as literal
                return template_object

        else:

            # The update process returned a result, so we should use that
            # Note that this might be None, as in the case when missing data files
            # are allowed
            return processed

class OperationParameterUpdate(object):
    """Base class for operations that replace placeholder parmeters of runmodule class."""

    class NotProcessed:
        """Return this class to indicate that an item was not processed."""
        pass

    def process(self, template_object):
        """Override this to do update."""
        
        raise NotImplementedError

class OperationResolveInputReferences(OperationParameterUpdate):

    def __init__(self, step, request_skip, catalogue):

        self.step = step
        self.request_skip = request_skip
        self.catalogue = catalogue

    def process(self, template_object):

        if isinstance(template_object, placeholder.Input):
                
            return template_object.operation_input_resolve(self.request_skip, self.catalogue, self.step)

        else:

            return OperationParameterUpdate.NotProcessed()

class OperationResolveOutputReferences(OperationParameterUpdate):

    def __init__(self, step, catalogue, restrict_to_operations):

        self.step = step
        self.catalogue = catalogue
        self.restrict_to_operations = restrict_to_operations

    def process(self, template_object):

        if isinstance(template_object, OperationParameterList):

            # Restrict any existing references with the requested subset
            return template_object.restrict_to_subset(self.restrict_to_operations)

        elif isinstance(template_object, placeholder.Output):

            # Generate output
            return template_object.operation_output_resolve(self.catalogue, self.step, self.restrict_to_operations)

        else:

            return OperationParameterUpdate.NotProcessed()

class OperationResolveSingleOperation(OperationParameterUpdate):
   
    def __init__(self, check_count, cataloguepathname, operationindex, timestamp):

        self.check_count = check_count
        self.cataloguepathname = cataloguepathname
        self.operationindex = operationindex
        self.timestamp = timestamp

    def process(self, template_object):

        if isinstance(template_object, OperationTimeStamp):

            # Output timestamp if given one by calling process
            return self.timestamp

        elif isinstance(template_object, OperationParameterList):

            # sanity-check: ensure number of items in input list matches the number of operations
            if self.check_count != template_object.count():
                raise ValueError('operation count mismatch')

            # resolve operation
            return template_object.resolve_single_operation(self.cataloguepathname, self.operationindex)

        elif isinstance(template_object, OperationParameter):

            # resolve operation
            return template_object.resolve_single_operation(self.cataloguepathname, self.operationindex)

        else:

            return OperationParameterUpdate.NotProcessed()

class RunModule(object):
    """Base class suitable for EUMOPPS module classes."""

    def run(self):
        """Derived classes should implement the run method."""

        raise NotImplementedError

