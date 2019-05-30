"""
Placeholder classes which can be specified in JSON configuration files and are
then expanded to describe per-operation parameters.
"""
   
from operationparameters import OperationParameterList
from operationparameters import OperationFileListReference
from operationparameters import OperationFileReference
import numpy
import itertools

class OperationException(ValueError):
    """This input exception is generated if inputs don't match the criteria described."""
    
    def __init__(self, message):
        
        super(OperationException, self).__init__(message)

class Input(object):
    """Base class for placeholders that represent inputs."""

    def operation_input_resolve(self, request_skip, catalogue, step):
        """Resolve placeholders for steps specified by step class.
           Optionally request to skip some steps."""

        raise NotImplementedError

class InputFileList(Input):
    """Describes the mapping from our position in the sequence to a required input.
       This includes a description of what happens if the input is missing."""

    MISSING_DATA_NOT_ALLOWED = 'not_allowed'
    MISSING_DATA_ALLOWED = 'allowed'
    MISSING_DATA_SKIP = 'skip'

    def __init__(self, datasetname, subsetindex=0, missing_data=MISSING_DATA_NOT_ALLOWED):
        """Construct reference to specified data set and subset.
           Optionally set policy for missing input data."""
        
        self.datasetname = datasetname
        self.subsetindex = subsetindex
        self.missing_data = missing_data

    def operation_input_resolve(self, request_skip, catalogue, step):
        
        return OperationFileListReference(self.find_references(request_skip, catalogue, step))

    def find_references(self, request_skip, catalogue, step):
        """Retrieve a list of DataReference objects that point to inputs on the given catalogue for each step defined by the OperationStep class.
           The request_skip list has elements set to True for those operation indices which should be skipped according to missing data policy."""
        
        # Retrieve subset
        subset = catalogue.datasubset(self.datasetname, self.subsetindex)

        # Must exist
        if not subset:
            raise OperationException('Missing input dataset {0} subset {1}'.format(self.datasetname, self.subsetindex))

        # List of empty lists (one per operation step)
        operation_data_references = [ [ ] for step_index in range(step.count()) ]

        # Need dataset index value for reference
        datasetindex = catalogue.datasetindex(self.datasetname)

        # Require the patterns attribute (as in DataStorageFiles)
        if not hasattr(subset.layout, 'patterns'):
            raise OperationException('Require subset layout with patterns attribute')

        # Must create map from operation steps to the relevant indices of inputs in the catalogue
        if step.is_uniquely_defined_by(subset.layout.patterns):

            # This is for the case where the relevant step can be uniqely identified from the input name
            #
            # For example, suppose data has 3 consecutive days crossing a year boundary and step is in years.
            # Suppose matches are like:
            #
            #   match 0: 2001-12-31.data
            #   match 1: 2002-01-01.data
            #   match 2: 2002-01-02.data
            # 
            # And our steps go like:
            #
            #   step 0: 2000
            #   step 1: 2001
            #   step 2: 2002
            #   step 3: 2003
            #
            #
            # Then match_with_step = [ ( 0, 1 ), ( 1, 2 ), ( 2, 2 ) ]
            # because the op group identifier is the index of the year.
            #
            match_with_step = [ (index, step.index_at_time(match.time))  for index, match in enumerate(subset.matches) ]

            # Lookup table where the key is an identifier for a group of operation steps that share the same matches
            # and the value is a list of the indices of files that match the pattern
            #
            # In the above example should get:
            #
            #     step 1 --> [ match 0 ]
            #     step 2 --> [ match 1, match 2 ]
            #
            #
            match_with_step = sorted(match_with_step, key=lambda(matchstep): matchstep[1])
            step_to_matchgroup = { step_index : [ matchstep[0] for matchstep in matchstep_iterator ] for step_index, matchstep_iterator in itertools.groupby(match_with_step, key=lambda(matchstep): matchstep[1]) }

        else:

            # This is for alternative case where knowing a match does not uniquely identify the step it belongs to
            #
            # Simple example is a fixed file that is used by time steps.
            #
            # Another example would be data in years and time steps are in days:
            #
            #   match 0: 2000.data
            #   match 1: 2001.data
            #   match 2: 2002.data
            #   match 3: 2003.data
            #
            #   step 0: 2001-12-31
            #   step 1: 2002-01-01
            #   step 2: 2002-01-02
            #
            # We use name as unique ID for such cases.
            #

            # Lookup of match index based on name
            matchname_to_index = { match.name: index  for index, match in enumerate(subset.matches) }

            # And find the operations that match this
            step_to_match = [ matchname_to_index.get(step.create_output_entry(subset.layout.patterns, step_index).name, None) for step_index in range(step.count()) ]

            # The same thing but as list object per operation to be consistent with the previous case above
            step_to_matchgroup = { step_index: [ match_index ] for step_index, match_index in enumerate(step_to_match) if match_index is not None }

        # Use the lookup to create references
        if self.missing_data == InputFileList.MISSING_DATA_ALLOWED:
            # indicator of -1 will later indicate that the file was missing
            default_entryindex = [-1]
        else:
            # empty entryindex indicates that there is no matching file and will be caught as an operation with missing data
            default_entryindex = []
            
        for step_index, data_reference_list in enumerate(operation_data_references):
            for entryindex in step_to_matchgroup.get(step_index, default_entryindex):
                data_reference_list.extend([ datasetindex, self.subsetindex, entryindex ])

        # Operation steps with no data
        operations_with_missing_data = [ operation_index for operation_index, references in enumerate(operation_data_references) if len(references) == 0 ]
        
        # Apply input policy
        if self.missing_data == InputFileList.MISSING_DATA_NOT_ALLOWED:

            # Missing data is not allowed - in this case every operation step must have one or more inputs
            if len(operations_with_missing_data) > 0:
                raise OperationException(
                    'Missing {0} inputs {1} (subset {2})'.format(
                        len(operations_with_missing_data),
                        self.datasetname,
                        self.subsetindex))

        elif self.missing_data == InputFileList.MISSING_DATA_SKIP:

            # Print warning
            print 'eumopps.catalogue.operation.OperationInput.resolve WARNING: missing input data from {0} for operations {1}'.format(
                self.datasetname, str(operations_with_missing_data))

            # Operation indices with missing data should be skipped
            for operation_index in operations_with_missing_data:
                request_skip[operation_index] = True

        return operation_data_references
        
class InputFile(InputFileList):
    """Describes the mapping from our position in the sequence to a required input.
       This includes a description of what happens if the input is missing."""

    def __init__(self, datasetname, subsetindex=0, missing_data=InputFileList.MISSING_DATA_NOT_ALLOWED):
        """Construct reference to specified data set and subset.
           Optionally set policy for missing input data.
           Optionally expect a list of multiple items per step."""

        super(InputFile, self).__init__(datasetname, subsetindex, missing_data)

    def operation_input_resolve(self, request_skip, catalogue, step):
        """Retrieve a list of DataReference objects that point to inputs on the given catalogue for each step defined by the OperationStep class.
           The request_skip list has elements set to True for those operation indices which should be skipped according to missing data policy."""
        
        # Evaluate references
        operation_data_references = self.find_references(request_skip, catalogue, step)

        # Can only have one item per step
        if  max([ len(ref) for ref in operation_data_references ]) > 3:
            raise OperationException('multiple items per step not allowed (if multiple files are expected, use class OperationInputFileList instead) [occurred on reference to {0}]'.format(self.datasetname))
        
        # Resolve
        return OperationFileReference(operation_data_references)

class StepTime(Input):

    def operation_input_resolve(self, request_skip, catalogue, step):
        
        return OperationParameterList([ step.time_at_index(operationindex) for operationindex in range(step.count()) ])

class StepIndex(Input):

    def operation_input_resolve(self, request_skip, catalogue, step):
        
        return OperationParameterList(numpy.array(range(step.count()), numpy.int32))

class Output(object):
    """Describes an operation output to create at each operation step."""

    def operation_output_resolve(self, catalogue, step, operation_indices):

        raise NotImplementedError

class OutputFile(Output):
    """Describes an operation output file to produce at each operation step."""
    
    def __init__(self, datasetname, subsetindex=0):
        
        self.datasetname = datasetname
        self.subsetindex = subsetindex

    def operation_output_resolve(self, catalogue, step, operation_indices):
        """Check output dataset has been defined (though should currently be empty) and append output entries to it."""

        # Retrieve subset
        subset = catalogue.datasubset(self.datasetname, self.subsetindex)

        # Must exist
        if not subset:
            raise OperationException('Missing output dataset {0} subset {1}'.format(self.datasetname, self.subsetindex))

        # Should be empty because each output can be made only once
        if len(subset.matches) > 0:
            raise OperationException('Operation requested additional outputs to a data subset that has already been defined')

        # Require the patterns attribute (as in DataStorageFiles)
        if not hasattr(subset.layout, 'patterns'):
            raise OperationException('Require subset layout with patterns attribute')

        # Record to avoid duplicates
        existing_name_index = { }

        # List of references to make (one trio of numbers per operation)
        references = [ ]

        # Need dataset index value for reference
        datasetindex = catalogue.datasetindex(self.datasetname)

        # New entry for each operation, but avoid duplicate names
        for operation_index in operation_indices:

            # Make new entry
            newentry = step.create_output_entry(subset.layout.patterns, operation_index)
        
            # Check if it already existed
            try:

                entryindex = existing_name_index[newentry.name]

            except KeyError:

                # Wasn't there - make a new match
                entryindex = len(subset.matches)
                subset.matches.append(newentry)
                existing_name_index[newentry.name] = entryindex

            # Append reference
            references.append([ datasetindex, self.subsetindex, entryindex ])
        
        # And return corresponding reference list
        return OperationFileReference(references)
  
