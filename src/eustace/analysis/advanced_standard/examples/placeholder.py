from eumopps.catalogue.placeholder import Input, InputFileList, Output
from eumopps.catalogue.operationparameters import OperationParameterList
from eustace.timeutils.epoch import days_since_epoch
import itertools
import numpy

import os.path

import datetime
import eustace.timeutils.epoch
import eumopps.timeutils.timebase
from eumopps.timeutils import datetime_numeric

from eumopps.catalogue.dataset import CatalogueFileEntry
from eumopps.catalogue.operationparameters import OperationFileReference, OperationFileListReference

def day_numbers_in_year(annual_step, index):
    
    d1 = int( eustace.timeutils.epoch.days_since_epoch( annual_step.time_at_index(index) ) )
    d2 = int( eustace.timeutils.epoch.days_since_epoch( annual_step.time_at_index(index+1) ) )
    
    return range( d1, d2 )

class AnalysisStepIndex(Input):
    """
    
    Allows operations to take integer valued time indices as days relative to 
    the analsis epoch.
    
    """

    def operation_input_resolve(self, request_skip, catalogue, step):
        
        return OperationParameterList( numpy.arange( int(days_since_epoch(step.start)),
                                                     int(days_since_epoch(step.end)+1),
                                                     1,
                                                     numpy.int32 ) )

class AnalysisAnnualBatchDayIndices(Input):
    """Return operation parameters for lists of day numbers within StepAnnual operations"""
    
    def operation_input_resolve(self, request_skip, catalogue, step):
        
        start_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.start) ) )
        end_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.end) ) )
        
        parameters = []
        for operationindex in range(step.count()):
            #day_numbers = numpy.array([day for day in day_numbers_in_year(step, operationindex) if day >= start_day_number and day <= end_day_number ], numpy.int32 )
            day_numbers = [day for day in day_numbers_in_year(step, operationindex) if day >= start_day_number and day <= end_day_number ]
            parameters.append( day_numbers )
        
        return OperationParameterList(parameters)

class AnalysisAnnualBatchDayTimes(Input):
    """Return operation parameters for lists of day numbers within StepAnnual operations"""
    
    def operation_input_resolve(self, request_skip, catalogue, step):
        
        timebase = eumopps.timeutils.timebase.TimeBaseDays(eustace.timeutils.epoch.EPOCH)
        
        start_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.start) ) )
        end_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.end) ) )
        
        parameters = []
        for operationindex in range(step.count()):
            #day_numbers = numpy.array([day for day in day_numbers_in_year(step, operationindex) if day >= start_day_number and day <= end_day_number ], numpy.int32 )
            day_numbers = [day for day in day_numbers_in_year(step, operationindex) if day >= start_day_number and day <= end_day_number ]
            
            datetimes = [timebase.number_to_datetime(daynumber) for daynumber in day_numbers]
            
            parameters.append( datetimes )
        
        return OperationParameterList(parameters)

class AnnualBatchDays(InputFileList):
    """
    
    Acts as InputFileList but returns references as lists of time ordered input files
    
    """

    def find_references(self, request_skip, catalogue, step):
        """Retrieve a list of DataReference objects that point to inputs on the given catalogue for each step defined by the OperationStep class.
           The request_skip list has elements set to True for those operation indices which should be skipped according to missing data policy."""
        
        # Retrieve subset
        subset = catalogue.datasubset(self.datasetname, self.subsetindex)
        
        print step
        print step.count()
        
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
            print "STEP UNIQUELY DEFINED BY PATTERN"
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
            match_with_step = [ (index, step.index_at_time(match.time), int( numpy.floor(eustace.timeutils.epoch.days_since_epoch(match.time))))  for index, match in enumerate(subset.matches) ]
            #print match_with_step
            
            #print match_with_step
            #print [ (index, int(eustace.timeutils.epoch.days_since_epoch(match.time)))  for index, match in enumerate(subset.matches) ]
            # Lookup table where the key is an identifier for a group of operation steps that share the same matches
            # and the value is a list of the indices of files that match the pattern
            #
            # In the above example should get:
            #
            #     step 1 --> [ match 0 ]
            #     step 2 --> [ match 1, match 2 ]
            #
            #
            match_with_step = sorted(match_with_step, key=lambda(matchstep): matchstep[1]) # presort so that itertools.groupby takes contiguous step indices for grouping
            step_to_matchgroup = { step_index : [ (matchstep[0], matchstep[2]) for matchstep in matchstep_iterator ] for step_index, matchstep_iterator in itertools.groupby(match_with_step, key=lambda(matchstep): matchstep[1]) }
            
            #print "step_to_matchgroup"
            #print step_to_matchgroup
            
            start_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.start) ) )
            end_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.end) ) )
            
            # For each step, map each matching file to its day in the year. Fill with missing_input_indicator for days with no matching input.
            missing_input_indicator = -1
            for step_index, matchdaypairs in step_to_matchgroup.iteritems():
                
                # unpack pair of file indices and their corresponding day numbers for this step
                matches, day_numbers = zip(*matchdaypairs)
                
                valid_indices = [i for i in range(len(matches)) if day_numbers[i] >= start_day_number and day_numbers[i] <= end_day_number ]
                
                matches = [matches[i] for i in valid_indices ]
                day_numbers = [day_numbers[i] for i in valid_indices ]
                
                #matches, day_numbers = [match, day for match, day in zip(*matchdaypairs) if day >= start_day_number and day <= end_day_number]
                
                # all valid day numbers in this step
                available_day_numbers = [day for day in day_numbers_in_year(step, step_index) if day >= start_day_number and day <= end_day_number   ]
                
                # indices to available_day_numbers matched by each input file
                index_map_to_days = numpy.searchsorted( available_day_numbers, day_numbers )
                match_with_day = zip(index_map_to_days, matches)

                # construct list of file indices for each day in this step with -1 indicating no input file for days with no match
                matches_this_year = [-1] * len(available_day_numbers)
                
                #print match_with_day
                #print matches_this_year
                
                for day, file_ind in match_with_day:
                    
                    #print day, file_ind
                    matches_this_year[day] = file_ind
                
                step_to_matchgroup[step_index] = matches_this_year
            #print "step_to_matchgroup"
            #print step_to_matchgroup

        else:
            print "STEP NOT UNIQUELY DEFINED BY PATTERN"
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
            #print matchname_to_index
            # And find the operations that match this
            step_to_match = [ matchname_to_index.get(step.create_output_entry(subset.layout.patterns, step_index).name, None) for step_index in range(step.count()) ]
            #print "step_to_match", step_to_match
            step_to_matchgroup = {}
            
            start_day_number = eustace.timeutils.epoch.days_since_epoch(step.start)
            end_day_number = eustace.timeutils.epoch.days_since_epoch(step.end)
            
            #print start_day_number, end_day_number
            
            for step_index, match_index in enumerate(step_to_match):
                if match_index is not None:
                    # all valid day numbers in this step
                    available_day_numbers = [day for day in day_numbers_in_year(step, step_index) if day >= start_day_number and day <= end_day_number   ]
                    
                    # map the file match_index to every day
                    step_to_matchgroup[step_index] = [match_index for day in available_day_numbers] 
            #print "step_to_matchgroup", step_to_matchgroup
                    
        for step_index, data_reference_list in enumerate(operation_data_references):
            # Use the lookup to create references
            if self.missing_data == InputFileList.MISSING_DATA_ALLOWED:
                # indicator of -1 will later indicate that the file was missing
                available_day_numbers = [day for day in day_numbers_in_year(step, step_index) if day >= start_day_number and day <= end_day_number   ]
                default_entryindex = [-1] * len(available_day_numbers)
            else:
                # empty entryindex indicates that there is no matching file and will be caught as an operation with missing data
                default_entryindex = []
            
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





"""

Generate lists of ouptut file names for operation defined by a time batch job.

Work around the single output file requirement of the standard EUMOPPS Output 
and step objects.

"""

class OutputFileList(Output):
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
            operation_references = [ ]
            # Make new entry
            #newentry = step.create_output_entry(subset.layout.patterns, operation_index)
            newentries = self.create_output_entries(step, subset.layout.patterns, operation_index)
            # Check if it already existed
            for newentry in newentries:
                try:

                    entryindex = existing_name_index[newentry.name]

                except KeyError:

                    # Wasn't there - make a new match
                    entryindex = len(subset.matches)
                    subset.matches.append(newentry)
                    existing_name_index[newentry.name] = entryindex

                # Append reference
                operation_references = operation_references + [ datasetindex, self.subsetindex, entryindex ]
             
            references.append(operation_references)
        print references
        # And return corresponding reference list
        return OperationFileListReference(references)

    def create_output_entries(self, step, patterns, operation_index):
        # Mimics behavious of step.create_output_entry for multiple outputs as defined by the placeholder
        
        raise NotImplementedError()

class AnnualBatchDaysOutput(OutputFileList):
    """
    
    Defines daily output files replacing some functionality of the EUMOPPs step object in order to do so.
    
    """
        
    def create_output_entries(self, step, patterns, operation_index):
        # Mimics behavious of step.create_output_entry for multiple outputs as defined by the placeholder
        
        timebase = eumopps.timeutils.timebase.TimeBaseDays(eustace.timeutils.epoch.EPOCH)
        
        if not self.substep_is_uniquely_defined_by(patterns):
            raise RuntimeError('Patterns do not uniquely define output substeps')

        start_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.start) ) )
        end_day_number = int( numpy.floor( eustace.timeutils.epoch.days_since_epoch(step.end) ) )
        
        day_numbers = day_numbers_in_year(step, operation_index)
        day_numbers = [day_numbers[i] for i in range(len(day_numbers)) if day_numbers[i] >= start_day_number and day_numbers[i] <= end_day_number ]
        print day_numbers
        # convert day numbers into datetime objects
        datetimes = [timebase.number_to_datetime(daynumber) for daynumber in day_numbers]
        
        # get list of filenames for each substep day within the operation
        filenames = [self.filename_from_patterns(patterns, t) for t in datetimes]
        print filenames
        # Build the entries
        cataloguefileentries = [CatalogueFileEntry(filename, t) for (filename, t) in zip(filenames, datetimes)]
        print cataloguefileentries
        # Build the entry
        return cataloguefileentries
        
    def substep_is_uniquely_defined_by(self, patterns):
        """Return True if the patterns are sufficient to uniquely identify the step they correspond to, False otherwise"""

        concatpattern = ''.join(patterns)
        return (('%Y' in concatpattern) and ('%m' in concatpattern) and ('%d' in concatpattern))
        
    def filename_from_patterns(self, patterns, t):

        paths = [ datetime_numeric.build_from_pattern(pattern, t) for pattern in patterns ]
        name = os.path.join(*paths)
        return name