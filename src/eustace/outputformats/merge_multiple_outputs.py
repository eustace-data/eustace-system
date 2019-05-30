"""Merges outputs obtained from different sources into a unique set of primary and ancillary fields netCDF file"""

import ast
import importlib
import numpy
import os
from eustace.outputformats.filesmerger import clean_list_of_outputs_to_merge
from eustace.outputformats.filesmerger import FilesMerger

def return_daily_outputs(catalogue_id,institution, 
                         output_main, output_ancillary, 
                         correlation_indexes, comment_prefix, definitions_module,
                         primary_fields_list, uncertainty_fields_list, ancillary_fields_list,
                         correlation_ranges,
                         input_main, input_ancillary):
    """Collect single-day relevant information into a nested-keys dictionary"""

    return {'catalogue_id':catalogue_id, 'institution':institution, 
            'output_main':output_main, 'output_ancillary':output_ancillary, 
            'correlation_indexes':correlation_indexes, 'comment_prefix':comment_prefix, 'definitions_module':definitions_module,
            'primary_fields_list':primary_fields_list, 'uncertainty_fields_list':uncertainty_fields_list, 'ancillary_fields_list':ancillary_fields_list,
            'correlation_ranges':correlation_ranges,
            'input_main':input_main, 'input_ancillary':input_ancillary}


def merge_single_day_outputs(catalogue_id, institution, 
			     output_main, output_ancillary, 
			     correlation_indexes, comment_prefix, definitions_module,
			     primary_fields_list, uncertainty_fields_list, ancillary_fields_list,
			     correlation_ranges,
			     input_main, input_ancillary):
    """Merge output files from different sources into a unique product, for a specific processed day"""

    # import definitions module
    module = importlib.import_module(definitions_module)
    
    # import primary, uncertainty, ancillary fields objects
    primary_fields = [getattr(module, name) for name in primary_fields_list]
    primary_uncertainties = [getattr(module, name) for name in uncertainty_fields_list]
    ancillary_fields = [getattr(module, name) for name in ancillary_fields_list]

    #update information about correlation ranges    
    for variable in ancillary_fields:
      if variable.name in correlation_ranges.keys:
	correlation_ranges.update_correlated_uncertainty_ranges(variable)

    # correlation indexes for ancillary_fields along (all possibly available) files direction
    list_of_correlation_indexes = ast.literal_eval(correlation_indexes)
    number_of_daily_sources = numpy.zeros([1], dtype=numpy.int64)
          
    
    # removing None values (referring to missing files) from input list, readapting the list of correlation indexes
    clean_list_of_outputs_to_merge(input_main, input_ancillary, number_of_daily_sources, list_of_correlation_indexes)

    if number_of_daily_sources.sum():
      # if we have some output to merge, then merge it
      merger=FilesMerger(input_main, input_ancillary, number_of_daily_sources, primary_fields, primary_uncertainties, ancillary_fields, list_of_correlation_indexes)
      merger.comment_prefix = comment_prefix
      merger.merge_outputs_from_multiple_daily_sources(__name__, institution,catalogue_id)
      merger.check_produced_output(output_main,output_ancillary)
    else:
      print("eustace.outputformats.merge_single_day_outputs WARNING:No input files to merge")

def merge_multiple_days_outputs(list_of_daily_outputs):
    """Merge a output files from different sources into a unique product, for a certain number of processed days"""

    for dictionary in list_of_daily_outputs:
	merge_single_day_outputs(**dictionary)