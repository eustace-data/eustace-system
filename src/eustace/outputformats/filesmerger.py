"""Contains a class that implements the merging procedure"""

import netCDF4
import os
import os.path
import numpy
from eustace.outputformats.outputvariable import OutputVariable
from eustace.surfaceairmodel.fileio.consistent import ConsistentModelOutputNetCDF
from eustace.surfaceairmodel.fileio.consistent import DatasetAttributesConsistentModelOutput
import copy
import warnings

def clean_list_of_outputs_to_merge(outputs_main, outputs_ancillary,number_of_sources,list_of_indexes):
    """Identify missing outputs and their corresponding position in the list of outputs, remove their correlation indexes, update number of sources"""

    # 1st sanity check
    if len(outputs_main) != len(outputs_ancillary):
	raise ValueError('Number of primary and ancillary files must be the same')

    original_number_of_outputs = len(outputs_main)

    indexes_to_remove = []
    for index, elements in enumerate(zip(outputs_main,outputs_ancillary)):
        print(elements)
        if elements[0] is not None and os.path.isfile(elements[0]):
	    number_of_sources += 1
	else:
	    outputs_main.remove(elements[0])
	    indexes_to_remove.append(index)

	if elements[1] is None or not(os.path.isfile(elements[1])):
	    outputs_ancillary.remove(elements[1])

    print(outputs_main, outputs_ancillary)

    # 2nd sanity check
    if len(outputs_main) != len(outputs_ancillary):
	raise ValueError('Number of primary and ancillary files must be the same')

    # Labeling elements to be removed as None
    for index_to_remove in indexes_to_remove:
      for sub_list in list_of_indexes:
	sub_list_index=0
	for index in range(0,original_number_of_outputs):
	  for jndex in range(index+1,original_number_of_outputs):
	    if index_to_remove in (index,jndex):
	      sub_list[sub_list_index]=None
	    sub_list_index+=1

    for sub_list in list_of_indexes:
	while None in sub_list:
	     sub_list.remove(None)

class FilesMerger(object):
  """Merges output files, produced by multiple sources, into a unique netCDF file, for both fundamental and ancillary fields."""

  def __init__(self,outputs_primary,outputs_ancillary,list_of_daily_sources,primary_fields,primary_uncertainties,ancillary_fields,list_of_correlation_indexes):
    # Basic sanity checks: make it easier to debug when problems occur

    if not isinstance(outputs_primary,list) or not isinstance(outputs_ancillary,list) or not isinstance(list_of_daily_sources,numpy.ndarray) or\
       not isinstance(primary_fields,list) or not isinstance(primary_uncertainties,list) or not isinstance(ancillary_fields,list) or not isinstance(list_of_correlation_indexes,list):
      raise ValueError("Class constructor requires list (outputs_primary,outputs_ancillary,list_of_correlation_indexes,primary_fields,primary_uncertainties,ancillary_fields) and numpy.ndarray (list_of_daily_sources) objects.")

    total = list_of_daily_sources.sum()
    if (len(outputs_primary) != len(outputs_ancillary)) or (len(outputs_ancillary) != total) or ((len(outputs_primary) !=total)):
      raise ValueError("Number of fundamental/ancillary files and total sources has to be the same")

    for name_1,name_2 in zip(outputs_primary,outputs_ancillary):
      if (not isinstance(name_1,unicode) and not isinstance(name_1,str))or (not isinstance(name_2,unicode) and not isinstance(name_2,str)):
	raise ValueError("Fundamental and ancillary files have to be string or unicode types")
    for number in list_of_daily_sources:
      if not isinstance(number,numpy.int64):
	raise ValueError("Number of sources per day has to be numpy.int64 type")
 
    for name in primary_fields+primary_uncertainties+ancillary_fields:
      if not isinstance(name,OutputVariable):
	raise ValueError("Primary_fields, primary_uncertainties and ancillary_field must contain Outputvariable objects")

    if len(ancillary_fields) != len(list_of_correlation_indexes):
      raise ValueError("Number of correlation indexes must macth number of ancillary fields")

    for sub_list in list_of_correlation_indexes:
      if not isinstance(sub_list,list):
	raise ValueError("Correlation indexes list must contain sublists of floats")
      for number in sub_list:
	if not isinstance(number,float):
	  raise ValueError("Correlation indexes sublists must contain floats")

  
    self.outputs_primary = outputs_primary
    self.outputs_ancillary = outputs_ancillary
    self.list_of_daily_sources = list_of_daily_sources
    self.primary_fields = primary_fields
    self.primary_uncertainties = primary_uncertainties
    self.ancillary_fields = ancillary_fields
    self.list_of_correlation_indexes = list_of_correlation_indexes
    self.comment_prefix = 'input data preprocessed from '
    self.list_of_merged_main_outputs = []
    self.list_of_merged_ancillary_outputs = []

  def merge_outputs_from_multiple_daily_sources(self,modulename, institution, catalogue_id):
    """Merges and write output files in final format for a certain number of days"""

    start = 0
    for sources in self.list_of_daily_sources:
	self.merge_outputs_from_single_day(modulename, institution, catalogue_id,self.outputs_primary[start:start+sources],self.outputs_ancillary[start:start+sources])
	start+=sources

  def merge_outputs_from_single_day(self,modulename, institution, catalogue_id,outputs_primary,outputs_ancillary):
    """Merges and write output files in final format for a single day"""
 
   # 1st case: single source
    if len(outputs_primary) == 1:
      self.update_single_file(outputs_primary,outputs_ancillary)
    # 2nd case: multiple sources
    else:
      # collecting the results of the merging
      results=self.merge_multiple_files(outputs_primary,outputs_ancillary)

      # creating list of ancillary variables fields
      ancillary_variables = copy.copy(self.ancillary_fields)

      for output_variable in self.primary_fields+self.primary_uncertainties+self.ancillary_fields:
	for filename in outputs_primary:
	  source_label = os.path.basename(filename).split('_')[4]
	  new_output_variable = copy.copy(output_variable)
	  new_output_variable.name = new_output_variable.name+'_'+source_label
	  new_output_variable.long_name = new_output_variable.long_name+' from '+source_label+' satellite.'
	  ancillary_variables.append(new_output_variable)

      # Output files names
      source_label = os.path.basename(outputs_primary[0]).split('_')[4]
      output_main = outputs_primary[0].replace('_'+source_label,'')
      output_ancillary = outputs_ancillary[0].replace('_'+source_label,'')

      self.list_of_merged_main_outputs.append(output_main)
      self.list_of_merged_ancillary_outputs.append(output_ancillary)

      # Writing information about original sources
      sources = [os.path.basename(filename).split('_')[4] for filename in outputs_primary]
      sources_string = ''
      for name in sources:
	sources_string += name+', '
      sources_string=sources_string[:-2]

      # Print to output files
      attributes = DatasetAttributesConsistentModelOutput(modulename, institution, catalogue_id)
      attributes.comment = self.comment_prefix+sources_string+' retrievals.'
      writer = ConsistentModelOutputNetCDF(attributes,ancillary_variables)
      writer.write_primary(output_main, results)
      writer.write_ancillary(output_ancillary, results)
   
      #deleting old files
      for filename in outputs_primary+outputs_ancillary:
	os.remove(filename) 

  def update_single_file(self,outputs_primary,outputs_ancillary):
    """Updates output filename to adhere SATSTACE convention and adds information about original satellite source into the 'comment' attribute. It is called when just only one output file is available"""  

    source_label = os.path.basename(outputs_primary[0]).split('_')[4]	

    for i in outputs_primary:
      datafile = netCDF4.Dataset(i,'r+')
      datafile.setncattr('comment',self.comment_prefix+source_label+' retrievals.')
      datafile.close()
      new_name = i.replace('_'+source_label,'')
      self.list_of_merged_main_outputs.append(new_name)
      os.rename(i,new_name)

    for i in outputs_ancillary:
      datafile = netCDF4.Dataset(i,'r+')
      datafile.setncattr('comment',self.comment_prefix+source_label+' retrievals.')
      datafile.close()
      new_name = i.replace('_'+source_label,'')
      self.list_of_merged_ancillary_outputs.append(new_name)
      os.rename(i,new_name)


  def merge_multiple_files(self,outputs_primary,outputs_ancillary):
    """Merges multiple output sources into a unique output. Merged fundamental fields are averaged grid-box-by-grid-box. 
       Fundamental and ancillary uncertainties are propagated grid-box-by-grid-box by taking into account of proper correlation among the different output sources.
       Fundamental and ancillary fields from original output sources are appended into the final ancillary file."""
        
    primary_field_names = [ field.name for field in self.primary_fields]
    primary_uncertainties_names = [ field.name for field in self.primary_uncertainties]
    ancillary_field_names = [ field.name for field in self.ancillary_fields]
    results = {}

    # opening all primary files
    handles_primary = []
    for name in outputs_primary:
      handles_primary.append(netCDF4.Dataset(name,'r'))
    
     # opening all ancillary files
    handles_ancillary = []
    for name in outputs_ancillary:
      handles_ancillary.append(netCDF4.Dataset(name,'r'))
    
    # checking that all files have the same date
    times = []
    for j in handles_primary+handles_ancillary:
      times.append(j.variables['time'][...][0])
    
    times=set(times)
    if len(times) == 1:
      results['daynumber'] = times.pop()
    else:
      raise ValueError("Multiple sources have mismatching date.")

    # primary uncertainties are uncorrelated among the instruments
    list_of_primary_correlation_indexes=[[0. for index in range(len(outputs_primary)*(len(outputs_primary)-1)//2)] for jndex in range(len(primary_uncertainties_names))]
    
    # averaging, propagating, appending original results with short and long names changed
    self.append_original_fields(primary_field_names+primary_uncertainties_names,handles_primary,results)
    self.append_original_fields(ancillary_field_names,handles_ancillary,results)
    self.compute_mean_fields(primary_field_names,handles_primary,results)
    self.compute_propagated_uncertainties(primary_uncertainties_names,handles_primary,list_of_primary_correlation_indexes,results)
    self.compute_propagated_uncertainties(ancillary_field_names,handles_ancillary,self.list_of_correlation_indexes,results)
   
    for handle in handles_primary+handles_ancillary:
      handle.close()

    return results

  def check_produced_output(self,expected_output_main,expected_output_ancillary):
    """Sanity check on the names of produced output"""

    check_main = [False for element in self.list_of_merged_main_outputs if element not in expected_output_main]
    check_ancillary = [False for element in self.list_of_merged_ancillary_outputs if element not in expected_output_ancillary]
    
    if False in check_main:
      message = 'Naming of merged {0} field files went wrong: expected output names = {1} \n obtained output names = {2}'.format('fundamental',expected_output_main,self.list_of_merged_main_outputs)
      raise ValueError(message)

    elif False in check_ancillary:
      message = 'Naming of merged {0} field files went wrong: expected output names = {1} \n obtained output names = {2}'.format('fundamental',expected_output_ancillary,self.list_of_merged_ancillary_outputs)
      raise ValueError(message)


  @classmethod
  def compute_mean_fields(cls,names,handles,results):
    """Compute grid-box-by-grid-box mean of primary fields for a given number of input files"""

    for name in names:
  
    # looping through variables names
      old_fill_value = handles[0].variables[name][...].fill_value
      new_fill_value = numpy.nan
      final_values,final_mask = cls.extract_and_modify_masked_array_data_and_mask(handles[0],name,new_fill_value)
      
      for i in range(1,len(handles)): 

	# looping through files
	temporary_values,temporary_mask = cls.extract_and_modify_masked_array_data_and_mask(handles[i],name,new_fill_value)
	final_values = numpy.row_stack((final_values,temporary_values))
	final_mask = numpy.logical_and(final_mask,temporary_mask)

      #averaging: mean of empty slice could occur
      with warnings.catch_warnings():
	warnings.simplefilter("ignore", category=RuntimeWarning)
	mean_field_data = numpy.nanmean(final_values,axis = 0,keepdims = True)
      
      #converting numpy.nans to original filling value
      mean_field_data[numpy.isnan(mean_field_data)] = old_fill_value

      results[name] = numpy.ma.MaskedArray(mean_field_data,mask = final_mask,fill_value = old_fill_value)
    
  @classmethod
  def compute_propagated_uncertainties(cls,names,handles,list_of_correlation_indexes,results):
     """Grid-boxby-grid-box propagation of uncertainties for a given number of input files"""
      
     number_of_files = len(handles)
     for index,name in enumerate(names):
      # looping through variables names

      # building correlation matrix
      sub_list_of_correlation_indexes = list_of_correlation_indexes[index]
      triangular_correlation_matrix = cls.build_triangular_correlation_matrix(number_of_files,sub_list_of_correlation_indexes)
      
      old_fill_value = handles[0].variables[name][...].fill_value
      new_fill_value = 0
      final_values,final_mask = cls.extract_and_modify_masked_array_data_and_mask(handles[0],name,new_fill_value)
      stack_of_masks = ~final_mask # for computing the correct weight factor for each grid-box
      
      for i in range(1,len(handles)): 
	
	# looping through files
	temporary_values,temporary_mask = cls.extract_and_modify_masked_array_data_and_mask(handles[i],name,new_fill_value)
	final_values = numpy.row_stack((final_values,temporary_values))
	stack_of_masks = numpy.row_stack((stack_of_masks,~temporary_mask))
	final_mask = numpy.logical_and(final_mask,temporary_mask)

      # propagating
      propagated_uncertainty = cls.propagate_uncertainty(final_values,triangular_correlation_matrix)
      
      # updating with weights: division by zero could sometimes occur if alla files for a given grid-box do not have data
      with warnings.catch_warnings():
	warnings.simplefilter("ignore", category=RuntimeWarning)
	weights = 1./stack_of_masks.sum(axis=0,keepdims=True) 
	propagated_uncertainty *= weights

      # converting potential numpy.nans (due to division by zero) to original filling values 
      propagated_uncertainty[numpy.isnan(propagated_uncertainty)] = old_fill_value

      results[name] = numpy.ma.MaskedArray(propagated_uncertainty,mask = final_mask,fill_value = old_fill_value)	
 
  @staticmethod
  def extract_and_modify_masked_array_data_and_mask(handle,name,new_fill_value):
    """Extracts masked array data and mask, sets original filling value to a new one"""

    masked_array_field = handle.variables[name][...]
    masked_array_field.data[masked_array_field.data==masked_array_field.fill_value]=new_fill_value

    return masked_array_field.data,masked_array_field.mask

  @staticmethod
  def build_triangular_correlation_matrix(number_of_files,sub_list_of_correlation_indexes):
    """Builds upper triangular matrix containing elements of covariance matrix: upper diagonal elements are multiplied by a factor of 2"""

    if (number_of_files*(number_of_files-1.))/2. != len(sub_list_of_correlation_indexes):
      raise ValueError("Number of correlation indexes [r(file_1,file_2),r(file_1,file_3),...,r(file_{n-1},file_{n})] must be equal to [(number_of_files*(number_of_files-1.))/2].")

    matrix=numpy.zeros([number_of_files,number_of_files])
    numpy.fill_diagonal(matrix, 1)
    sub_list_index = 0
    
    # Run just over upper triangle, discarding the diagonal
    for index in range(0,number_of_files):
      for jndex in range(index+1,number_of_files):
	matrix[index,jndex] = 2.*sub_list_of_correlation_indexes[sub_list_index]
	sub_list_index+=1

    return matrix

  @staticmethod
  def propagate_uncertainty(stack_of_uncertainties,triangular_correlation_matrix):
    """Propagate single uncertainty along the ''files direction'', according to a given covariance matrix"""

    result=0.
    size=triangular_correlation_matrix.shape[0]
    
    # Run just over upper triangle, taking the diagonal
    for i in range(size):
      for j in range(i,size):	
	result += stack_of_uncertainties[i]*triangular_correlation_matrix[i,j]*stack_of_uncertainties[j]

    return numpy.sqrt(numpy.array([result]))

  @staticmethod
  def append_original_fields(names, handles, results):
    """Appends original fundamental and ancillary fields to the dictionary containing merged data"""
    for handle in handles:
      for name in names:
	data=handle.variables[name][...]
	filename=os.path.basename(handle.filepath())
	source = filename.split('_')[4]
	newname = name+"_"+source
	results[newname] = data
