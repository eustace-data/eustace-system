"""Contains correlation length/time scales definitions, and a method for converting them into attributes of the correlated uncertainties field"""

from eustace.outputformats.outputvariable import OutputVariable

class CorrelationRanges(object):
  """Base class for storing information about correlation ranges and feeding them into attributes of the correlated uncertainties field"""
  def __init__(self,keys,units,length_values,time_values):
      if not isinstance(keys,list) or not isinstance(units,list) or not isinstance(length_values,list) or not isinstance(time_values,list):
	raise ValueError('Class attributes must be defined in terms of list of values')
      else:
	self.keys=keys
	self.units=units
	self.length_values=length_values
	self.time_values=time_values
      
  def update_correlated_uncertainty_ranges(self,outputvariable):
    """Takes an output variable object as input and update the information about lenght/time scales and their relative units"""
    key=outputvariable.name
    ranges_dictionary=self.length_time_scales_dictionary()[key]
    outputvariable.length_scale=ranges_dictionary['length scale']
    outputvariable.length_scale_units=ranges_dictionary['length scale units']
    outputvariable.time_scale=ranges_dictionary['time scale']
    outputvariable.time_scale_units=ranges_dictionary['time scale units']
    
  def length_time_scales_dictionary(self):
    """Returns a dictionary containing information about lenght/time scales and their relative units"""
    return {uncertainty_key:{'length scale':length,'time scale':time, 'length scale units':self.units[0], 'time scale units':self.units[1]} for (uncertainty_key,length,time) in zip(self.keys,self.length_values,self.time_values) }
