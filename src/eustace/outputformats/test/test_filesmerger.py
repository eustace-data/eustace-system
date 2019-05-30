"""Test of multiple files merging procedure"""

import json
import numpy
import numpy.ma
import numpy.testing 
import os
import tempfile
import unittest
from ..filesmerger import clean_list_of_outputs_to_merge
from ..filesmerger import FilesMerger
from datetime import datetime
from eumopps.timeutils.timebase import TimeBaseDays
from eustace.outputformats.outputvariable import OutputVariableTemplate
from eustace.outputformats.outputvariable import OutputVariable
from eustace.outputformats.definitions import TAS, TASUNCERTAINTY
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_RANDOM
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED2
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC2
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER0
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER1
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER2
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER3
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER4
from netCDF4 import Dataset
from shutil import copyfile

# Output variable template
t = OutputVariableTemplate(numpy.float32, -3200, units='mm')

# Output variable objects created fro testing purposes
f1 = OutputVariable.from_template(t, 'field_1', 'first field', ancillary_variables='Bob\'s Shed')
f2 = OutputVariable.from_template(t, 'field_2', 'second field', ancillary_variables='Bob\'s Shed')
u1 = OutputVariable.from_template(t, 'total_uncertainty_field_1', 'first total uncertainty field')
u2 = OutputVariable.from_template(t, 'total_uncertainty_field_2', 'second total uncertainty field')
a11 = OutputVariable.from_template(t, 'ancillary_1_for_field_1', 'first ancillary variable for first field')
a12 = OutputVariable.from_template(t, 'ancillary_1_for_field_2', 'second ancillary variable for first field')
a21 = OutputVariable.from_template(t, 'ancillary_2_for_field_1', 'first ancillary variable for second field')
a22 = OutputVariable.from_template(t, 'ancillary_2_for_field_2', 'second ancillary variable for second field')

# LIST OF CORRELATION INDEXES

list_of_correlation_indexes = [[1.,0.,.5],[0.,.25,1.],[.1,0.2,.5],[.03,.5,.3]]
#   a11        a12          a21         a22

# CONSTANT ARRAYS USED FOR TESTING PURPOSES
  
TEST_MASKS = [[numpy.array([[[False,   True, False],[  True,False,False],[False,False,False],[False,False, True]]]),
               numpy.array([[[False,False,False],[False,False,False],[False,False,False],[  True,  True,False]]])],

              [numpy.array([[[False,  True, False],[False,  True, False],[  True,False,False],[  True,  True, True]]]),
               numpy.array([[[False,False,True],[False,False,False],[  True,False,False],[False,False, True]]])],

              [numpy.array([[[False,  True,False],[False,  True,False],[  True,False,False],[  True,  True, True]]]),
               numpy.array([[[False,False, True],[False,False,False],[  True,False,False],[ False,False, True]]])]]

TEST_FUNDAMENTAL_FIELDS = [[numpy.array([[[    3.,   -111.,    0.],[  -111.,     0.,  -10.],[    4.,     18.,    4.],[   -7.,     -7., -111.]]]),
          	            numpy.array([[[     3.,     0.,    0.],[    77.,    30.,   22.],[   -11.,   -22.,   66.],[  -111.,  -111.,   -1.]]])],

                           [numpy.array([[[     3.,  -111.,    0.],[     0.,  -111.,   22.],[  -111.,    18.,   6.],[  -111.,  -111., -111.]]]),
                            numpy.array([[[     0.,    19., -111.],[    20.,   -10.,   30.],[  -111.,    18.,   66.],[     0.,     8., -111.]]])],

                           [numpy.array([[[     3.,  -111.,    0.],[     0.,  -111.,   22.],[  -111.,    18.,   6.],[  -111.,  -111., -111.]]]),
                            numpy.array([[[     0.,    19., -111.],[    20.,   -10.,   30.],[  -111.,    18.,   66.],[     0.,     8., -111.]]])]]

TEST_FUNDAMENTAL_UNCERTAINTIES = [[numpy.array([[[  .3,  -111.,     0.1],[-111.,     0.9,  .1],[    .4,   .18,   .4],[     .7,      7., -111.]]]),
                                   numpy.array([[[  .3,      0.2,     0.2],[    .7,    .3, .22],[   .11,   .22,  .66],[ -111.,  -111.,     1.]]])],

                                  [numpy.array([[[  .3,  -111.,     0.3],[    0.9, -111., .22],[-111.,   .18,  .66],[ -111.,  -111., -111.]]]),
                                   numpy.array([[[  0.4,     .19, -111.],[    .2,     .1,  .3],[-111.,   .18,  .66],[     0.,      .8, -111.]]])],

                                  [numpy.array([[[ .3,  -111.,     0.2],[    0.9, -111., .22],[-111.,   .18,  .16],[ -111.,  -111., -111.]]]),
                                   numpy.array([[[  0.5,     .29, -111.],[    .2,     .1,  .3],[-111.,   .18,  .46],[     0.,      .8, -111.]]])]]

TEST_ANCILLARY_11 = [numpy.array([[[    3.,   -111.,    0.],[  -111.,     0.,  +10.],[    4.,     18.,    4.],[   +7.,     +7., -111.]]]),
                     numpy.array([[[     3.,  -111.,    0.],[     0.,  -111.,   22.],[  -111.,    18.,   66.],[  -111.,  -111., -111.]]]),
                     numpy.array([[[     3.,  -111.,    0.],[     0.,  -111.,   22.],[  -111.,    18.,   66.],[  -111.,  -111., -111.]]])]

TEST_ANCILLARY_12 = [numpy.array([[[     3.,     0.,    .24],[    77.,    .13,   22.],[   +11.,   +22.,   66.],[  -111.,  -111.,   +1.]]]),
                     numpy.array([[[     0.,    19., -111.],[    .2,   +10.,   .13],[  -111.,    18.,   66.],[     .24,     8., -111.]]]),
                     numpy.array([[[     0.,    19., -111.],[    .2,   +10.,   .13],[  -111.,    18.,   66.],[     .24,     8., -111.]]])]

TEST_ANCILLARY_21 = [numpy.array([[[    3.,   -111.,    0.],[  -111.,     0.,  +10.],[    4.,     18.,    4.],[   +7.,     +7., -111.]]]),
                     numpy.array([[[     3.,  -111.,    0.],[     0.,  -111.,   22.],[  -111.,    18.,   .66],[  -111.,  -111., -111.]]]),
                     numpy.array([[[     3.,  -111.,    0.],[     0.,  -111.,   22.],[  -111.,    18.,   .66],[  -111.,  -111., -111.]]])]

TEST_ANCILLARY_22 = [numpy.array([[[     3.,     0.,    0.],[    77.,    .13,   22.],[   +11.,   +22.,   66.],[  -111.,  -111.,   +1.]]]),
                     numpy.array([[[     0.,    19., -111.],[    .2,   +10.,   .13],[  -111.,    .18,   66.],[     0.,     8., -111.]]]),
                     numpy.array([[[     0.,    19., -111.],[    .2,   +10.,   .13],[  -111.,    .18,   66.],[     0.,     8., -111.]]])]

    
EXPECTED_MERGED_FIELDS = { 'field_1' : numpy.array([[[3., -111., 0. ],[0.,0.,34./3.],[4., 18., 16./3.],[-7.,-7.,-111.]]]) , 
			    'field_2' : numpy.array([[[1., 38./3., 0.],[ 117./3., 10./3., 82./3.],[-11., 14./3., 66.],[ 0., 8., -1.]]]),
			    'total_uncertainty_field_1' : numpy.array([[[ 0.17320508, -111.,0.12472191],[ 0.6363961 ,0.9, 0.10893423],[ 0.4,0.10392305,0.26272081],[0.7,7.,-111.]]]),
			    'total_uncertainty_field_2' : numpy.array([[[ 0.23570226,  0.13341664,  0.2 ],[ 0.25166115,  0.11055416,  0.15930404],[ 0.11 , 0.11215069,  0.34685892],[ 0.,0.56568542,1.]]]), 
			    'ancillary_1_for_field_1' : numpy.array([[[  2.44948974, -111.,0.], [0.,0.,  14.87727574],[4.,14.69693846,38.89015871],[7.,7.,-111.]]]),
			    'ancillary_1_for_field_2' : numpy.array([[[  1.,12.66666667,0.24 ],[ 25.683676,   6.67222185,   7.34467003],[ 11.,  14.82490397,  51.59457336],[0.24, 8., 1.]]]),
			    'ancillary_2_for_field_1' : numpy.array([[[  2.14476106, -111.,0.],[0.,0.,13.67885635],[  4.,  12.86856635,   1.44878493],[7.,7., -111.]]]),
			    'ancillary_2_for_field_2' : numpy.array([[[  1.,10.21219315,0.],[ 25.7022,   5.38923722,   7.35659598],[ 11,   7.36570009,  47.49147],[  0.,6.4498062 ,1.]]])}

EXPECTED_MERGED_FIELDS_MASKS = { 'field_1' : numpy.array([[[False, True,  False],[False,False,False],[False,False,False],[False,False,True]]]) , 
				  'field_2' : numpy.array([[[False, False, False],[False,False,False],[False,False,False],[False,False,False]]]),
				  'total_uncertainty_field_1' : numpy.array([[[False, True,  False],[False,False,False],[False,False,False],[False,False,True]]]) , 
				  'total_uncertainty_field_2' : numpy.array([[[False, False, False],[False,False,False],[False,False,False],[False,False,False]]]),
				  'ancillary_1_for_field_1' : numpy.array([[[False, True,  False],[False,False,False],[False,False,False],[False,False,True]]]) , 
				  'ancillary_1_for_field_2' : numpy.array([[[False, False, False],[False,False,False],[False,False,False],[False,False,False]]]),
				  'ancillary_2_for_field_1' : numpy.array([[[False, True,  False],[False,False,False],[False,False,False],[False,False,True]]]) , 
				  'ancillary_2_for_field_2' : numpy.array([[[False, False, False],[False,False,False],[False,False,False],[False,False,False]]])} 

APPENDED_TEST_FIELDS = [{ 'field_1_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_FUNDAMENTAL_FIELDS[index][0], mask = TEST_MASKS[index][0], fill_value = -111),
			  'field_2_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_FUNDAMENTAL_FIELDS[index][1],mask = TEST_MASKS[index][1], fill_value = -111),
			  'total_uncertainty_field_1_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_FUNDAMENTAL_UNCERTAINTIES[index][0], mask = TEST_MASKS[index][0], fill_value = -111),
			  'total_uncertainty_field_2_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_FUNDAMENTAL_UNCERTAINTIES[index][1], mask = TEST_MASKS[index][1], fill_value = -111),
			  'ancillary_1_for_field_1_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_ANCILLARY_11[index], mask = TEST_MASKS[index][0], fill_value = -111),
			  'ancillary_1_for_field_2_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_ANCILLARY_12[index], mask = TEST_MASKS[index][1], fill_value = -111),
			  'ancillary_2_for_field_1_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_ANCILLARY_21[index], mask = TEST_MASKS[index][0], fill_value = -111),
			  'ancillary_2_for_field_2_LABEL'+str(index) : numpy.ma.MaskedArray(TEST_ANCILLARY_22[index], mask = TEST_MASKS[index][1], fill_value = -111),} for index in range(3)]


# Testing

class TestFilesMerger(unittest.TestCase):
  
  def test_clean_list_of_outputs_to_merge(self):

     test_input = {'outputs_main': ['A','B','C'],'outputs_ancillary': ['A','B'],'number_of_sources':None, 'list_of_indexes':None}
     self.assertRaises(ValueError,clean_list_of_outputs_to_merge,**test_input)

     test_input = {'outputs_main': ['A','B',None,'C'],'outputs_ancillary': ['A','B'],'number_of_sources':None, 'list_of_indexes':None}
     self.assertRaises(ValueError,clean_list_of_outputs_to_merge,**test_input)

     test_input = {'outputs_main': ['A','B',None],'outputs_ancillary': ['A','B',None, 'D'],'number_of_sources':numpy.zeros(1,dtype=numpy.int64), 'list_of_indexes':[[1.,.5,-.2],[.2,.3,-.46]]}
     self.assertRaises(ValueError, clean_list_of_outputs_to_merge,**test_input)

     test_input = {'outputs_main': ['A','B'],'outputs_ancillary': ['A','B'],'number_of_sources':numpy.zeros(1,dtype=numpy.int64), 'list_of_indexes':[[1.,-.2],[.2,-.46]]}
     
     self.assertItemsEqual(test_input['outputs_main'],['A','B'])
     self.assertItemsEqual(test_input['outputs_ancillary'],['A','B'])
     
     numpy.testing.assert_array_equal(test_input['number_of_sources'],numpy.array([0]),err_msg='Testing correct sources counting')
     self.assertItemsEqual(test_input['list_of_indexes'],[[1.,-.2],[.2,-.46]])

     test_input = {'outputs_main': [None,'B',None,'C'],'outputs_ancillary': [None,'B',None,'C'],'number_of_sources':numpy.zeros(1,dtype=numpy.int64).tolist(), 'list_of_indexes':[[1.,.8,.6,.2,.3,-.6],[1.,.5,-.2,.4,.1,-.3],[.7,-.8,.9,.2,.3,-.46]]}

     clean_list_of_outputs_to_merge(**test_input)

     self.assertItemsEqual(test_input['outputs_main'],['B','C'])
     self.assertItemsEqual(test_input['outputs_ancillary'],['B','C'])
     numpy.testing.assert_array_equal(test_input['number_of_sources'],numpy.array([2]),err_msg='Testing correct sources counting')
     self.assertItemsEqual(test_input['list_of_indexes'],[[.8,.2],[.5,.4],[-.8,.2]])	
     
     test_input = {'outputs_main': ['A','B',None,'C'],'outputs_ancillary': ['A','B',None,'C'],'number_of_sources':numpy.zeros(1,dtype=numpy.int64), 'list_of_indexes':[[1.,.8,.6,.2,.3,-.6],[1.,.5,-.2,.4,.1,-.3],[.7,-.8,.9,.2,.3,-.46]]}
     clean_list_of_outputs_to_merge(**test_input)

     self.assertItemsEqual(test_input['outputs_main'],['A','B','C'])
     self.assertItemsEqual(test_input['outputs_ancillary'],['A','B','C'])
     numpy.testing.assert_array_equal(test_input['number_of_sources'],numpy.array([3]),err_msg='Testing correct sources counting')
     self.assertItemsEqual(test_input['list_of_indexes'],[[1.,.6,.3],[1.,-.2,.1],[.7,.9,.3]])	


  def  test_init(self):
    #Test the correcteness of the init procedure

    v = OutputVariable.from_template(t, 'bob', 'Bob\'s House', ancillary_variables='Bob\'s Shed')

    self.assertRaises(ValueError,FilesMerger,1,['a'],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],1,numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],['a'],[1],[v],[v],[v],[v])
    self.assertRaises(ValueError,FilesMerger,[1],['a'],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],[1],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],['a'],['a'],[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],['a'],['a'],[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],['a','b'],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],['a','b'],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a','b'],['a','b'],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a'],['a','b'],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a','b'],['a','b'],numpy.array([1]),[v],[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a','b'],['a','b'],numpy.array([2]),1,[v],[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a','b'],['a','b'],numpy.array([2]),[v],1,[v],[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a','b'],['a','b'],numpy.array([2]),[v],[v],1,[[1.]])
    self.assertRaises(ValueError,FilesMerger,['a','b'],['a','b'],numpy.array([2]),[v],[v],[v],1.)
    self.assertRaises(ValueError,FilesMerger,['a','b'],['a','b'],numpy.array([2]),[v],[v],[v],[[True]])
 
    merger=FilesMerger(['a','b','c'],['d','e','f'],numpy.array([1,2]),[v],[v],[v],[[1.]])
    self.assertEqual(merger.outputs_primary,['a','b','c'])
    self.assertEqual(merger.outputs_ancillary,['d','e','f'])
    self.assertEqual(merger.list_of_daily_sources[0],1)
    self.assertEqual(merger.list_of_daily_sources[1],2)
    self.assertEqual(merger.primary_fields,[v])
    self.assertEqual(merger.primary_uncertainties,[v])
    self.assertEqual(merger.ancillary_fields,[v])
    self.assertEqual(merger.list_of_correlation_indexes,[[1.]])
    self.assertEqual(merger.list_of_merged_main_outputs,[])
    self.assertEqual(merger.list_of_merged_ancillary_outputs,[])

  def test_update_single_file(self):
    # Testing the update of single files

    outputs = [ tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.FIELD_SURFACE_eustace_PRODUCTIONUMBER_LABEL_XXXX', suffix='.nc',delete=False) for index in range(2)]
    dataset = [ Dataset(output.name,'w','NETCDF4') for output in outputs]

    for files in dataset:
      files.comment=''
      files.close()
    v = OutputVariable.from_template(t, 'bob', 'Bob\'s House', ancillary_variables='Bob\'s Shed')
    merger=FilesMerger(['a','b','c'],['d','e','f'],numpy.array([1,2]),[v],[v],[v],[[1.]])
    merger.update_single_file([outputs[0].name],[outputs[1].name])

    new_outputs = [ output.name.replace('_LABEL','') for output in outputs]
    dataset = [ Dataset(name,'r','NETCDF4') for name in new_outputs]

    for files in dataset:
	self.assertEqual(files.ncattrs(),['comment'])
	self.assertEqual(files.getncattr(files.ncattrs()[0]),u'input data preprocessed from LABEL retrievals.')
	files.close()
    
    for name in new_outputs:
      os.remove(name)
    

  def test_extract_and_modify_masked_array_data_and_mask(self):
    # Testing correct extraction and modification of masked array data and masks

    output = tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.FIELD_SURFACE_eustace_PRODUCTIONUMBER_LABEL_XXXX', suffix='.nc',delete=True)
    A=numpy.ma.MaskedArray([[-100.,1.,2.],[3.,-100.,-100.],[6.,7.,-100.]],mask=[[True,False,False],[False,True,True],[False,False,True]],fill_value=-100.)
    B=numpy.ma.MaskedArray([[200.,1.,2.],[200.,-10.,200.],[200.,7.,-100.]],mask=[[True,False,False],[True,False,True],[True,False,False]],fill_value=200.)
    C=numpy.ma.MaskedArray([[-1.,1.,2.],[45.,-1.,1.],[6.,-1.,9.]],mask=[[True,False,False],[False,True,False],[False,True,False]],fill_value=-1.)
    
    dataset=Dataset(output.name,'w','NETCDF4')
    dataset.createDimension('x',3)
    dataset.createDimension('y',3)
    field_A=dataset.createVariable('field_A',numpy.float32,dimensions=('x','y'),fill_value=A.fill_value)
    field_B=dataset.createVariable('field_B',numpy.float32,dimensions=('x','y'),fill_value=B.fill_value)
    field_C=dataset.createVariable('field_C',numpy.float32,dimensions=('x','y'),fill_value=C.fill_value)
    field_A[:]=A
    field_B[:]=B
    field_C[:]=C
    dataset.close()

    dataset=Dataset(output.name,'r','NETCDF4')
    variables_names=dataset.variables.keys()
    new_fill_values=[numpy.nan,18.,44.]
    original_arrays=[A,B,C]

    for name,new_fill_value,array in zip(variables_names,new_fill_values,original_arrays):
      test_array_data,test_array_mask = FilesMerger.extract_and_modify_masked_array_data_and_mask(dataset,name,new_fill_value)
      array.data[array.data==array.fill_value]=new_fill_value
      numpy.testing.assert_array_equal(array.data,test_array_data,err_msg="Masked array data values mismatch")
      numpy.testing.assert_array_equal(array.mask,test_array_mask,err_msg="Masked array mask values mismatch")

    dataset.close()

  def test_compute_mean_fields(self):
    # Testing computation of the mean along multiple source files

    outputs = [ tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.FIELD_SURFACE_eustace_PRODUCTIONUMBER_LABEL_XXXX', suffix='.nc',delete=True) for index in range(3)]
    dataset = [ Dataset(output.name,'w','NETCDF4') for output in outputs]
    
    for index,files in enumerate(dataset):
      files.createDimension('t',1)
      files.createDimension('x',4)
      files.createDimension('y',3)
      field_1=files.createVariable('field_1',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      field_2=files.createVariable('field_2',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      field_1.set_auto_maskandscale(True)
      field_2.set_auto_maskandscale(True)
      field_1[:]=TEST_FUNDAMENTAL_FIELDS[index][0]
      field_2[:]=TEST_FUNDAMENTAL_FIELDS[index][1]
      files.close()
   
    handles = []
    for output in outputs:
      handles.append(Dataset(output.name,'r'))
    
    primary_names=['field_1','field_2']
    results={}	

    FilesMerger.compute_mean_fields(primary_names,handles,results)

    for name in primary_names:
      numpy.testing.assert_array_equal(EXPECTED_MERGED_FIELDS_MASKS[name],results[name].mask,err_msg="Testing mask correctness for "+name)
      numpy.testing.assert_array_almost_equal(EXPECTED_MERGED_FIELDS[name],results[name].data,err_msg="Testing averaged values for "+name)

    for files in handles:
      files.close()

  def test_build_triangular_correlation_matrix(self):
    # Testing the generation of upper triangular matrix containing covariances information
    self.assertRaises(ValueError,FilesMerger.build_triangular_correlation_matrix,3,[1,2])

    matrices = {'3' : FilesMerger.build_triangular_correlation_matrix(3,range(1,4)),
                '4' : FilesMerger.build_triangular_correlation_matrix(4,range(1,7)),
                '6' : FilesMerger.build_triangular_correlation_matrix(6,range(1,16))}

    test_matrices = { '3' : numpy.array([[1.,2.,4.],
			                 [0.,1.,6.],
			                 [0.,0.,1.]]),
		      '4' : numpy.array([[1.,2.,4.,6. ],
					 [0.,1.,8.,10.],
					 [0.,0.,1.,12.],
					 [0.,0.,0.,1.]]),
		      '6' : numpy.array([[1.,2., 4.,6., 8.,10.],
					 [0.,1.,12.,14.,16.,18.],
					 [0.,0., 1.,20.,22.,24.],
					 [0.,0., 0., 1.,26.,28.],
					 [0.,0., 0., 0., 1.,30.],
					 [0.,0., 0., 0., 0., 1.]]) }

    for i in ['3','4','6']:
      numpy.testing.assert_array_equal(matrices[i],test_matrices[i],err_msg="Testing generation of triangular matrix for "+i+" source files")

  def test_propagate_uncertainty(self):
    # Testing the propagation of a single uncertainty

    A=numpy.array([[[1.,2.,3.],[4.,5.,6.]]])
    B=numpy.array([[[3.,1.,2.],[1.,1.,6.]]])
    C=numpy.array([[[0.,1.,3.],[5.,9.,6.]]])

    D=numpy.row_stack((numpy.row_stack((A,A)),A))
    E=numpy.row_stack((numpy.row_stack((A,B)),C))
    
    input_data_stack = [D,D,D,E]

    input_matrices = [numpy.array([[1.,0,0],[0,1.,0],[0,0,1.]]),
                     numpy.array([[1.,0.5,0],[0,1.,0.5],[0,0,1.]]),
		     numpy.array([[1.,0.5,-2],[0,1.,1.5],[0.,0.,1.]]),
		     numpy.array([[1.,0.5, -2],[0, 1.,1.5],[0.,0.,1.]])]

    expected_results = [numpy.sqrt(3)*A,2.*A,numpy.sqrt(3)*A,numpy.array([[[  3.39116499, 2.12132034,4. ],[3.39116499,5.74456265,10.39230485]]])]
    error_messages = ["Testing propagation of single uncertainty, equal arrays, diagonal covariance matrix.",
                      "Testing propagation of single uncertainty, equal arrays, non diagonal covariance matrix with positive entries",
                      "Testing propagation of single uncertainty, equal arrays, non diagonal covariance matrix with negative entries",
                      "Testing propagation of single uncertainty, equal arrays, non diagonal covariance matrix with negative entries"]
    for index in range(4):
      numpy.testing.assert_array_almost_equal(FilesMerger.propagate_uncertainty(input_data_stack[index],input_matrices[index]),expected_results[index],err_msg=error_messages[index])
      
  def test_compute_propagated_uncertainties(self):
    # Testing propagation of uncertainties along multiple source files

    outputs = [ tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.FIELD_SURFACE_eustace_PRODUCTIONUMBER_LABEL_XXXX', suffix='.nc',delete=True) for index in range(3)]
    dataset = [ Dataset(output.name,'w','NETCDF4') for output in outputs]
    
    for index,files in enumerate(dataset):
      files.createDimension('t',1)
      files.createDimension('x',4)
      files.createDimension('y',3)
      field_1=files.createVariable('ancillary_1_for_field_1',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      field_2=files.createVariable('ancillary_1_for_field_2',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      field_1.set_auto_maskandscale(True)
      field_2.set_auto_maskandscale(True)
      field_1[:]=TEST_ANCILLARY_11[index]
      field_2[:]=TEST_ANCILLARY_12[index]
      files.close()
  
    handles = []
    for output in outputs:
      handles.append(Dataset(output.name,'r'))
    
    primary_names=['ancillary_1_for_field_1','ancillary_1_for_field_2']
    results={}
    FilesMerger.compute_propagated_uncertainties(primary_names,handles,[list_of_correlation_indexes[0],list_of_correlation_indexes[1]],results) 
    
    for name in primary_names:
      numpy.testing.assert_array_equal(EXPECTED_MERGED_FIELDS_MASKS[name],results[name].mask,err_msg="Testing mask correctness for "+name)
      numpy.testing.assert_array_almost_equal(EXPECTED_MERGED_FIELDS[name],results[name].data,err_msg="Testing propagation of "+name)

  def test_merge_multiple_files(self):
    # Testing propagation of uncertainties along multiple source files
    
    outputs_primary = [ tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.FUNDAMENTALFIELD_SURFACE_eustace_PRODUCTIONUMBER_LABEL'+str(index)+'_XXXX', suffix='.nc',delete=True) for index in range(3)]
    daynumber=34093

    # filling fundamental fields files
    dataset = [ Dataset(output.name,'w','NETCDF4') for output in outputs_primary]
    for index,files in enumerate(dataset):
      files.createDimension('t',1)
      files.createDimension('x',4)
      files.createDimension('y',3)
      time=files.createVariable('time',numpy.float32,dimensions=('t'),fill_value=-111)
      field_1=files.createVariable('field_1',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      field_2=files.createVariable('field_2',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      total_uncertainty_field_1=files.createVariable('total_uncertainty_field_1',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      total_uncertainty_field_2=files.createVariable('total_uncertainty_field_2',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      field_1.set_auto_maskandscale(True)
      field_2.set_auto_maskandscale(True)
      total_uncertainty_field_1.set_auto_maskandscale(True)
      total_uncertainty_field_2.set_auto_maskandscale(True)
      time[:]=numpy.array([daynumber])
      field_1[:]=TEST_FUNDAMENTAL_FIELDS[index][0]
      field_2[:]=TEST_FUNDAMENTAL_FIELDS[index][1]
      total_uncertainty_field_1[:]=TEST_FUNDAMENTAL_UNCERTAINTIES[index][0]
      total_uncertainty_field_2[:]=TEST_FUNDAMENTAL_UNCERTAINTIES[index][1]
      files.close()

    outputs_ancillary = [ tempfile.NamedTemporaryFile(prefix='eustace.outputformats.test.ANCILLARYFIELD_SURFACE_eustace_PRODUCTIONUMBER_LABEL'+str(index)+'_XXXX', suffix='.nc',delete=True) for index in range(3)]
    # filling ancillary fields files
    dataset = [ Dataset(output.name,'w','NETCDF4') for output in outputs_ancillary]
    for index,files in enumerate(dataset):
      files.createDimension('t',1)
      files.createDimension('x',4)
      files.createDimension('y',3)
      time=files.createVariable('time',numpy.float32,dimensions=('t'),fill_value=-111)
      time[:]=numpy.array([daynumber])
      ancillary_1_for_field_1=files.createVariable('ancillary_1_for_field_1',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      ancillary_1_for_field_1.set_auto_maskandscale(True)
      ancillary_1_for_field_1[:]=TEST_ANCILLARY_11[index]
      ancillary_1_for_field_2=files.createVariable('ancillary_1_for_field_2',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      ancillary_1_for_field_2.set_auto_maskandscale(True)
      ancillary_1_for_field_2[:]=TEST_ANCILLARY_12[index]
      ancillary_2_for_field_1=files.createVariable('ancillary_2_for_field_1',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      ancillary_2_for_field_1.set_auto_maskandscale(True)
      ancillary_2_for_field_1[:]=TEST_ANCILLARY_21[index]
      ancillary_2_for_field_2=files.createVariable('ancillary_2_for_field_2',numpy.float32,dimensions=('t','x','y'),fill_value=-111)
      ancillary_2_for_field_2.set_auto_maskandscale(True)
      ancillary_2_for_field_2[:]=TEST_ANCILLARY_22[index]
      files.close()

    merger=FilesMerger([output.name for output in outputs_primary],[output.name for output in outputs_ancillary],numpy.array([3]),[f1,f2],[u1,u2],[a11,a12,a21,a22],list_of_correlation_indexes)   
    results=merger.merge_multiple_files(merger.outputs_primary,merger.outputs_ancillary)
 
    self.assertEqual(results['daynumber'],daynumber)
    for name in ['field_1','field_2','total_uncertainty_field_1','total_uncertainty_field_2','ancillary_1_for_field_1','ancillary_1_for_field_2','ancillary_2_for_field_1','ancillary_2_for_field_2']:
      numpy.testing.assert_array_almost_equal(results[name].data,EXPECTED_MERGED_FIELDS[name],err_msg="Testing correct mergin of "+name+".")
      numpy.testing.assert_array_equal(results[name].mask,EXPECTED_MERGED_FIELDS_MASKS[name],err_msg="Testing correct mask mergin of "+name+".")
    
    # Testing the appending procedure 

    for dictionary in APPENDED_TEST_FIELDS:
      for name in dictionary.keys():
	numpy.testing.assert_array_almost_equal(results[name].data,dictionary[name].data,err_msg="Testing correct appending of "+name+".")
	numpy.testing.assert_array_equal(results[name].mask,dictionary[name].mask,err_msg="Testing correct mask mergin of "+name+".")
	self.assertEqual(results[name].fill_value,dictionary[name].fill_value,"Testing correct fille value")

  def test_merge_outputs_from_multiple_daily_sources(self):
    # Testing correct formats for output files

    basepath = '/gws/nopw/j04/eustace/data/internal/satgrid_lst/test_0.25_0.25/'
    outputs_primary=['tas_ocean_eustace_0_AATSR_20020714.nc','tas_ocean_eustace_0_ATSR2_20020714.nc','tas_ocean_eustace_0_AATSR_20030623.nc']
    outputs_ancillary=['ancillary_ocean_eustace_0_AATSR_20020714.nc','ancillary_ocean_eustace_0_ATSR2_20020714.nc','ancillary_ocean_eustace_0_AATSR_20030623.nc']
    
    suffix='copy.'

    outputs_primary_copy = []
    outputs_ancillary_copy = []
    
    for name_primary, name_ancillary in zip(outputs_primary,outputs_ancillary):
      oldname_primary = basepath+name_primary
      oldname_ancillary = basepath+name_ancillary
      newname_primary = basepath+suffix+name_primary
      newname_ancillary = basepath+suffix+name_ancillary
      copyfile(oldname_primary, newname_primary)
      copyfile(oldname_ancillary, newname_ancillary)
      outputs_primary_copy.append(newname_primary)
      outputs_ancillary_copy.append(newname_ancillary)

    # merging outputs derived from multiple daily sources
    primary_fields=[TAS]
    primary_uncertainties=[TASUNCERTAINTY]
    ancillary_fields= [SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_RANDOM,
		       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED2,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC2,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER0,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER1,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER2,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER3,
                       SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER4]
    # correlation indexes for ancillary_fields along file direction
    list_of_correlation_indexes=[[0.],[1.],[0.],[0.],[0.],[0.],[0.],[0.],[0.],[0.]]
    merger=FilesMerger(outputs_primary_copy, outputs_ancillary_copy, numpy.array([2,1]),primary_fields,primary_uncertainties,ancillary_fields,list_of_correlation_indexes)
    merger.merge_outputs_from_multiple_daily_sources(__name__, 'Met Office',str(1))
    
    source_labels=['AATSR','ATSR2']
    final_primary = ['copy.tas_ocean_eustace_0_20020714.nc','copy.tas_ocean_eustace_0_20030623.nc']
    final_ancillary = ['copy.ancillary_ocean_eustace_0_20020714.nc','copy.ancillary_ocean_eustace_0_20030623.nc']

    # First case, we had merging
    dataset = Dataset(basepath+final_primary[0],'r')

    date = TimeBaseDays(datetime(1850,01,01)).datetime_to_number(datetime(2002,07,14))
    self.assertEqual(dataset.variables['time'][...][0],date)
    self.assertEqual(dataset.institution,'Met Office')
    self.assertEqual(dataset.comment,'input data preprocessed from '+source_labels[0]+', '+source_labels[1]+' retrievals.')
    self.assertEqual(dataset.source,'EUSTACE Catalogue 1')
    for variable in primary_fields+primary_uncertainties:
      self.assertTrue(dataset.variables.has_key(variable.name))
    dataset.close()

    dataset = Dataset(basepath+final_ancillary[0],'r')

    self.assertEqual(dataset.variables['time'][...][0],date)
    self.assertEqual(dataset.institution,'Met Office')
    self.assertEqual(dataset.comment,'input data preprocessed from '+source_labels[0]+', '+source_labels[1]+' retrievals.')
    self.assertEqual(dataset.source,'EUSTACE Catalogue 1')
    for variable in primary_fields+primary_uncertainties:
	for source in source_labels:
	  self.assertTrue(dataset.variables.has_key(variable.name+'_'+source))
    for variable in ancillary_fields:
	self.assertTrue(dataset.variables.has_key(variable.name))
	for source in source_labels:
	  self.assertTrue(dataset.variables.has_key(variable.name+'_'+source))
    dataset.close()

    # Second case, we just had renaming

    dataset = Dataset(basepath+final_primary[1],'r')

    date = TimeBaseDays(datetime(1850,01,01)).datetime_to_number(datetime(2003,06,23))
    self.assertEqual(dataset.variables['time'][...][0],date)
    self.assertEqual(dataset.institution,'Met Office')
    self.assertEqual(dataset.comment,'input data preprocessed from '+source_labels[0]+' retrievals.')
    for variable in primary_fields+primary_uncertainties:
      self.assertTrue(dataset.variables.has_key(variable.name))
    dataset.close()

    dataset = Dataset(basepath+final_ancillary[1],'r')

    self.assertEqual(dataset.variables['time'][...][0],date)
    self.assertEqual(dataset.institution,'Met Office')
    self.assertEqual(dataset.comment,'input data preprocessed from '+source_labels[0]+' retrievals.')
    for variable in ancillary_fields:
	self.assertTrue(dataset.variables.has_key(variable.name))
    dataset.close()

    for name in final_primary+final_ancillary:
      os.remove(basepath+name)

    # Last sanity check: what if we expect an output with a name different from that obtained from the merging?
    final_primary = [basepath+'copy.tas_ocean_eustace_gramnbodk_0_20020714.nc',basepath+'copy.tas_ocean_faafaf_eustace_0_20030623.nc']
    final_ancillary = [basepath+'copy.ancillary_ocean_122eustace_0_20020714.nc',basepath+'copy.ancillary_fiuuu_ocean_eustace_0_20030623.nc']

    self.assertRaises(ValueError,merger.check_produced_output,final_primary,final_ancillary)
