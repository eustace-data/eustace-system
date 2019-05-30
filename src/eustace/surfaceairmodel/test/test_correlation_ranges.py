import unittest
from eustace.outputformats.outputvariable import OutputVariable
from ..correlation_ranges import CorrelationRanges

class TestCorrelationRanges(unittest.TestCase):
  
  def test_init(self):

    self.assertRaises(ValueError,CorrelationRanges,3,['c','d'],['d','e'],['f','g'])
    self.assertRaises(ValueError,CorrelationRanges,['a','b'],3,['d','e'],['f','g'])
    self.assertRaises(ValueError,CorrelationRanges,['a','b'],['c','d'],3,['f','g'])
    self.assertRaises(ValueError,CorrelationRanges,['a','b'],['c','d'],['d','e'],3)
    
    A=CorrelationRanges(['a','b'],['c','d'],[1,2],['f','g'])
    self.assertEqual(A.keys,['a','b'])
    self.assertEqual(A.units,['c','d'])
    self.assertEqual(A.length_values,[1,2])
    self.assertEqual(A.time_values,['f','g'])

  def test_length_time_scales_dictionary(self):

    A=CorrelationRanges(['a','b'],[1,2],[4,'e'],['f','g'])
    dictionary=A.length_time_scales_dictionary()
    self.assertEqual(A.keys,dictionary.keys())
    
    for i,j,z in zip(A.keys,A.length_values,A.time_values):
      sub_dictionary=dictionary[i]
      self.assertEqual(j,sub_dictionary['length scale'])
      self.assertEqual(z,sub_dictionary['time scale'])
      self.assertEqual(A.units[0],sub_dictionary['length scale units'])
      self.assertEqual(A.units[1],sub_dictionary['time scale units'])

  def test_update_correlated_uncertainty_ranges(self):
    
    A=CorrelationRanges(['test','b'],['g','h'],[1.2,3.67],[3,4.4])
    TEST_VARIABLE=TIME = OutputVariable(name='test',dtype=int, fill_value=0)

    A.update_correlated_uncertainty_ranges(TEST_VARIABLE)
    self.assertEqual(TEST_VARIABLE.length_scale,1.2)
    self.assertEqual(TEST_VARIABLE.time_scale,3)
    self.assertEqual(TEST_VARIABLE.length_scale_units,'g')
    self.assertEqual(TEST_VARIABLE.time_scale_units,'h')
