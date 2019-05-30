"""Round floating point to nearest higher integer and return as integer."""

__version__ = "$Revision: 787 $"
__author__ = "Joel R. Mitchelson"

import unittest
import math

def roundup(x):
  """Compute maximum"""
  
  return int(math.ceil(x))


class TestRoundup(unittest.TestCase):
  """Test that roundup works."""
  
  def test_roundup(self):    
    self.assertEqual(10978, roundup(10977.6))
    self.assertEqual(3, roundup(3.0))
    self.assertEqual(3, roundup(2.0001))
    self.assertEqual(8, 8)
  
