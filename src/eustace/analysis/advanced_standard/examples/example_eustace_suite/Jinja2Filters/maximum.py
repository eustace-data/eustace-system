"""Maximum of a list."""

__version__ = "$Revision: 787 $"
__author__ = "Joel R. Mitchelson"

import unittest

def maximum(x):
  """Compute maximum"""
  
  return max(x)


class TestMaximum(unittest.TestCase):
  """Test that maximum works."""
  
  def test_maximum(self):    
    self.assertEqual(7, maximum([ 2, 7, 1, 6.0 ]))
    self.assertEqual(3, maximum([ 3 ]))
  
