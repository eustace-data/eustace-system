"""Run the given command line when preprocessing a rose script, and return the result as a string."""

__version__ = "$Revision: 787 $"
__author__ = "Joel R. Mitchelson"

import subprocess
import unittest

def runcommand(commandline):
  """Run the specified command line."""

  return subprocess.Popen(commandline.split(), stdout=subprocess.PIPE).communicate()[0]


class TestRunCommand(unittest.TestCase):
  """Test that runcommand works."""

  def test_runcommand(self):    
    self.assertEqual('534\n', runcommand('echo 534'))
