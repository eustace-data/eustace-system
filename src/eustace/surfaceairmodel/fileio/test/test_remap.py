"""Test NetCDF remapping."""

import unittest
import tempfile
import netCDF4
import numpy
from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eustace.outputformats.outputvariable import OutputVariable
from eustace.outputformats import definitions
from ..remap import RemapNetCDFSpecCopy
from ..remap import RemapNetCDFSpecDayNumber
from ..remap import RemapNetCDF

class TestRemapNetCDFSpecification(unittest.TestCase):

    def test_init(self):

        a = RemapNetCDFSpecCopy('thatone', 'thisone', 27.9)
        self.assertEqual('thisone', a.inputname)
        self.assertEqual('thatone', a.outputname)
        self.assertAlmostEqual(27.9, a.offset)

        b = RemapNetCDFSpecCopy('good', 'bad')
        self.assertEqual(b.inputname, 'bad')
        self.assertEqual(b.outputname, 'good')
        self.assertIsNone(b.offset)

class TestRemapNetCDF(unittest.TestCase):

    def test_read_fields(self):

        # Make test data
        testfile = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmodel.fileio.test.test_remap.', suffix='.nc')
        testdata = FileBuilderGlobalField(testfile.name, 10000, 'testdata', 'NOVERSION', 'jam', 'someplace', '', '', '')
        testdata.add_global_field(OutputVariable('jam', numpy.float32, -2000), 182.2 * numpy.ones(definitions.GLOBAL_FIELD_SHAPE))
        testdata.add_global_field(OutputVariable('apples', numpy.float32, -2000), 37.8 * numpy.ones(definitions.GLOBAL_FIELD_SHAPE))

        # Remapper
        remapper = RemapNetCDF([ RemapNetCDFSpecDayNumber(), RemapNetCDFSpecCopy('cream', 'jam', -29.1), RemapNetCDFSpecCopy('oranges', 'apples') ] )

        # Get results
        results = remapper.read_fields(testfile.name)

        # Should have 3 fields (oranges, cream, and the daynumber)
        self.assertEqual(3, len(results))

        # Check daynumber
        self.assertEqual(10000, results['daynumber'])

        # Check global fields are as expected
        numpy.testing.assert_almost_equal(results['cream'], 153.1 * numpy.ones(definitions.GLOBAL_FIELD_SHAPE[1:]), decimal=4)
        numpy.testing.assert_almost_equal(results['oranges'], 37.8 * numpy.ones(definitions.GLOBAL_FIELD_SHAPE[1:]), decimal=4)
        
