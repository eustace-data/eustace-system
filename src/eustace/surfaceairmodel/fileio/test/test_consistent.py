import unittest
import numpy
import tempfile
import netCDF4
from ..consistent import ConsistentModelOutputNetCDF, DatasetAttributesConsistentModelOutput, ObservationSourceSatstace
from eustace.analysis.observationsource import ObservationSource

class TestConsistentModelOutputNetCDF(unittest.TestCase):

    def test_write_consistent(self):

        # Make data
        testfile = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmode.fileio.test.test_consistent.', suffix='.nc')
        testfield_data = numpy.ones((720,1440)) * 7.12
        testfield_mask = numpy.zeros((720,1440), dtype=numpy.bool)
        testfield_mask[1,3] = True
        testfield_data[7,8] = 4.2
        testfield = numpy.ma.masked_array(data=testfield_data, mask=testfield_mask)
        testdata = { 'daynumber': 576, 'tasuncertainty': testfield }
        attributes = DatasetAttributesConsistentModelOutput('My Module', 'This Place', '5672')
        writer = ConsistentModelOutputNetCDF(attributes, [ ])
        writer.write_primary(testfile.name, testdata)
        
        # Check it
        netcdf = netCDF4.Dataset(testfile.name, 'r')
        self.assertEqual(netcdf.institution, 'This Place')
        numpy.testing.assert_equal(testfield_mask, netcdf.variables['tasuncertainty'][0,:,:].mask)
        result_field = netcdf.variables['tasuncertainty'][0,:,:]
        self.assertAlmostEqual(7.12, result_field.data[   0,   0])
        self.assertAlmostEqual(7.12, result_field.data[   0,   1])
        self.assertAlmostEqual( 4.2, result_field.data[   7,   8])
        self.assertAlmostEqual(7.12, result_field.data[   0,1439])
        self.assertAlmostEqual(7.12, result_field.data[   0,1439])
        self.assertAlmostEqual(7.12, result_field.data[ 719,1438])
        self.assertAlmostEqual(7.12, result_field.data[ 719,1439])
        netcdf.close()
        
    def test_post_process(self):
        # Check correct error raise
        attributes = DatasetAttributesConsistentModelOutput('My Module', 'This Place', '5672')
        writer = ConsistentModelOutputNetCDF(attributes, [ ])
	input_dict = {'results':{}, 'surface':'Dada'}
	self.assertRaises(ValueError, writer.post_process, **input_dict)

class TestObservationSourceSatstace(unittest.TestCase):

    def test_read(self):

        # Make data using writer
        testfile = tempfile.NamedTemporaryFile(prefix='eustace.surfaceairmode.fileio.test.test_consistent.', suffix='.nc')

        testtas_data = numpy.ones((720,1440)) * 277.12
        testtas_mask = numpy.zeros((720,1440), dtype=numpy.bool)
        testtas_mask[1,3] = True
        testtas_data[7,8] = 274.2
        testtas = numpy.ma.masked_array(data=testtas_data, mask=testtas_mask)

        testuncertainty_data = numpy.ones((720,1440)) * 7.12
        testuncertainty_mask = numpy.zeros((720,1440), dtype=numpy.bool)
        testuncertainty_mask[1,3] = True
        testuncertainty_data[7,8] = 4.2
        testuncertainty = numpy.ma.masked_array(data=testuncertainty_data, mask=testuncertainty_mask)

        testdata = { 'daynumber': 576, 'tas': testtas, 'tasuncertainty': testuncertainty }
        attributes = DatasetAttributesConsistentModelOutput('My Module', 'This Place', '5672')
        writer = ConsistentModelOutputNetCDF(attributes, [ ])
        writer.write_primary(testfile.name, testdata)

        # Read back in
        source = ObservationSourceSatstace([ ObservationSource.TMEAN ], testfile.name)
        obs = source.observations(ObservationSource.TMEAN)

        # Check same as starting point
        numpy.testing.assert_equal(testtas_mask.ravel(), obs.mask)
        numpy.testing.assert_equal(testtas_data[~testtas_mask].ravel(), obs.measurement[~obs.mask])
        numpy.testing.assert_equal(testuncertainty_mask.ravel(), obs.mask)
        numpy.testing.assert_equal(testuncertainty_data[~testuncertainty_mask].ravel(), obs.uncorrelatederror[~obs.mask])
