import unittest
import numpy
from eustace.preprocess.fileio.height_adjustment import AltitudeAdjustment, LAPSE_RATE, LAPSE_RATE_UNCERTAINTY

class TestObservationSourceInsituLand(unittest.TestCase):
    
    def setUp(self):
    
        # setup using the analysis' DEM covariate file
        filename = "/gws/nopw/j04/eustace/data/internal/climatology_covariates/DEM_global_0.25_0.25.nc"
        latitude_label = "lat"
        longitude_label = "lon" 
        covariate_label = "dem"
        rescale_factor = 0.001 # rescale altitude map to to K km^-1

        # initialise the model
        self.adjustment_model = AltitudeAdjustment(filename, latitude_label, longitude_label, covariate_label, rescale_factor)
        self.adjustment_model.load_reference_map()
    
    def test_zero_adjustment(self):
        
        # arbitrary test locations
        latitudes  = numpy.array([52.0, 52.0, 52.0, 52.0])
        longitudes = numpy.array([0.0, 1.0, 2.0, 3.0])
        locations = numpy.array(zip( latitudes, longitudes ))

        # set altitudes to equal that returned from the reference map
        altitudes = self.adjustment_model.reference_altitudes(locations)

        # check that the output adjustment equals zero
        numpy.testing.assert_almost_equal(numpy.zeros(4), self.adjustment_model.adjustment(altitudes, locations))
        
        # check that the adjustment uncertainty equals zero
        numpy.testing.assert_almost_equal(numpy.zeros(4), self.adjustment_model.adjustment_uncertainty(altitudes, locations))
        
        
    def test_adjustment_from_sealevel(self):

        # arbitrary test locations
        latitudes  = numpy.array([52.0, 52.0, 52.0, 52.0])
        longitudes = numpy.array([0.0, 1.0, 2.0, 3.0])
        locations = numpy.array(zip( latitudes, longitudes ))

        # set altitudes to equal that returned from the reference map
        altitudes = numpy.zeros((4))
        reference_altitudes = self.adjustment_model.reference_altitudes(locations)

        # check that the output adjustment is the reference surface altitude times the lapse rate
        desired_adjustment = -LAPSE_RATE * (altitudes - reference_altitudes)
        numpy.testing.assert_almost_equal( desired_adjustment, self.adjustment_model.adjustment(altitudes, locations))
        
        # check that the adjustment uncertainty equals the absolute value of the reference surface altitude times the lapse rate uncertainty
        desired_uncertainty = LAPSE_RATE_UNCERTAINTY * numpy.abs(altitudes - reference_altitudes)
        numpy.testing.assert_almost_equal(desired_uncertainty, self.adjustment_model.adjustment_uncertainty(altitudes, locations))

    def test_unit_adjustment(self):
        
        # arbitrary test locations
        latitudes  = numpy.array([52.0, 52.0, 52.0, 52.0])
        longitudes = numpy.array([0.0, 1.0, 2.0, 3.0])
        locations = numpy.array(zip( latitudes, longitudes ))

        # set altitudes to 1km lower than the reference altitudes returned from the reference map
        altitudes = self.adjustment_model.reference_altitudes(locations) - 1.0

        # check that the output is adjustmented upwards by 1km times the lapse rate
        desired_adjustment = LAPSE_RATE * numpy.ones(4)
        numpy.testing.assert_almost_equal(desired_adjustment, self.adjustment_model.adjustment(altitudes, locations))
        
        # check that the adjustment uncertainty equals the lapse rate uncertainty
        desired_uncertainty = LAPSE_RATE_UNCERTAINTY * numpy.ones(4)
        numpy.testing.assert_almost_equal(desired_uncertainty, self.adjustment_model.adjustment_uncertainty(altitudes, locations))
        
        
        # set altitudes to 1km higher than the reference altitudes returned from the reference map
        altitudes = self.adjustment_model.reference_altitudes(locations) + 1.0

        # check that the output is adjustment downwards by one km time the lapse rate
        desired_adjustment = LAPSE_RATE * -1.0 * numpy.ones(4)
        numpy.testing.assert_almost_equal(desired_adjustment, self.adjustment_model.adjustment(altitudes, locations))
        
        # check that the adjustment uncertainty equals the lapse rate uncertainty
        desired_uncertainty = LAPSE_RATE_UNCERTAINTY * numpy.ones(4)
        numpy.testing.assert_almost_equal(desired_uncertainty, self.adjustment_model.adjustment_uncertainty(altitudes, locations))
