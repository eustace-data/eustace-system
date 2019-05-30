"""Tests for covariate element."""

import json
import numpy
import os
import scipy.sparse
import tempfile
import unittest

from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariatePrior
from eustace.analysis.advanced_standard.elements.covariate import LoadCovariateElement
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT
from eustace.analysis.advanced_standard.elements.geography_based import GeographyBasedElement
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement

class TestCovariateHyperparameters(unittest.TestCase):

    def test_init(self):
        
        h = CovariateHyperparameters(-386.26)
        self.assertEqual(-386.26, h.value)

    def test_get_array(self):

        numpy.testing.assert_equal(CovariateHyperparameters(222.2).get_array(), [ 222.2 ])

    def test_set_array(self):

        h = CovariateHyperparameters(9.0)
        h.set_array(numpy.array([ 999.0 ]))
        self.assertEqual(999.0, h.value)
        numpy.testing.assert_equal(h.get_array(), [ 999.0 ])

class TestCovariatePrior(unittest.TestCase):
                
    def test_init(self):

        p = CovariatePrior(CovariateHyperparameters(8.8), 7)
        self.assertEqual(8.8, p.hyperparameters.value)
        self.assertEqual(7, p.number_of_state_parameters)

    def test_prior_number_of_state_parameters(self):

        self.assertEqual(7, CovariatePrior(CovariateHyperparameters(0.1), 7).prior_number_of_state_parameters())

    def test_prior_precision(self):

        Q = CovariatePrior(CovariateHyperparameters(-0.5 * numpy.log(1.1)), 7).prior_precision()
        self.assertEqual(SPARSEFORMAT, Q.getformat())
        self.assertEqual((7,7), Q.shape)
        self.assertEqual(7, Q.nnz)
        numpy.testing.assert_almost_equal(Q.todense(), numpy.diag([ 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1 ]))

    def test_prior_precision_derivative(self):

        dQ = CovariatePrior(CovariateHyperparameters(-0.5 * numpy.log(1.1)), 7).prior_precision_derivative(0)
        self.assertEqual(SPARSEFORMAT, dQ.getformat())
        self.assertEqual((7, 7), dQ.shape)
        self.assertEqual(7, dQ.nnz)
        self.assertAlmostEqual(-2.2, dQ[1,1])
        with self.assertRaises(ValueError):
            CovariatePrior(CovariateHyperparameters(-0.5 * numpy.log(1.1)), 7).prior_precision_derivative(1)

class TestLoadCovariateElement(unittest.TestCase):
      def setUp(self):
	  self.json_descriptor_wrong = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.covariate.wrong_json', suffix='.json',delete=False)
	  self.json_descriptor_good = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.covariate.good_json', suffix='.json',delete=False)
	  
	  self.val = 1.15
	  self.wrong_dictionary = {"sbrambod":0}
	  self.good_dictionary = {"latitude_harmonics":{"element":{"python_class":"eustace.analysis.advanced_standard.elements.latitudeharmonics.LatitudeHarmonicsElement"},
			                                "hyperparameters":{"hyperparameter_1":self.val,
                                                                           "hyperparameter_2":self.val,
                                                                           "hyperparameter_3":self.val,
                                                                           "hyperparameter_4":self.val}},
                                  "altitude":{"element":{"python_class":"eustace.analysis.advanced_standard.elements.geography_based.GeographyBasedElement",
                                              "parameters":{"filename":"/gws/nopw/j04/eustace/data/internal/climatology_covariates/DEM_global_0.25_0.25.nc", 
                                                            "latitude_label":"lat", 
                                                            "longitude_label":"lon", 
                                                            "covariate_label":"dem",
                                                            "rescale_factor": 1.0}},
                                              "hyperparameters":{"value":self.val}}}

	  json.dump(self.wrong_dictionary, self.json_descriptor_wrong)
 	  json.dump(self.good_dictionary, self.json_descriptor_good)
	  self.json_descriptor_wrong.close()
	  self.json_descriptor_good.close()

      def tearDown(self):
          self.wrong_dictionary = None
          self.good_dictionary = None
	  os.remove(self.json_descriptor_wrong.name)
	  os.remove(self.json_descriptor_good.name)


      def test_init(self):
	    A = LoadCovariateElement(self.json_descriptor_wrong.name)
	    self.assertDictEqual(self.wrong_dictionary, A.data, 'Init testing: loading of wrong json file')

	    B = LoadCovariateElement(self.json_descriptor_good.name)
	    self.assertDictEqual(self.good_dictionary, B.data, 'Init testing: loading of good json file')

      def test_check_keys(self):
	    A = LoadCovariateElement(self.json_descriptor_wrong.name)
	    self.assertRaises(A.check_keys,msg='Testing the \"check_keys\" method')

      def test_load_covariate(self):
	    B = LoadCovariateElement(self.json_descriptor_good.name)
	    latitude_harmonics_element = B.load_covariate('latitude_harmonics')
	    self.assertTrue(isinstance(latitude_harmonics_element, LatitudeHarmonicsElement ), msg='Testing correct instantiation of LatitudeHarmonicsElement class')
	    self.assertEqual(4, len(latitude_harmonics_element.HARMONICS))

	    altitude_element = B.load_covariate('altitude')
	    self.assertTrue(isinstance(altitude_element, GeographyBasedElement ), msg='Testing correct instantiation of GeographyBasedElement class')
	    self.assertEqual(altitude_element.filename, self.good_dictionary['altitude']['element']['parameters']['filename'], msg='Testing correct filename attribute of GeographyBasedElement object')
	    self.assertEqual(altitude_element.latitude_label, self.good_dictionary['altitude']['element']['parameters']['latitude_label'], msg='Testing correct latitude_label attribute of GeographyBasedElement object')
	    self.assertEqual(altitude_element.longitude_label, self.good_dictionary['altitude']['element']['parameters']['longitude_label'], msg='Testing correct longitude_label attribute of GeographyBasedElement object')
	    self.assertEqual(altitude_element.covariate_label, self.good_dictionary['altitude']['element']['parameters']['covariate_label'], msg='Testing correct covariate_label attribute of GeographyBasedElement object')
	    self.assertEqual(altitude_element.NUMBER_OF_STATE_PARAMETERS, 1)
	    self.assertTrue(hasattr(altitude_element, 'latitude'), msg='Testing correct latitude loading from physical covariate file')
	    self.assertTrue(hasattr(altitude_element, 'longitude'), msg='Testing correct longitude loading from physical covariate file')
	    self.assertTrue(hasattr(altitude_element, 'covariate'), msg='Testing correct covariate loading from physical covariate file')

      def test_load_hyperparameters(self):
  	    B = LoadCovariateElement(self.json_descriptor_good.name)
	    latitude_harmonics_hyperparameters = B.load_hyperparameters('latitude_harmonics')

	    self.assertTrue(isinstance(latitude_harmonics_hyperparameters, CombinationHyperparameters), msg='Testing correct instantiation of CombinationHyperparameters class')
	    self.assertEqual(4, len(latitude_harmonics_hyperparameters.elementparameters), msg='Testing correct number of hyperparameter objects')

	    for index,hyperparameter_object in enumerate(latitude_harmonics_hyperparameters.elementparameters):
		self.assertTrue(isinstance(hyperparameter_object, CovariateHyperparameters), msg='Testing correct instantiation of CovariateHyperparameters class, object '+str(index))
	    
	    numpy.testing.assert_array_equal(numpy.log(self.val), latitude_harmonics_hyperparameters.get_array(), err_msg='Testing latitude harmonics hyperparameters values')

	    altitude_hyperparameters = B.load_hyperparameters('altitude')

	    self.assertTrue(isinstance(altitude_hyperparameters, CovariateHyperparameters ), msg='Testing correct instantiation of CovariateHyperparameters class')	    
	    numpy.testing.assert_array_equal(numpy.log(self.val), altitude_hyperparameters.get_array(), err_msg='Testing altitude hyperparameter value')

      def test_load_covariates_and_hyperparameters(self):
	    B = LoadCovariateElement(self.json_descriptor_good.name)
	    list_of_covariate_elements, list_of_covariate_hyperparameters = B.load_covariates_and_hyperparameters()

	    self.assertEqual(len(list_of_covariate_elements), 2, msg='Testing number of loaded covariates')
	    self.assertEqual(len(list_of_covariate_hyperparameters), 2, msg='Testing number of loaded hyperparameters')