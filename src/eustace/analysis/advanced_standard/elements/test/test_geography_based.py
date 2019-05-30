"""Tests for gegraphy-based element."""

import unittest
import numpy
import os
import scipy.sparse
import tempfile

from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.advanced_standard.elements.element import SPARSEFORMAT
from eustace.analysis.advanced_standard.elements.geography_based import GeographyBasedCovariateFunction
from eustace.analysis.advanced_standard.elements.geography_based import GeographyBasedElement
from eustace.analysis.advanced_standard.elements.geography_based import GeographyBasedElementDesign
from eustace.analysis.advanced_standard.elements.geography_based import GeographyBasedPrior
from eustaceconfig import WORKSPACE_PATH
from netCDF4 import Dataset

ALTITUDE_RELATIVE_PATH = 'data/internal/climatology_covariates/DEM_global_0.25_0.25.nc'
COASTAL_INFLUENCE_RELATIVE_PATH = 'data/internal/climatology_covariates/coastal_influence.test.0.25_0.25.nc'
class SimulatedObservationStructure(ObservationStructure):
	    
    def location_polar_coordinates(self):

	return numpy.array([[0.5, 1.],[0.2, 1.1],[0.25, 0.5],[0.7, -0.2]])


class TestGeographyBasedCovariateFunction(unittest.TestCase):

    def setUp(self):

        nx, ny = (5, 3)
        x = numpy.linspace(0, 1, nx) 
        y = numpy.linspace(0, 0.5, ny)
        self.longitude, self.latitude = numpy.meshgrid(x, y)
        self.covariate=numpy.array([[1009., 2091, 300., 10. , 24.],
                                    [1415., 121., 2231., 342., 5550.],
                                    [223., 4317., 1819., 949., 3426.]])
                                    
        self.rescale_factor = 1.0
        
    def tearDown(self):
        self.latitude = None
        self.longitude = None
        self.covariate = None

    def test_init(self):
        """Testing class initialization"""

        a = GeographyBasedCovariateFunction('a','b','c','d')
        self.assertEqual(a.latitude,'a')
        self.assertEqual(a.longitude,'b')
        self.assertEqual(a.covariate,'c')
        self.assertEqual(a.rescale_factor, 'd')
        
    def test_bilinear_interpolation(self):
        """Testing bilinear interpolation"""

        nx, ny = (5, 3)
        x = numpy.linspace(0, 1, nx) 
        y = numpy.linspace(0, 1, ny)
        longitude, latitude = numpy.meshgrid(x, y)
        covariate, locations =None, None
  
        # Testing correct exception raise if nodes spacing not uniform
        self.assertRaises(GeographyBasedCovariateFunction.bilinear_interpolation,covariate, latitude, longitude, locations)

        # 2x2 grid => 4x4 extendend: 
        #
        # extended covariate array is          extended grid array is
            #
        #       m | m  m | m                          (-1,-1)| (-1,0) (-1,1)| (-1,2)
        #      --------------                         -----------------------------
        #       b | a  b | a                          (0,-1) | (0,0)  (0,1) | (0,2)
        #       d | c  d | c                          (1,-1) | (1,0)  (1,1) | (1,2)
        #      --------------                         -----------------------------
        #       M | M  M | M                          (2,-1) | (2,0)  (2,1) | (2,2)
        #
        #locations points at:
        # 
        # [0., 0. ]  :  node
        # [0., 1.]   :  node
        # [1., 0.]   :  node
        # [1., 1.]   :  node
        # [0.5, 0.5] :  not node, contributing nodes are [0., 0.], [1., 0.], [0., 1.], [1., 1.], each contributing with 0.25
        # [0., 0.5]  :  not node, contributing nodes are [0., 0.], [0., 1.], with weights 0.5, 0.5
        # [0.5, 0.]  :  not node, contributing nodes are [0., 0.], [1., 0.], with weights 0.5, 0.5
        # [0.3, 0.7] :  not node, contributing nodes are [0., 0.], [0., 1.], [1., 0.], [1., 1.], with weights 0.3*0.7, 0.7*0.7, 0.3*0.3, 0.3*0.7
        # [-0.4, 0.6]:  not node, contributing nodes are [0., 0.], [-1., 0.], [-1., 1.], [0., 1.], with weights 0.6*0.4, 0.4*0.4, 0.6*0.4, 0.6*0.6
        # [1.8, -0.2]:  not node, contributing nodes are [1., 0.], [2., 0.], [1., -1.], [2., -1.], with weights 0.8*0.2, 0.8*0.8, 0.2*0.2, 0.8*0.2
        # [0.9, 1.9] :  not node, contributing nodes are [1., 1.], [1., 2.], [0., 1.], [0., 2.], with weights 0.1*0.9, 0.9*0.9, 0.1*0.1, 0.9*0.1

        nx, ny = (2, 2)
        x = numpy.linspace(0, 1, nx) 
        y = numpy.linspace(0, 1, ny)
        longitude, latitude = numpy.meshgrid(x, y)
        locations = numpy.array([[0., 0.],[0., 1.],[1., 0.],[1., 1.],[0.5, 0.5],[0.5, 0.],[0., 0.5],[0.3, 0.7],[-0.4, 0.6],[1.8, -.2],[0.9, 1.9]])

        # 1st test case: uniform covariate
        # extended covariate array is          
        #
        #       1 | 1  1 | 1  
        #      -------------- 
        #       1 | 1  1 | 1  
        #       1 | 1  1 | 1  
        #      -------------- 
        #       1 | 1  1 | 1  

        covariate = numpy.zeros([2,2])+1
        result_vector = GeographyBasedCovariateFunction.bilinear_interpolation(covariate, latitude, longitude, locations)
        numpy.testing.assert_array_almost_equal( numpy.array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]), result_vector, err_msg  ="Testing covariate design vector, trivial case")

        # 2nd test case: non uniform covariate
        # extended covariate array is          
            #
        #       1.5 | 1.5  1.5 | 1.5  
        #      --------------------- 
        #        2  |  1    2  | 1  
        #        4  |  3    4  | 3  
        #      --------------------- 
        #      3.5  | 3.5  3.5 | 3.5

        covariate = numpy.array([[1., 2.],[ 3., 4.]])
        result_vector = GeographyBasedCovariateFunction.bilinear_interpolation(covariate, latitude, longitude, locations)
        numpy.testing.assert_array_almost_equal( numpy.array([1., 2., 3., 4., 2.5, 2.0, 1.5, 2.3, 1.56, 3.44, 2.9]), result_vector, err_msg  ="Testing covariate design vector, non trivial first case")

    def test_compute(self):
        """Testing the compute method on a fairly more real dataset"""

        # extended grid array is
            #
        #   (-.25,-.25)|   (-.25,0)   (-.25,.25)  (-.25,0.5)  (-.25,0.75)  (-.25,1.) |   (-.25,1.25)
        #   --------------------------------------------------------------------------------------
        #     (0,-.25) |     (0,0)      (0,.25)     (0,0.5)     (0,0.75)     (0,1.) |     (0,1.25)
        #   (.25,-.25) |   (.25,0)    (.25,.25)   (.25,0.5)   (.25,0.75)   (.25,1.) |   (.25,1.25)
            #   (0.5,-.25) |   (0.5,0)    (0.5,.25)   (0.5,0.5)   (0.5,0.75)   (0.5,1.) |   (0.5,1.25)
        #   --------------------------------------------------------------------------------------
        #  (0.75,-.25) |  (0.75,0)   (0.75,.25)  (0.75,0.5)  (0.75,0.75)  (0.75,1.) |  (0.75,1.25)
        
        # extended covariate array is          
        #
        #                     638.1428833  .... 
        #                     -------------------------------------------------
        #                       24. | 1009.   2091.   300.   10.    24. | 1009.
        #                     5550. | 1415.    121.  2231.  342.  5550. | 1415.
        #                     3426. |  223.   4317.  1819.  949.  3426. |  223.
        #                     -------------------------------------------------
        #                     2054.71435547  ....

        a = GeographyBasedCovariateFunction(self.latitude.astype(numpy.float32), self.longitude.astype(numpy.float32), self.covariate.astype(numpy.float32), self.rescale_factor)
        locations = numpy.array([[0., 0.],[.5, 0.],[0.5, 1.],[0.2, 1.1],[0.25, 0.5],[0.7, -0.2],[-0.16, 0.2],[.75, .75],[-0.25, 0.],[0.5, 1.25],[-0.19, .1]])
        result_vector = a.compute(locations)
        numpy.testing.assert_array_almost_equal( numpy.array([1009., 223., 3426., 3200.4, 2231, 2200.851484, 1083.267445, 2054.714355, 638.142883, 223., 831.020591]).reshape(-1,1), result_vector,  decimal = 6, err_msg  ="Testing covariate design matrix, general case")

class TestGeographyBasedElement(unittest.TestCase):
    def setUp(self):
        self.covariate_file = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.elements.test_geography_based', suffix='.nc',delete=True)

        nx, ny = (5, 3)
        x = numpy.linspace(0, 1, nx) 
        y = numpy.linspace(0, 0.5, ny)
        self.longitude, self.latitude = numpy.meshgrid(x, y)
        self.covariate=numpy.array([[1009., 2091, 300., 10. , 24.],
                                    [1415., 121., 2231., 342., 5550.],
                                    [223., 4317., 1819., 949., 3426.]])
        
        dataset=Dataset(self.covariate_file.name,'w','NETCDF4')
        dataset.createDimension('x',3)
        dataset.createDimension('y',5)
        field_latitude=dataset.createVariable('lat',numpy.float32,dimensions=('x','y'))
        field_longitude=dataset.createVariable('lon',numpy.float32,dimensions=('x','y'))
        field_covariate=dataset.createVariable('covariate',numpy.float32,dimensions=('x','y'))
        field_latitude[:]=self.latitude
        field_longitude[:]=self.longitude
        field_covariate[:]=self.covariate
        dataset.close()

    def tearDown(self):
        self.latitude = None
        self.longitude = None
        self.covariate = None
    
    def test_init(self):
        
        a = GeographyBasedElement('a', 'b', 'c', 'd', 'e')
        self.assertFalse(a.isnonlinear())

    def test_load(self):

        a = GeographyBasedElement(self.covariate_file.name,'lat','lon','covariate',1.0)
        a.load()
       
        for original_array, file_array, name in zip([self.latitude, self.longitude, self.covariate],[a.latitude, a.longitude, a.covariate],['latitude', 'longitude', 'covariate']):
          numpy.testing.assert_array_equal(original_array, file_array, err_msg="Testing vaues of "+name+" array")

    def test_element_design(self):

        geography = GeographyBasedElement(self.covariate_file.name,'lat','lon','covariate',1.0)
        geography.load()

        design = geography.element_design(SimulatedObservationStructure())
        self.assertTrue(isinstance(design, GeographyBasedElementDesign))
        self.assertTrue(isinstance(design.function_class, GeographyBasedCovariateFunction))
        for original_array, file_array, name in zip([self.latitude, self.longitude, self.covariate],[design.function_class.latitude, design.function_class.longitude, design.function_class.covariate],['latitude', 'longitude', 'covariate']):
            numpy.testing.assert_array_equal(original_array, file_array, err_msg="Testing vaues of "+name+" array")

        self.assertEqual(design.design_number_of_state_parameters(),GeographyBasedElement.NUMBER_OF_STATE_PARAMETERS)

    def test_element_prior(self):

        geography = GeographyBasedElement(self.covariate_file.name,'lat','lon','covariate',1.0)
        value = 3.333
        prior = geography.element_prior(CovariateHyperparameters(value))
        self.assertIsInstance(prior, GeographyBasedPrior)
        
        # As an example check precision - should be diagonals with exp(-2 x hyperparameter)
        precision = prior.prior_precision()
        self.assertEqual(SPARSEFORMAT, precision.getformat())
        self.assertEqual(1, precision.nnz)
        self.assertAlmostEqual(numpy.exp(-2.0*value), precision[0,0])

        precision_derivative = prior.prior_precision_derivative(0)
        self.assertEqual(SPARSEFORMAT, precision_derivative.getformat())
        self.assertEqual(1, precision_derivative.nnz)
        self.assertAlmostEqual(-2.0 * precision[0,0], precision_derivative[0,0])

class TestGeographyBasedElementDesign(unittest.TestCase):

    def setUp(self):
        self.covariate_file = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.elements.test_geography_based', suffix='.nc',delete=True)

        nx, ny = (5, 3)
        x = numpy.linspace(0, 1, nx) 
        y = numpy.linspace(0, 0.5, ny)
        self.longitude, self.latitude = numpy.meshgrid(x, y)
        self.covariate=numpy.array([[1009., 2091, 300., 10. , 24.],
                                    [1415., 121., 2231., 342., 5550.],
                                    [223., 4317., 1819., 949., 3426.]])

        self.rescale_factor = 1.0

        self.function_class =  GeographyBasedCovariateFunction(self.latitude.astype(numpy.float32), self.longitude.astype(numpy.float32), self.covariate.astype(numpy.float32), self.rescale_factor)

    def tearDown(self):
        self.latitude = None
        self.longitude = None
        self.covariate = None
        self.function_class = None

    def test_design_number_of_state_parameters(self):

        self.assertEqual(GeographyBasedElement.NUMBER_OF_STATE_PARAMETERS, GeographyBasedElementDesign(self.function_class, SimulatedObservationStructure()).design_number_of_state_parameters())

    def test_isnonlinear(self):
        
        self.assertFalse(GeographyBasedElementDesign(self.function_class, SimulatedObservationStructure()).isnonlinear())

    def test_design_matrix(self):

        design = GeographyBasedElementDesign(self.function_class, SimulatedObservationStructure())
        m = design.design_matrix()
        self.assertEqual(SPARSEFORMAT, m.getformat())
        numpy.testing.assert_almost_equal(numpy.array([3426., 3200.4, 2231, 2200.8514844]).reshape(-1,1), m.todense())

class TestRealCovariate(unittest.TestCase):
    """Testing the Geography-based module on a real case"""
    
    class SimulatedObservationStructureAltitude(ObservationStructure):
            
        def location_polar_coordinates(self):

            return numpy.array([[-89.875, -179.875], [-89.625, -179.625],[-89.75, -179.750], [-89.627, -179.75]])

    class SimulatedObservationStructureCoastalInfluence(ObservationStructure):
            
        def location_polar_coordinates(self):

            return numpy.array([[-85.125, -159.375], [-85.125, -159.125],[-84.875, -159.125], [-85., -159.2]])


    def setUp(self):
        
        self.altitude_datafile = os.path.join(WORKSPACE_PATH, ALTITUDE_RELATIVE_PATH)
        self.coastal_influence_datafile = os.path.join(WORKSPACE_PATH, COASTAL_INFLUENCE_RELATIVE_PATH)
    
    def test_Altitude_Covariate(self):

        altitude_element = GeographyBasedElement(self.altitude_datafile, 'lat', 'lon', 'dem', 1.0)
        altitude_element.load()

        numpy.testing.assert_array_equal(-89.875, altitude_element.latitude[0,0])
        numpy.testing.assert_array_equal(-84.875, altitude_element.latitude[20,0])
        numpy.testing.assert_array_equal(-89.875, altitude_element.latitude[0,10])
        numpy.testing.assert_array_equal(-84.875, altitude_element.latitude[20,10])

        numpy.testing.assert_array_equal(-179.875, altitude_element.longitude[0,0])
        numpy.testing.assert_array_equal(-179.875, altitude_element.longitude[20,0])
        numpy.testing.assert_array_equal(-177.375, altitude_element.longitude[0,10])
        numpy.testing.assert_array_equal(-177.375, altitude_element.longitude[20,10])

        numpy.testing.assert_array_equal(2779, altitude_element.covariate[0,0])
        numpy.testing.assert_array_equal(2861, altitude_element.covariate[1,1])

        design = altitude_element.element_design(TestRealCovariate.SimulatedObservationStructureAltitude())
        altitude_design = design.design_matrix()
        numpy.testing.assert_almost_equal(numpy.array([2779, 2861, 2820, 2860.344]).reshape(-1,1), altitude_design.todense())
    
    def test_CoastalInfluence_Covariate(self):

        coastalinfluence_element = GeographyBasedElement(self.coastal_influence_datafile, 'lat', 'lon', 'coastal_influence', 1.0)
        coastalinfluence_element.load()

        numpy.testing.assert_array_equal(-89.875, coastalinfluence_element.latitude[0,0])
        numpy.testing.assert_array_equal(-84.875, coastalinfluence_element.latitude[20,0])
        numpy.testing.assert_array_equal(-89.875, coastalinfluence_element.latitude[0,10])
        numpy.testing.assert_array_equal(-84.875, coastalinfluence_element.latitude[20,10])

        numpy.testing.assert_array_equal(-179.875, coastalinfluence_element.longitude[0,0])
        numpy.testing.assert_array_equal(-179.875, coastalinfluence_element.longitude[20,0])
        numpy.testing.assert_array_equal(-177.375, coastalinfluence_element.longitude[0,10])
        numpy.testing.assert_array_equal(-177.375, coastalinfluence_element.longitude[20,10])

        numpy.testing.assert_array_almost_equal(2.6666667, coastalinfluence_element.covariate[19,82])
        numpy.testing.assert_array_almost_equal(100.0, coastalinfluence_element.covariate[20,83])

        design = coastalinfluence_element.element_design(TestRealCovariate.SimulatedObservationStructureCoastalInfluence())
        coastalinfluence_design = design.design_matrix()
        numpy.testing.assert_almost_equal(numpy.array([2.6666667, 14.9382706, 100.0, 55.604321]).reshape(-1,1), coastalinfluence_design.todense())