
import unittest
import numpy
import scipy.sparse
from StringIO import StringIO

from eustace.timeutils.epoch import epoch_plus_days

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystemInputLoader

from eustace.analysis.advanced_standard.fileio.output_projector import Projector
from eustace.analysis.advanced_standard.elements.element import ObservationStructure
from eustace.analysis.observationsource import ObservationSource

from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory

class SimulatedObservationStructure(ObservationStructure):

    def __init__(self, time_index):

        # Store time index
        self.t = time_index

        # Inspect the same mesh as used for output
        # mesh = MeshIcosahedronSubdivision.build(number_of_subdivisions=1)
        # print cartesian_to_polar2d(mesh.points)
        #
        # Index 12: [   0,    0  ]
        # Index 17: [   0,  180  ]
        # Index 41: [   0,   -90 ]
        
        # Simulated locations
        self.locations = numpy.array( [ [ 0.0,   0.0 ],
                                        [ 0.0, 180.0 ],
                                        [ 0.0, -90.0 ] ] )

        # Simulated measurements
        self.measurement = numpy.array( [ 20.0, -15.0, 5.0 ] )

        # Simulated errors
        self.uncorrelatederror = 0.1*numpy.ones(self.measurement.shape)


    def time_index(self):

        return self.t
        
    def time_datetime(self):

        return epoch_plus_days(self.t)

    def number_of_observations(self):

        return self.locations.shape[0]
        
    def location_polar_coordinates(self):

        return self.locations
                
    def observation_vector(self):
        """Array of observations."""

        return self.measurement

    def observation_precision(self):
        """Observation precision matrix (sparse)."""

        return scipy.sparse.diags(1.0 / (self.uncorrelatederror ** 2), format='csc')

class SimulatedInputLoader(AnalysisSystemInputLoader):
    """Simulated loading of input data."""

    def __init__(self):

        pass

    def datetime_at_time_index(self, time_index):
        """Use EUSTACE epoch."""

        return epoch_plus_days(time_index)
        
    def load_observation_structure(self, observable, time_index, log=None):
        """Load simulated data at specified time index."""

        if observable == ObservationSource.TMEAN:

            return SimulatedObservationStructure(time_index)

def mini_world_spatial():

    # Local component
    local_component = SpatialComponent(ComponentStorage_InMemory(LocalElement(n_triangulation_divisions=4),
                                        LocalHyperparameters(log_sigma=numpy.log(5.0), log_rho=numpy.log( 15. * numpy.pi / 180. )  )),
                                        SpatialComponentSolutionStorage_InMemory(), compute_uncertainties=False, method='APPROXIMATED', compute_sample = True, sample_size=100)

    # Analysis system using the specified components, for the Tmean observable
    analysis_system = AnalysisSystem([ local_component ], ObservationSource.TMEAN, log=StringIO())
    
    # Simulated inputs
    simulated_input_loader = SimulatedInputLoader()

    # Simulate evaluation of this time index
    simulated_time_indices = [ 0 ]

    # Update with data
    analysis_system.update([ simulated_input_loader ], simulated_time_indices)

    # Check state vector directly
    statevector = analysis_system.components[0].solutionstorage.partial_state_read(0)

    return analysis_system, statevector


class TestProjector(unittest.TestCase):
  
    def setUp(self):
        
        analysis_system, statevector = mini_world_spatial()
        self.analysis_system = analysis_system
        self.statevector = statevector
    
        centre_latitudes = numpy.array([-1.0, 0.0, 1.0])
        centre_longitudes = numpy.array([-1.0, 0.0, 1.0])
        
        delta_lat = 2.0
        delta_lon = 2.0
        centre_latitudes = numpy.linspace( -90.+delta_lat, 90.-delta_lat, int(180./delta_lat) )
        centre_longitudes = numpy.linspace( -180.+delta_lon, 180.-delta_lon, int(360./delta_lon) )
        
        grid_resolution = [delta_lat, delta_lon]
        cell_sampling = [3, 3]
        blocking = 1
        time_index = 0
        self.projector = Projector(centre_latitudes, centre_longitudes, grid_resolution, time_index, cell_sampling, blocking)
    
        self.projector.set_component(self.analysis_system.components[0])

    def tearDown(self):
        self.analysis_system = None
    
    def test_read_state(self):

        my_state = self.projector.read_state_solution(self.analysis_system.components[0].solutionstorage)
        numpy.testing.assert_almost_equal(my_state, self.statevector)
    
    def atest_design_matrix(self):
        self.projector.evaluate_design_matrix()
        print self.projector.design_matrix.shape
        
        expected_values =  self.projector.project_expected_value()
        uncertainty_values =  self.projector.project_sample_std()
        uncertainty_values2 =  self.projector.project_sample_deviation()
        import matplotlib.pyplot as plt
        
        plt.figure()
        plt.imshow( expected_values.reshape( (len(self.projector.centre_latitudes), len(self.projector.centre_longitudes) ) ) )
        
        plt.figure()
        plt.imshow( uncertainty_values.reshape( (len(self.projector.centre_latitudes), len(self.projector.centre_longitudes) ) ) )
        
        plt.figure()
        plt.plot( uncertainty_values.ravel() )
        plt.plot( uncertainty_values2.ravel() )
        
        plt.figure()
        plt.imshow( uncertainty_values.reshape( (len(self.projector.centre_latitudes), len(self.projector.centre_longitudes) ) ) )

        plt.show()
        
    
    def test_project_expectation(self):
        # Check simple projection of state through identity matrix
        self.projector.design_matrix = scipy.sparse.eye( len(self.statevector.ravel()) ).tocsr()

        expected_values = self.projector.project_expected_value()
        numpy.testing.assert_almost_equal(expected_values, self.statevector)
    
    def test_project_sample(self):
        # Check simple projection of sample variates through identity matrix
        self.projector.design_matrix = scipy.sparse.eye( len(self.statevector.ravel()) ).tocsr()
        
        sample_values = self.projector.project_sample_values()
        numpy.testing.assert_almost_equal(sample_values, self.projector.state_samples )
    
    def test_project_deviation(self):
        # Check deviations of sample from mean for known deviation
        
        FIXED_EXPECTATION = 0.0
        FIXED_STD = 2.0
        self.projector.state_solution[:,:] = FIXED_EXPECTATION
        self.projector.state_samples[:,:] = FIXED_STD
        
        # test with model design matrix
        self.projector.evaluate_design_matrix()
        numpy.testing.assert_almost_equal(self.projector.project_sample_deviation(), numpy.ones((self.projector.design_matrix.shape[0],1)) * FIXED_STD )
        
        # test with identity design matrix
        self.projector.design_matrix = scipy.sparse.eye( len(self.statevector.ravel()) ).tocsr()
        numpy.testing.assert_almost_equal(self.projector.project_sample_deviation(), numpy.ones(self.projector.state_solution.shape) * FIXED_STD )
     
        
        
        