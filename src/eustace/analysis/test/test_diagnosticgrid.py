"""Tests for diagnostic gridding."""

import unittest
import numpy
from eustace.satgrid_iris.gridbin.gridfieldlist import GridFieldAxis
from ..observationsource import ObservationSource
from ..observationsource import Observations
from ..diagnosticgrid import DiagnosticGridInput
from ..diagnosticgrid import DiagnosticGridObservationConnector
from ..diagnosticgrid import DiagnosticGridBins

class SimulatedObservationSource(ObservationSource):
    """Helper class providing simulated observation test data."""

    def number_of_observations(self):
        """Total number of observations."""
        return 4

    def observables(self):
        """The names of variables estimated from this source."""
        raise [ 'pretend' ]

    def observation_location_lookup(self):
        """NumPy array of of latitude, longitude.
           Note extra unused location in middle to check this works."""

        return numpy.array( [ [  45.0, 38.0, 42.0, 58.0, 41.0 ],
                              [  -8.9, 22.1,  0.0, -5.0, 10.0 ] ], numpy.float32)

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""

        if observable == 'pretend':
            mask = numpy.array([ False, False, False, False ], numpy.bool)
            time = numpy.array([ 57500, 57499, 57500, 57500 ], numpy.float32)
            location = numpy.array([ 0, 1, 3, 4 ], numpy.uint64)
            measurements = numpy.array([ 275.0, 288.0, 204.0, 310.0 ], numpy.float32)
            uncorrelatederror = numpy.array([5.0, 5.0, 7.0, 8.0], numpy.float32)
            locallycorrelatederror = numpy.array([ [ 2.3, 3.3,  1.2, 8.9 ], [ 0.222, 7.6, 8.888, 9.999 ] ], numpy.float32)
            return Observations(mask, time, location, measurements, uncorrelatederror, locallycorrelatederror)

    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty (two components in this case)."""

        if observable == 'pretend':
            return numpy.array([ 6.6, 5.5 ], dtype=numpy.float64)

class TestDiagnosticGridInput(unittest.TestCase):

    def test_init(self):
        s = DiagnosticGridInput(
            'madeup',
            133077, 
            numpy.array([ 9, 0, 23, 22], numpy.int32),
            numpy.array([100.0, 500.0, 222.33, 89.99], numpy.float64),
            numpy.array([0.22, 0.33, 0.111, 0.0004], numpy.float64))
        self.assertEqual('madeup', s.observable)
        self.assertEqual(133077, s.daynumber)
        numpy.testing.assert_equal([ 9, 0, 23, 22], s.indices)
        numpy.testing.assert_equal([100.0, 500.0, 222.33, 89.99], s.observation)
        numpy.testing.assert_equal([0.22, 0.33, 0.111, 0.0004], s.uncertainty)


class TestDiagnosticGridObservationConnector(unittest.TestCase):

    def test_init(self):

        # Simulated axes
        axes = [ GridFieldAxis( 32.5, 5.0, 5, 'latitude', circular=False), GridFieldAxis(-15.0, 10.0, 4, 'longitude', circular=False) ]

        # Simulated data
        source = SimulatedObservationSource()

        # Make connector
        result = DiagnosticGridObservationConnector(axes, source)

        # Source and axes variables should be stored as-is
        self.assertEqual(axes, result.axes)
        self.assertEqual(source, result.source)

        # Check mapping to indices
        # Grid box start points are: [ 35, 40, 45, 50, 55 ] and [ -10, 0, 10, 20 ]
        # - should get  ( 2 ), ( 0 ), ( 1 ), ( 4 ), ( 1 )
        #               ( 0 ), ( 3 ), ( 1 ), ( 0 ), ( 2 )
        #
        # which means ( 2*4 + 0 ), (0*4 + 3), (1*4 + 1), (4*4 + 0), (1*4 + 2)
        # 
        # hence:          8      ,     3    ,      5   ,     16   ,     6
        #
        self.assertEqual(numpy.int32, result.grid_indices_for_location_ids.dtype)
        numpy.testing.assert_equal([8, 3, 5, 16, 6], result.grid_indices_for_location_ids)


    def test_get_day(self):

        # Simulated axes
        axes = [ GridFieldAxis( 32.5, 5.0, 5, 'latitude', circular=False), GridFieldAxis(-15.0, 10.0, 4, 'longitude', circular=False) ]

        # Make connector
        connector = DiagnosticGridObservationConnector(axes, SimulatedObservationSource())

        # Process a couple of days
        resultA = connector.get_day('pretend', 57500)
        resultB = connector.get_day('pretend', 57501)

        # Should have non-empty results for three items here
        self.assertEqual('pretend', resultA.observable)
        self.assertEqual(57500, resultA.daynumber)
        numpy.testing.assert_equal([8, 16, 6], resultA.indices)
        numpy.testing.assert_equal([275, 204, 310], resultA.observation)
        numpy.testing.assert_equal([5, 7, 8], resultA.uncertainty)
        
        # But this should be empty
        self.assertEqual('pretend', resultB.observable)
        self.assertEqual(57501, resultB.daynumber)
        self.assertEqual((0,), resultB.indices.shape)
        self.assertEqual((0,), resultB.observation.shape)
        self.assertEqual((0,), resultB.uncertainty.shape)


class TestDiagnosticGridBins(unittest.TestCase):
    """Tests for gridding of this data."""

    def test_create_fields_from_sparse_observations(self):
        
        # Axes to grid onto
        axes = [ GridFieldAxis( 32.5, 5.0, 5, 'latitude', circular=False), GridFieldAxis(-15.0, 10.0, 4, 'longitude', circular=False) ]

        # Simulated data
        connector = DiagnosticGridObservationConnector(axes, SimulatedObservationSource())
        dayA = connector.get_day('pretend', 57500)

        # The grid
        grid = DiagnosticGridBins(axes)

        # Aggregate the first data source
        grid.create_fields_from_sparse_observations('justtesting', dayA)

        # Retrieve results
        pretend_observation_weightedsum = grid.get_field('justtesting_pretend_observation_weightedsum')
        pretend_precision_sum = grid.get_field('justtesting_pretend_precision_sum')
        pretend_min = grid.get_field('justtesting_pretend_observation_min')
        pretend_max = grid.get_field('justtesting_pretend_observation_max')
        pretend_count = grid.get_field('justtesting_pretend_count')
        pretend_uncertainty_max = grid.get_field('justtesting_pretend_uncertainty_max')
        pretend_uncertainty_min = grid.get_field('justtesting_pretend_uncertainty_min')

        # Check we have all results as appropriate type of numpy array
        self.assertEqual(numpy.float64,  pretend_observation_weightedsum.dtype)
        self.assertEqual(numpy.float64,  pretend_precision_sum.dtype)
        self.assertTrue(numpy.float64, pretend_min.dtype)
        self.assertTrue(numpy.float64, pretend_max.dtype)
        self.assertTrue(numpy.int64, pretend_count.dtype)
        self.assertTrue(numpy.float64, pretend_uncertainty_min.dtype)
        self.assertTrue(numpy.float64, pretend_uncertainty_max.dtype)


        # Check result values

        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_count.data.mask)
        numpy.testing.assert_equal([ [ 0 , 0 , 0 , 0 ],
                                     [ 0 , 0 , 1 , 0 ],
                                     [ 1 , 0 , 0 , 0 ],
                                     [ 0 , 0 , 0 , 0 ],
                                     [ 1 , 0 , 0 , 0 ] ],
                                   pretend_count.data.data)

        self.assertEqual((5,4,), pretend_observation_weightedsum.data.shape)
        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_observation_weightedsum.data.mask)
        self.assertAlmostEqual(310.0/8.0, pretend_observation_weightedsum.data.data[1,2])
        self.assertAlmostEqual(275.0/5.0, pretend_observation_weightedsum.data.data[2,0])
        self.assertAlmostEqual(204.0/7.0, pretend_observation_weightedsum.data.data[4,0])

        self.assertEqual((5,4,), pretend_precision_sum.data.shape)
        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_precision_sum.data.mask)
        self.assertAlmostEqual(1.0/8.0, pretend_precision_sum.data.data[1,2])
        self.assertAlmostEqual(1.0/5.0, pretend_precision_sum.data.data[2,0])
        self.assertAlmostEqual(1.0/7.0, pretend_precision_sum.data.data[4,0])

        self.assertEqual((5,4,), pretend_max.data.shape)
        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_max.data.mask)
        self.assertAlmostEqual(310.0, pretend_max.data.data[1,2])
        self.assertAlmostEqual(275.0, pretend_max.data.data[2,0])
        self.assertAlmostEqual(204.0, pretend_max.data.data[4,0])

        self.assertEqual((5,4,), pretend_min.data.shape)
        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_min.data.mask)
        self.assertAlmostEqual(310.0, pretend_min.data.data[1,2])
        self.assertAlmostEqual(275.0, pretend_min.data.data[2,0])
        self.assertAlmostEqual(204.0, pretend_min.data.data[4,0])

        self.assertEqual((5,4,), pretend_uncertainty_max.data.shape)
        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_uncertainty_max.data.mask)
        self.assertAlmostEqual(8.0, pretend_uncertainty_max.data.data[1,2])
        self.assertAlmostEqual(5.0, pretend_uncertainty_max.data.data[2,0])
        self.assertAlmostEqual(7.0, pretend_uncertainty_max.data.data[4,0])

        self.assertEqual((5,4,), pretend_uncertainty_min.data.shape)
        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_uncertainty_min.data.mask)
        self.assertAlmostEqual(8.0, pretend_uncertainty_min.data.data[1,2])
        self.assertAlmostEqual(5.0, pretend_uncertainty_min.data.data[2,0])
        self.assertAlmostEqual(7.0, pretend_uncertainty_min.data.data[4,0])

    def test_compute_weighted_mean(self):

        
        # Axes to grid onto
        axes = [ GridFieldAxis( 32.5, 5.0, 5, 'latitude', circular=False), GridFieldAxis(-15.0, 10.0, 4, 'longitude', circular=False) ]

        # Simulated data
        connector = DiagnosticGridObservationConnector(axes, SimulatedObservationSource())
        dayA = connector.get_day('pretend', 57500)

        # The grid
        grid = DiagnosticGridBins(axes)

        # Aggregate from data source
        grid.create_fields_from_sparse_observations('justtesting', dayA)

        # Compute mean
        grid.compute_weighted_mean('justtesting', 'pretend')

        # Check results
        pretend_observation_weightedmean = grid.get_field('justtesting_pretend_observation_weightedmean')
        
        numpy.testing.assert_equal([ [ True , True , True , True ],
                                     [ True , True , False, True ],
                                     [ False, True , True , True ],
                                     [ True , True , True , True ],
                                     [ False, True , True , True ] ],
                                   pretend_observation_weightedmean.data.mask)

        self.assertAlmostEqual(310.0, pretend_observation_weightedmean.data.data[1,2])
        self.assertAlmostEqual(275.0, pretend_observation_weightedmean.data.data[2,0])
        self.assertAlmostEqual(204.0, pretend_observation_weightedmean.data.data[4,0])
