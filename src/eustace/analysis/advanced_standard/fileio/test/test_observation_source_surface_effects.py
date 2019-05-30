"""Testing the source surface covariate effects, for several different cases"""

import numpy
import tempfile
import unittest

from netCDF4 import Dataset

from eustace.analysis.advanced_standard.fileio.observation_source_surface_effects import insitu_land_covariate_effect
from eustace.analysis.advanced_standard.fileio.observation_source_surface_effects import global_satellite_effect
from eustace.preprocess.fileio.insitu_land_breakpoints import ObservationBreakPointSourceInsituLand

class TestSourceInsituLandEffets(unittest.TestCase):
  
    def setUp(self):
        """
        We will have the following setup:
        
        9 stations (locations)
        10 time measurements
        """
        self.times = numpy.array([23740, 24836, 27027, 38715, 39811, 47116, 47117, 51499, 53325, 54056])
        self.breakpoints_file = tempfile.NamedTemporaryFile(prefix='eustace.analysis.advanced_standard.elements.test_bias_insitu_land', suffix='.nc',delete=True)

        self.break_stations = numpy.array([4, 4, 4, 4, 6, 6, 7, 8 , 8 , 8], dtype=numpy.int32)
        self.break_times = numpy.array([47116, 51499, 53325, 54056, 38715, 39811, 47116, 23740, 24836, 27027], dtype=numpy.int32)
        self.break_likelihood = numpy.array([1, 8, 7, 10, 11, 2, 23, 4, 7, 10], dtype=numpy.int8)
        self.detection_feasibility = numpy.array(range(1,21), dtype=numpy.int8)
        dataset=Dataset(self.breakpoints_file.name,'w','NETCDF4')
        dataset.createDimension('merged_break',size=10)
        dataset.createDimension('station',size=len(self.detection_feasibility))
        merged_break_time=dataset.createVariable('merged_break_time',numpy.int32,dimensions=('merged_break'))
        merged_break_station=dataset.createVariable('merged_break_station',numpy.int32,dimensions=('merged_break'))
        merged_break_likelihood=dataset.createVariable('merged_break_likelihood',numpy.int8,dimensions=('merged_break'))
        detection_feasibility=dataset.createVariable('detection_feasibility',numpy.int8,dimensions=('station'))
        merged_break_time[:]=self.break_times
        merged_break_station[:]=self.break_stations
        merged_break_likelihood[:]=self.break_likelihood
        detection_feasibility[:]=self.detection_feasibility
        dataset.close()  
  
    def test_insitu_land_covariate_effect(self):
        reader = ObservationBreakPointSourceInsituLand(self.breakpoints_file.name)
        observed_breakpoints = reader.observations('merged_break')

        # Testing covariate for a certain number of mock cases
        # break points =  [47116, 51499, 53325, 54056, 38715, 39811, 47116, 23740, 24836, 27027]
        # stations     =  [4,     4,     4,     4,     6,     6,     7,     8 ,    8 ,    8]  # indexed from 1 in input file
        
        # measurement times = [27027, 39811, 47117]
        # valid observation indices = [0,3,5]                                                 # indexed from 0 in analysis
        # expected covariate effects = [[station index_i, breakpoints index_j]]
        # t = 27027 -> [[1, 0],[2, 4]]
        # t = 39811 -> [[1, 0],[2, 5]]
        # t = 47117 -> [[1, 1]]
 #       break_indices = [2, 4, 6]
 #       valid_indices = numpy.array([0, 3 ,5])
 #       expected_result = [numpy.array([[1, 0],[2, 4]]), numpy.array([[1, 0],[2,5]]), numpy.array([[1, 1]])]
 #       for index, result in zip(break_indices, expected_result):
 #           numpy.testing.assert_array_equal(insitu_land_covariate_effect(self.times[index], valid_indices, observed_breakpoints), result, err_msg='Testing example breakpoint '+str(self.times[index]))

        # measurement times = [38715, 51499, 54056]
        # valid observation indices = [3, 5, 6, 7]                                # indexed from 0 in analysis
        # expected covariate effects = [[station index_i, breakpoints index_j]]
        # t = 38715 -> [[0, 0],[1, 4],[2, 6]]
        # t = 51499 -> [[0, 1]]
        # t = 54056 -> [[0,3]]
        break_indices = [3, 7, 9]
        valid_indices = numpy.array([3, 5, 6, 7])
        expected_result = [numpy.array([[0, 0],[1, 4],[2, 6]]), numpy.array([[0, 1]]), numpy.array([[0, 3]])]
        for index, result in zip(break_indices, expected_result):
            numpy.testing.assert_array_equal(insitu_land_covariate_effect(self.times[index], valid_indices, observed_breakpoints), result, err_msg='Testing example breakpoint '+str(self.times[index]))

        # measurement times = [23740, 24836, 53325]
        # valid observation indices = [6, 7]                                      # indexed from 0 in analysis
        # expected covariate effects = [[station index_i, breakpoints index_j]]
        # t = 23740 -> [[0, 6],[1, 7]]
        # t = 24836 -> [[0, 6],[1, 8]]
        # t = 53325 -> None
        break_indices = [0, 1, 8]
        valid_indices = numpy.array([6, 7])
        expected_result = [numpy.array([[0, 6],[1, 7]]), numpy.array([[0, 6],[1, 8]]), None]
        for index, result in zip(break_indices, expected_result):
            numpy.testing.assert_array_equal(insitu_land_covariate_effect(self.times[index], valid_indices, observed_breakpoints), result, err_msg='Testing example breakpoint '+str(self.times[index]))

    def test_global_satellite_effect(self):
      
        numpy.testing.assert_array_equal(global_satellite_effect(numpy.array([0, 1, 2, 5])), numpy.array([[0, 0], [1, 0],[2, 0],[3, 0]]))
        numpy.testing.assert_array_equal(global_satellite_effect(numpy.array([0, 1, 5])), numpy.array([[0, 0],[1, 0],[2, 0]]))
