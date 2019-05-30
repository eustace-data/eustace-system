"""Tests for classes that represent station data."""

# pylint: disable=missing-docstring,invalid-name

__version__ = "$Revision: 552 $"
__author__ = "Joel R. Mitchelson"

import unittest
from eustace.outputformats.stationdata import StationData
from eustace.outputformats.stationdata import StationDescription
from eustace.outputformats.stationdata import StationDescriptionArray
from eustace.outputformats.stationdata import StationBreakDetectionQuality
from eustace.outputformats.stationdata import StationBreak
from eustace.outputformats.stationdata import StationTemperatureSample
from eustace.outputformats.stationdata import MeasurementTime
from eustace.outputformats import definitions
import numpy

class TestStationData(unittest.TestCase):

    def test_description(self):
        d = StationDescription('Bob', 57.0, 3.2, 2000.0)
        self.assertEqual('Bob', d.name)
        self.assertEqual(57.0, d.latitude)
        self.assertEqual(3.2, d.longitude)
        self.assertEqual(2000.0, d.elevation)

    def test_description_array(self):
        a = StationDescriptionArray(
            [StationDescription('OneStation', 33.3, 22.2, 88.8),
             StationDescription('TwoStations', -18.2, -40.0, 50.0),
             StationDescription('UnderTheSea', 4.0, 4.0, -3000.0)])
        self.assertEqual(3, a.get_number_of_stations())
        self.assertEqual(None, a.get_station_index('Nowhere'))
        self.assertEqual(2, a.get_station_index('UnderTheSea'))
        self.assertEqual(0, a.get_station_index('OneStation'))

    def test_temperature_sample(self):
        t = StationTemperatureSample(
            "Bob", 60000.22, 289.2, 27, 301.0, 5999)
        self.assertEqual("Bob", t.station_name)
        self.assertEqual(60000.22, t.time)
        self.assertEqual(289.2, t.tasmin)
        self.assertEqual(27, t.tasmin_qc)
        self.assertEqual(301.0, t.tasmax)
        self.assertEqual(5999, t.tasmax_qc)

    def test_break_detection_quality(self):
        self.assertEqual(0, StationBreakDetectionQuality.NOT_POSSIBLE)
        self.assertEqual(
            1, StationBreakDetectionQuality.POSSIBLE_BUT_UNRELIABLE)
        self.assertEqual(2, StationBreakDetectionQuality.RELIABLE)
        b = StationBreakDetectionQuality(
            33.5,
            StationBreakDetectionQuality.RELIABLE,
            StationBreakDetectionQuality.POSSIBLE_BUT_UNRELIABLE)
        self.assertEqual(33.5, b.detection_time)
        self.assertEqual(2, b.tasmin_detection_qc)
        self.assertEqual(1, b.tasmax_detection_qc)

    def test_break(self):
        b = StationBreak('ActionStation', 5.5, [60001.0, 60003.0], [60000.0, 60365.0])
        self.assertEqual('ActionStation', b.get_station_name())
        self.assertEqual(5.5, b.amplitude)
        self.assertEqual(2, len(b.time_bounds))
        self.assertEqual(2, len(b.time_affected_bounds))
        self.assertEqual(60001.0, b.time_bounds[0])
        self.assertEqual(60003.0, b.time_bounds[1])
        self.assertEqual(60000.0, b.time_affected_bounds[0])
        self.assertEqual(60365.0, b.time_affected_bounds[1])

    def test_data_init(self):
        sd = StationData([StationDescription('OneStation', 33.3, 22.2, 88.8),
                                   StationDescription('TwoStations', -18.2, -40.0, 50.0),
                                   StationDescription('UnderTheSea', 4.0, 4.0, -3000.0)])
        self.assertEqual(3, len(sd.descriptions))
        self.assertEqual('TwoStations', sd.descriptions[1].name)
        self.assertEqual(0, len(sd.temperature_samples))
        self.assertEqual(0, len(sd.break_detection_quality))
        self.assertEqual(0, len(sd.breaks['tasmin']))
        self.assertEqual(0, len(sd.breaks['tasmax']))

    def test_data_add_temperature_sample(self):
        sd = StationData([StationDescription('JustOne', 18, 78.0, 10.0)])
        sd.add_temperature_sample(StationTemperatureSample('JustOne', 60000.0, 272.89, 22, 290.22, 43))
        self.assertEqual(1, len(sd.temperature_samples))
        self.assertEqual('JustOne', sd.temperature_samples[0].station_name)
        self.assertEqual(60000.0, sd.temperature_samples[0].time)
        self.assertEqual(272.89, sd.temperature_samples[0].tasmin)
        self.assertEqual(22, sd.temperature_samples[0].tasmin_qc)
        self.assertEqual(290.22, sd.temperature_samples[0].tasmax)
        self.assertEqual(43, sd.temperature_samples[0].tasmax_qc)

    def test_data_add_break_detection_period(self):
        sd = StationData(
            [StationDescription('JustOne', 18, 78.0, 10.0)])
        self.assertEqual(0, len(sd.break_detection_quality))
        sd.add_break_detection_period(StationBreakDetectionQuality(
            60591.5, StationBreakDetectionQuality.RELIABLE, StationBreakDetectionQuality.POSSIBLE_BUT_UNRELIABLE))
        sd.add_break_detection_period(StationBreakDetectionQuality(
            60593.5, StationBreakDetectionQuality.NOT_POSSIBLE, StationBreakDetectionQuality.POSSIBLE_BUT_UNRELIABLE))
        sd.add_break_detection_period(StationBreakDetectionQuality(
            60599.0, StationBreakDetectionQuality.RELIABLE, StationBreakDetectionQuality.RELIABLE))
        self.assertEqual(3, len(sd.break_detection_quality))
        self.assertEqual(60591.5, sd.break_detection_quality[0].detection_time)
        self.assertEqual(2, sd.break_detection_quality[0].tasmin_detection_qc)
        self.assertEqual(1, sd.break_detection_quality[0].tasmax_detection_qc)
        self.assertEqual(60593.5, sd.break_detection_quality[1].detection_time)
        self.assertEqual(0, sd.break_detection_quality[1].tasmin_detection_qc)
        self.assertEqual(1, sd.break_detection_quality[1].tasmax_detection_qc)
        self.assertEqual(60599.0, sd.break_detection_quality[2].detection_time)
        self.assertEqual(2, sd.break_detection_quality[2].tasmin_detection_qc)
        self.assertEqual(2, sd.break_detection_quality[2].tasmax_detection_qc)

    def test_data_add_break(self):
        sd = StationData(
            [StationDescription('JustOne', 18, 78.0, 10.0)])
        sd.add_break(definitions.TASMAX, StationBreak('JustOne', 3.22, [60466.0, 60511.0], [60100.0, 60800.0]))
        self.assertEqual(0, len(sd.breaks['tasmin']))
        self.assertEqual(1, len(sd.breaks['tasmax']))
        self.assertEqual('JustOne', sd.breaks['tasmax'][0].station_name)
        self.assertEqual(3.22, sd.breaks['tasmax'][0].amplitude)
        self.assertEqual(60800.0, sd.breaks['tasmax'][0].time_affected_bounds[1])

    def test_get_unique_temperature_sample_times(self):

        sd = StationData([
            StationDescription('StationOne', 18.0, 78.0, 10.0),
            StationDescription('StationTwo', 53.0, 22.0, 500.0)])

        sd.add_temperature_sample(StationTemperatureSample('StationOne', 60000.0001, 288.15, 0, 288.15, 0))
        sd.add_temperature_sample(StationTemperatureSample('StationOne', 60003.0000, 288.15, 0, 288.15, 0))
        sd.add_temperature_sample(StationTemperatureSample('StationOne', 60001.0000, 288.15, 0, 288.15, 0))
        sd.add_temperature_sample(StationTemperatureSample('StationTwo', 60003.0000, 288.15, 0, 288.15, 0))
        sd.add_temperature_sample(StationTemperatureSample('Stationtwo', 60000.0002, 288.15, 0, 288.15, 0))

        result = sd.get_unique_temperature_sample_times()

        self.assertEqual([60000.0001, 60001.0000, 60003.0000], result)

    def test_get_temperature_samples(self):

        sd = StationData(StationDescriptionArray([
            StationDescription('StationOne', 18.0, 78.0, 10.0),
            StationDescription('StationTwo', 53.0, 22.0, 500.0)]))

        sd.add_temperature_sample(StationTemperatureSample('StationOne', 60000.0001, 288.15, 11, 298.15, 20))
        sd.add_temperature_sample(StationTemperatureSample('StationOne', 60003.0000, 287.15, 13, 297.15, 0))
        sd.add_temperature_sample(StationTemperatureSample('StationOne', 60001.0000, 286.15, 15, 296.15, 23))
        sd.add_temperature_sample(StationTemperatureSample('StationTwo', 60003.0000, 299.15, 17, 301.05, 18))
        sd.add_temperature_sample(StationTemperatureSample('StationTwo', 60000.0002, 300.15, 53, 302.05, 1))

        tasmax = sd.get_temperature_samples([60000.0, 60001.0, 60003.0], definitions.TASMAX)
        tasmin_qc = sd.get_temperature_samples([60000.0, 60001.0, 60003.0], definitions.TASMIN_QC)

        numpy.testing.assert_array_almost_equal([[298.15, 302.05], [296.15, 0.0], [297.15, 301.05]], tasmax, decimal=4)

        numpy.testing.assert_array_equal([[11, 53], [15, 255], [13, 17]], tasmin_qc)

    def test_get_break_detection_quality_variable(self):

        sd = StationData(StationDescriptionArray([]))

        sd.add_break_detection_period(StationBreakDetectionQuality(
            60591.5, StationBreakDetectionQuality.RELIABLE, StationBreakDetectionQuality.POSSIBLE_BUT_UNRELIABLE))
        sd.add_break_detection_period(StationBreakDetectionQuality(
            60593.5, StationBreakDetectionQuality.NOT_POSSIBLE, StationBreakDetectionQuality.POSSIBLE_BUT_UNRELIABLE))
        sd.add_break_detection_period(StationBreakDetectionQuality(
            60599.0, StationBreakDetectionQuality.RELIABLE, StationBreakDetectionQuality.RELIABLE))

        numpy.testing.assert_array_equal(
            numpy.array([60591.5, 60593.5, 60599.0], numpy.float32),
            sd.get_break_detection_quality_variable(definitions.DETECTION_TIME))

        numpy.testing.assert_array_equal(
            numpy.array([2, 0, 2], numpy.int32),
            sd.get_break_detection_quality_variable(definitions.TASMIN_DETECTION_QC))

        numpy.testing.assert_array_equal(
            numpy.array([1, 1, 2], numpy.int32),
            sd.get_break_detection_quality_variable(definitions.TASMAX_DETECTION_QC))

    def test_get_break_variable(self):

        sd = StationData(StationDescriptionArray([
            StationDescription('StationOne', 18.0, 78.0, 10.0),
            StationDescription('StationTwo', 53.0, 22.0, 500.0)]))

        sd.add_break(definitions.TASMAX, StationBreak('StationTwo', 3.11, [60466.5, 60511.5], [60100.5, 60800.5]))
        sd.add_break(definitions.TASMAX, StationBreak('StationOne', 3.33, [60511.5, 60599.5], [60400.5, 60700.5]))
        sd.add_break(definitions.TASMAX, StationBreak('StationTwo', 3.34, [60599.5, 61000.5], [60700.0, 61100.0]))

        numpy.testing.assert_array_equal(
            numpy.array([1, 0, 1], numpy.int32),
            sd.get_break_variable(definitions.TASMAX_BREAK_STATION))

        numpy.testing.assert_array_equal(
            numpy.array([3.11, 3.33, 3.34], numpy.float32),
            sd.get_break_variable(definitions.TASMAX_BREAK_AMPLITUDE))

        numpy.testing.assert_array_equal(
            numpy.array([[60466.5, 60511.5], [60511.5, 60599.5],
                         [60599.5, 61000.5]], numpy.float32),
            sd.get_break_variable(definitions.TASMAX_BREAK_TIME_BOUNDS))

        numpy.testing.assert_array_equal(
            numpy.array([[60100.5, 60800.5], [60400.5, 60700.5],
                         [60700.0, 61100.0]], numpy.float32),
            sd.get_break_variable(definitions.TASMAX_BREAK_TIME_AFFECTED_BOUNDS))

        numpy.testing.assert_array_equal(
            numpy.array([], numpy.int32),
            sd.get_break_variable(definitions.TASMIN_BREAK_STATION))

        numpy.testing.assert_array_equal(
            numpy.array([], numpy.float32),
            sd.get_break_variable(definitions.TASMIN_BREAK_AMPLITUDE))

        numpy.testing.assert_array_equal(
            numpy.array([], numpy.float32),
            sd.get_break_variable(definitions.TASMIN_BREAK_TIME_BOUNDS))

        numpy.testing.assert_array_equal(
            numpy.array([], numpy.float32),
            sd.get_break_variable(definitions.TASMIN_BREAK_TIME_AFFECTED_BOUNDS))


class TestMeasurementTime(unittest.TestCase):

    def test_init(self):
        m = MeasurementTime(239.22)
        self.assertEqual(239.22, m.time)

    def test_equal(self):
        self.assertTrue(MeasurementTime(290.0000) == MeasurementTime(289.9996))
        self.assertFalse(MeasurementTime(290.0000) == MeasurementTime(289.9994))

    def test_set(self):

        s = set([MeasurementTime(t) for t in [60000.0001, 60003.0, 60001.0, 60003.0, 60000.0002]])
        results = [m.time for m in s]

        self.assertEqual(3, len(results))
        self.assertTrue(60000.0001 in results)
        self.assertTrue(60003.0000 in results)
        self.assertTrue(60001.0000 in results)
