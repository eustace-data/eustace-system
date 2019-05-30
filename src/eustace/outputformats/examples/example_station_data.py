"""Generation of example files containing station data."""

__version__ = "$Revision: 562 $"
__author__ = "Joel R. Mitchelson"

from eustace.outputformats import definitions
from eustace.outputformats.stations import FileBuilderStationData
from eustace.outputformats.stationdata import StationDescriptionArray
from eustace.outputformats.stationdata import StationDescription
from eustace.outputformats.stationdata import StationData
from eustace.outputformats.stationdata import StationTemperatureSample
from eustace.outputformats.stationdata import StationBreakDetectionQuality
from eustace.outputformats.stationdata import StationBreak
import numpy
import time
import math


def write_example_station_data(source, version, outputdirectory, institution):
    """Make an example output format for station data."""

    # station descriptors for example data
    descriptions = StationDescriptionArray([
        StationDescription('Cholet', 47.059407, -0.879787, 85.0),
        StationDescription('Adelaide', -34.9286212, 138.5999594, 45.0),
        StationDescription('Tobermory', 56.6227813, -6.0723004, 53.0),
        StationDescription('Orinoco', 5.4329819, -64.9895156, 628.0),
        StationDescription('Bungo', 35.3078828, 139.2819831, 15.0),
        StationDescription('Wellington', 50.977071, -3.225954, 75.0),
        StationDescription('Tomsk', 56.5010397, 84.9924506, 144.0),
        StationDescription('Sofia', 42.660993, 23.333331, 598.0),
    ])

    # helper object to contain station data and put it onto a common time axis
    data = StationData(descriptions)

    # some sample times at which we will provide data
    # - measured in days since 01/01/1850 00:00 UTC
    simulated_times = numpy.arange(43342, 43442, 1, dtype=numpy.float32)

    # simulate missing data at one station
    wellington_missing_time_range = [43399.0, 43600.0]

    # simulate breaks at some stations
    sofia_break_time = 43360.0
    sofia_break_affected_time_range = [43350.0, 43370.0]
    sofia_break_amplitude = 3.5
    bungo_break_time = 43403.0
    bungo_break_affected_time_range = [43399.0, 43439.0]
    bungo_break_amplitude = 5.0

    # add some data
    # (example only, does not correspond in any way to actual temperatures)
    for t in simulated_times:    # pylint: disable=invalid-name

        # simulated breaks
        sofia_break_offset = 0.0 if t < sofia_break_time else sofia_break_amplitude
        bungo_break_offset = 0.0 if t < bungo_break_time else bungo_break_amplitude

        data.add_temperature_sample(StationTemperatureSample(
            station_name='Sofia',
            time=t,
            tasmin=273.15 + 5.0 *
            math.sin((t + 3) * 0.0086) + sofia_break_offset,
            tasmax=293.15 + 5.0 *
            math.sin((t + 10) * 0.0086) + sofia_break_offset,
            tasmin_qc=0,
            tasmax_qc=0))

        data.add_temperature_sample(StationTemperatureSample(
            station_name='Cholet',
            time=t,
            tasmin=285.15 + 4.0 * math.sin((t + 30) * 0.0086),
            tasmax=288.15 + 3.0 * math.sin((t + 32) * 0.0086),
            tasmin_qc=0,
            tasmax_qc=0))

        data.add_temperature_sample(StationTemperatureSample(
            station_name='Adelaide',
            time=t,
            tasmin=290.15 + 3.0 * math.sin((t + 60) * 0.0086),
            tasmax=299.15 + 3.0 * math.sin((t + 59) * 0.0086),
            tasmin_qc=0,
            tasmax_qc=0))

        data.add_temperature_sample(StationTemperatureSample(
            station_name='Tobermory',
            time=t,
            tasmin=270.15 + 5.0 * math.sin((t + 40) * 0.0086),
            tasmax=282.15 + 4.5 * math.sin((t + 39) * 0.0086),
            tasmin_qc=0,
            tasmax_qc=0))

        data.add_temperature_sample(StationTemperatureSample(
            station_name='Orinoco',
            time=t,
            tasmin=274.15 + 4.0 * math.sin((t + 1) * 0.0086),
            tasmax=288.15 + 3.0 * math.sin((t + 2) * 0.0086),
            tasmin_qc=0,
            tasmax_qc=0))

        data.add_temperature_sample(StationTemperatureSample(
            station_name='Tomsk',
            time=t,
            tasmin=265.15 + 8.0 * math.sin((t + 10) * 0.0086),
            tasmax=299.15 + 7.0 * math.sin((t + 11) * 0.0086),
            tasmin_qc=0,
            tasmax_qc=0))

        data.add_temperature_sample(StationTemperatureSample(
            station_name='Bungo',
            time=t,
            tasmin=280.15 + 4.5 * math.sin((t + 5) * 0.0086),
            tasmax=289.15 + 4.0 *
            math.sin((t + 10) * 0.0086) + bungo_break_offset,
            tasmin_qc=0,
            tasmax_qc=0))

        if (t < wellington_missing_time_range[0]) or (t >= wellington_missing_time_range[1]):
            data.add_temperature_sample(StationTemperatureSample(
                station_name='Wellington',
                time=t,
                tasmin=285.15 + 12.0 * math.sin((t + 20) * 0.0086),
                tasmax=300.15 + 12.0 * math.sin((t + 19) * 0.0086),
                tasmin_qc=0,
                tasmax_qc=0))

    # add descriptions of break detection periods
    data.add_break_detection_period(StationBreakDetectionQuality(
        43098.0,
        StationBreakDetectionQuality.RELIABLE,
        StationBreakDetectionQuality.RELIABLE))

    data.add_break_detection_period(StationBreakDetectionQuality(
        43464.0,
        StationBreakDetectionQuality.RELIABLE,
        StationBreakDetectionQuality.RELIABLE))

    # add breaks that may have been detected
    # At station Sofia both TASMIN and TASMAX were affected
    # At station Bungo only TASMAX was affected

    data.add_break(definitions.TASMIN,
                   StationBreak(
                       station_name='Sofia',
                       amplitude=sofia_break_amplitude,
                       time_bounds=[sofia_break_time -
                                    15.0, sofia_break_time + 15.0],
                       time_affected_bounds=sofia_break_affected_time_range))

    data.add_break(definitions.TASMAX,
                   StationBreak(
                       station_name='Sofia',
                       amplitude=sofia_break_amplitude,
                       time_bounds=[sofia_break_time -
                                    15.0, sofia_break_time + 15.0],
                       time_affected_bounds=sofia_break_affected_time_range))

    data.add_break(definitions.TASMAX,
                   StationBreak(
                       station_name='Bungo',
                       amplitude=bungo_break_amplitude,
                       time_bounds=[bungo_break_time -
                                    15.0, bungo_break_time + 15.0],
                       time_affected_bounds=bungo_break_affected_time_range))

    # could retrieve sample times like this if we didn't already have them
    # time_axis = data.get_unique_temperature_sample_times()

    # time axis known in advance
    time_axis = simulated_times.tolist()

    # object to build global field file at current time
    builder = FileBuilderStationData(
        outputdirectory,
        region='Example',
        source=source,
        version=version,
        institution=institution,
        comment='EUSTACE project example file format for station data',
        history='Created ' + time.strftime('%c'))

    # common axes for temperature variables which update daily for every
    # station
    builder.set_temperature_dimensions(
        descriptions, numpy.array(time_axis, numpy.float32))

    # add temperature variables
    builder.add_temperature_variable(
        definitions.TASMIN, data.get_temperature_samples(time_axis, definitions.TASMIN))
    builder.add_temperature_variable(
        definitions.TASMIN_QC, data.get_temperature_samples(time_axis, definitions.TASMIN_QC))
    builder.add_temperature_variable(
        definitions.TASMAX, data.get_temperature_samples(time_axis, definitions.TASMAX))
    builder.add_temperature_variable(
        definitions.TASMAX_QC, data.get_temperature_samples(time_axis, definitions.TASMAX_QC))

    # dimensions for detection statistics
    builder.add_status_dimension(definitions.DIMENSION_NAME_DETECTION_TIME)
    builder.add_status_dimension(definitions.DIMENSION_NAME_TASMIN_BREAK)
    builder.add_status_dimension(definitions.DIMENSION_NAME_TASMAX_BREAK)

    # add detection statistics information
    builder.add_status_variable(definitions.DIMENSION_NAME_DETECTION_TIME,
                                definitions.DETECTION_TIME, data.get_break_detection_quality_variable(definitions.DETECTION_TIME))
    builder.add_status_variable(definitions.DIMENSION_NAME_DETECTION_TIME,
                                definitions.TASMIN_DETECTION_QC, data.get_break_detection_quality_variable(definitions.TASMIN_DETECTION_QC))
    builder.add_status_variable(definitions.DIMENSION_NAME_DETECTION_TIME,
                                definitions.TASMAX_DETECTION_QC, data.get_break_detection_quality_variable(definitions.TASMAX_DETECTION_QC))

    # tasmin break information
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMIN_BREAK,
                                definitions.TASMIN_BREAK_STATION, data.get_break_variable(definitions.TASMIN_BREAK_STATION))
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMIN_BREAK,
                                definitions.TASMIN_BREAK_AMPLITUDE, data.get_break_variable(definitions.TASMIN_BREAK_AMPLITUDE))
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMIN_BREAK,
                                definitions.TASMIN_BREAK_TIME_BOUNDS, data.get_break_variable(definitions.TASMIN_BREAK_TIME_BOUNDS), bounds=True)
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMIN_BREAK,
                                definitions.TASMIN_BREAK_TIME_AFFECTED_BOUNDS, data.get_break_variable(definitions.TASMIN_BREAK_TIME_AFFECTED_BOUNDS), bounds=True)

    # tasmax break information
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMAX_BREAK,
                                definitions.TASMAX_BREAK_STATION, data.get_break_variable(definitions.TASMAX_BREAK_STATION))
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMAX_BREAK,
                                definitions.TASMAX_BREAK_AMPLITUDE, data.get_break_variable(definitions.TASMAX_BREAK_AMPLITUDE))
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMAX_BREAK,
                                definitions.TASMAX_BREAK_TIME_BOUNDS, data.get_break_variable(definitions.TASMAX_BREAK_TIME_BOUNDS), bounds=True)
    builder.add_status_variable(definitions.DIMENSION_NAME_TASMAX_BREAK,
                                definitions.TASMAX_BREAK_TIME_AFFECTED_BOUNDS, data.get_break_variable(definitions.TASMAX_BREAK_TIME_AFFECTED_BOUNDS), bounds=True)

    # store the result
    builder.save_and_close()
