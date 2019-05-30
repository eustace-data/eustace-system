import unittest
from datetime import datetime
import os.path

from eustace.analysis.advanced_standard.examples.inputloader_rawbinary import ObservationSourceSingleDayRawBinary_ComputedMean
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from eustace.timeutils.epoch import days_since_epoch
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystemInputLoader
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceBinaryFilenameGenerator
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceSingleDayRawBinary
from eustace.analysis.fileio.observationsource_rawbinary import LocalCorrelationRangeRawBinaryReader
from eustace.analysis.fileio.observationsource_rawbinary import LocationLookupRawBinaryReader
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceBinaryFilenameGenerator
from eustace.analysis.fileio.observationsource_rawbinary import ObservationSourceSingleDayRawBinary
from eustace.analysis.fileio.observationsource_rawbinary import ObservableFileSpec
from eustace.analysis.advanced_standard.examples.inputloader_rawbinary import AnalysisSystemInputLoaderRawBinary

class TestInputLoaderRawBinary(unittest.TestCase):

    BASEPATH = '/work/scratch/eustace/rawbinary'

    def test_computed_mean(self):
    
        source = 'insitu_land'
        time_index = int(days_since_epoch(datetime(2006, 02, 01)))
        location_lookup = LocationLookupRawBinaryReader().read( os.path.join(TestInputLoaderRawBinary.BASEPATH, 'locationlookup_insitu_land.bin') )
        namebuilder = ObservationSourceBinaryFilenameGenerator(source, TestInputLoaderRawBinary.BASEPATH)
        filespecs = {
            ObservationSource.TMIN: ObservableFileSpec(
                namebuilder.filename_observations(ObservationSource.TMIN, time_index),
                LocalCorrelationRangeRawBinaryReader().read(namebuilder.filename_local_correlation_ranges(ObservationSource.TMIN))),
            ObservationSource.TMAX: ObservableFileSpec(
                namebuilder.filename_observations(ObservationSource.TMAX, time_index),
                LocalCorrelationRangeRawBinaryReader().read(namebuilder.filename_local_correlation_ranges(ObservationSource.TMIN))) }
        result = ObservationSourceSingleDayRawBinary_ComputedMean(location_lookup, filespecs, time_index)

