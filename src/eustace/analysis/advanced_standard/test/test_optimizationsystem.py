"""Test of optimizationsysem using simulated data, low-resolution mesh, and dummy climatology and largescale models."""


import numpy
import scipy.sparse
import unittest

from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem
from eustace.analysis.advanced_standard.optimizationsystem import OptimizationSystem
from eustace.analysis.advanced_standard.examples import example_eustace
from eustace.analysis.advanced_standard.examples import example_optimization
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory, SpaceTimeComponentSolutionStorage_InMemory, SpatialComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.spatialdelayed import DelayedSpatialComponent
from eustace.analysis.advanced_standard.elements.combination import CombinationElement, CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanElement

from eumopps.catalogue.fileio.test.test_formatnetcdf import NamedTemporaryDirectory
import os.path
import tempfile
import shutil
import uuid 

from eustace.analysis.observationsource import ObservationSource, Observations
from eustace.analysis.fileio.observationsource_rawbinary import LocationLookupRawBinaryWriter, LocationLookupWithID, LocalCorrelationRangeRawBinaryWriter, ObservationSourceBinaryFilenameGenerator, ObservationRawBinaryWriter

import datetime
from eustace.timeutils import epoch

class TestOptimizationSystem(unittest.TestCase):
    
    def setUp(self):
        
        self.setup_path()
        self.generate_test_rawbinaries()
        
        # dummy climatology model
        climatology_element = CombinationElement( [GrandMeanElement(),])
        climatology_hyperparameters = CombinationHyperparameters( [CovariateHyperparameters(numpy.log(15.0)),] )
        climatology_definition = ComponentStorage_InMemory(climatology_element, climatology_hyperparameters)
        climatology_storage = SpaceTimeComponentSolutionStorage_InMemory()
        
        # dummy large scale model
        largescale_element = CombinationElement( [GrandMeanElement(),])
        largescale_hyperparameters = CombinationHyperparameters( [CovariateHyperparameters(numpy.log(15.0)),] )
        largescale_definition = ComponentStorage_InMemory(largescale_element, largescale_hyperparameters)
        largescale_storage = SpaceTimeComponentSolutionStorage_InMemory()
        
        # local spde+bias model for analysissystem
        global_biases, global_biases_group_list = False, []
        local_definition = example_eustace.LocalDefinition(global_biases, global_biases_group_list)
        local_storage = SpatialComponentSolutionStorage_InMemory()
        
        # bias component for split local components of optimizationsystem
        local_bias_definition = example_optimization.PureBiasComponentDefinition(global_biases, global_biases_group_list)
        local_bias_storage = SpatialComponentSolutionStorage_InMemory()
        
        # local spde component for split local components of optimizationsystem
        local_spde_definition = example_optimization.PureLocalComponentDefinition()
        local_spde_storage = SpatialComponentSolutionStorage_InMemory()
        
        neighbourhood_level = 0
        region_index = 0
        region_definition = example_optimization.LocalViewDefinition( neighbourhood_level, region_index, regionspec = 'LocalSubRegion' )
        region_storage = SpatialComponentSolutionStorage_InMemory()
    
        # analysis sytem init
        self.analysissystem = AnalysisSystem( [ SpaceTimeComponent(climatology_definition, climatology_storage),
                                                SpaceTimeComponent(largescale_definition, largescale_storage),
                                                DelayedSpatialComponent(local_definition, local_storage) ],
                                                observable = ObservationSource.TMEAN )
        
        # analysis sytem init
        self.optimizationsystem = OptimizationSystem( [ SpaceTimeComponent(climatology_definition, climatology_storage),
                                                        SpaceTimeComponent(largescale_definition, largescale_storage),
                                                        DelayedSpatialComponent(local_bias_definition, local_bias_storage),
                                                        DelayedSpatialComponent(region_definition, region_storage), ],
                                                        observable = ObservationSource.TMEAN )

        # need to setup state reading from disk for full model - does using the standard setup do this? Maybe not for local - separate jobs per day leads to missing initialisation of files.
        
        # need to test state split into regions for optimisationsystem
        
        # need to test merge of hyperparameters for system
        
        # need to test use of non-stationary model in system

    def setup_path(self):
        self.pathname = tempfile.mkdtemp()
        
    def tearDown(self):
        
        del self.analysissystem
        del self.optimizationsystem
        
        # delete temporary files
        shutil.rmtree(self.pathname)
        
    def test_process_inputs(self):
        # test loading multiple time steps from json descriptor of dummy input data
        
        print self.analysissystem
        print self.optimizationsystem
        
        self.setup_input_descriptors()
        component_index = 3
        #self.optimizationsystem.process_inputs(self.input_descriptors, component_index, ["20170601","20170602"])
        self.optimizationsystem.update_component(self.input_descriptors, component_index, ["20170601","20170602"])
    
        time_index = int( epoch.days_since_epoch( datetime.datetime(2017, 6, 1) ) )
        print self.optimizationsystem.components[component_index].solutionstorage.measurement_at_time
        print self.optimizationsystem.components[component_index].solutionstorage.state_at_time
        #print self.optimizationsystem.components[2].solutionstorage.partial_state_read(time_index)
        
        hyperparameter_file = os.path.join( self.pathname, "test_hyperparameters.npy")
        
        raise NotImplemented("Optimisation system tests requires small dummy system to optimize")
        #self.optimizationsystem.optimize_component(component_index, hyperparameter_storage_file=hyperparameter_file)
    

    def generate_test_rawbinaries(self):
        """Setup some test raw binary files"""
        
        datadir = NamedTemporaryDirectory()
        
        # Location lookup rawbinaries

        # --- make fixed location lookup
        fixed_location_lookup_filename = os.path.join( self.pathname, "fixed_location_lookup_filename.bin")
        testuuid = uuid.UUID('00000000-0000-0000-0000-000000244c00')
        testdata = numpy.array( [ [   52.0,  -9.0, 22.0 ],
                                  [ -179.0, 133.0,  8.2 ] ], numpy.float64)
        LocationLookupRawBinaryWriter().write(fixed_location_lookup_filename, LocationLookupWithID(testuuid, testdata))

        # --- make mobile location lookups
        mobile_location_lookup_filename_1 = os.path.join( self.pathname, "mobile_location_lookup_filename.20170601.bin")
        testuuid = uuid.UUID('00000000-0000-0000-0000-000000244c01')
        testdata = numpy.array( [ [   52.0,  -9.0, 22.0 ],
                                  [ -179.0, 133.0,  8.2 ] ], numpy.float64)
        LocationLookupRawBinaryWriter().write(mobile_location_lookup_filename_1, LocationLookupWithID(testuuid, testdata))

        mobile_location_lookup_filename_2 = os.path.join( self.pathname, "mobile_location_lookup_filename.20170602.bin")
        testuuid = uuid.UUID('00000000-0000-0000-0000-000000244c02')
        testdata = numpy.array( [ [   52.0,  -9.0, 22.0 ],
                                  [ -179.0, 133.0,  8.2 ] ], numpy.float64)
        LocationLookupRawBinaryWriter().write(mobile_location_lookup_filename_2, LocationLookupWithID(testuuid, testdata))
        
        # Correlation range raw binaries
        
        # --- fixed source
        source_1_correlation_range_filename = os.path.join( self.pathname, "fixed_location_correlation_ranges.bin")
        testdata = [ 0.99, 0.88 ]
        LocalCorrelationRangeRawBinaryWriter().write(source_1_correlation_range_filename, testdata)

        # --- mobile source
        source_2_correlation_range_filename = os.path.join( self.pathname, "mobile_location_correlation_ranges.bin")
        testdata = [ 1.99, 0.00 ]
        LocalCorrelationRangeRawBinaryWriter().write(source_2_correlation_range_filename, testdata)

        # Observation files
        
        writer =  ObservationRawBinaryWriter()
        
        # Store two days for two data sources - 4 files
        
        # we'll just use these multiple times
        location = numpy.array([ 0, 1, 2, 2 ], numpy.uint64)
        measurements = numpy.array([ 275.0, 288.0, 204.0, 310.0 ], numpy.float32)
        uncorrelatederror = numpy.array([5.0, 5.0, 7.0, 8.0], numpy.float32)
        locallycorrelatederror = numpy.array([ [ 2.3, 3.3,  1.2, 8.9 ],
                                               [ 0.222, 7.6, 8.888, 9.999 ] ], numpy.float32)
        mask = numpy.array([ False, False, False, False ], numpy.bool)
        
        # -- source 1
        time_index = int( epoch.days_since_epoch( datetime.datetime(2017, 6, 1) ) )
        time = numpy.array([ time_index, time_index, time_index, time_index ], numpy.float32)
        obs = Observations(mask, time, location, measurements, uncorrelatederror, locallycorrelatederror)
        writer.write_day(os.path.join( self.pathname, "source1.20170601.bin"), uuid.UUID('00000000-0000-0000-0000-000000244c00'), obs, time_index)
        
        time_index = int( epoch.days_since_epoch( datetime.datetime(2017, 6, 2) ) )
        time = numpy.array([ time_index, time_index, time_index, time_index ], numpy.float32)
        obs = Observations(mask, time, location, measurements, uncorrelatederror, locallycorrelatederror)
        writer.write_day(os.path.join( self.pathname, "source1.20170602.bin"), uuid.UUID('00000000-0000-0000-0000-000000244c00'), obs, time_index)
        
        # -- source 2
        
        
        time_index = int( epoch.days_since_epoch( datetime.datetime(2017, 6, 1) ) )
        time = numpy.array([ time_index, time_index, time_index, time_index ], numpy.float32)
        obs = Observations(mask, time, location, measurements, uncorrelatederror, locallycorrelatederror)
        writer.write_day(os.path.join( self.pathname, "source2.20170601.bin"), uuid.UUID('00000000-0000-0000-0000-000000244c01'), obs, time_index)
        
        time_index = int( epoch.days_since_epoch( datetime.datetime(2017, 6, 2) ) )
        time = numpy.array([ time_index, time_index, time_index, time_index ], numpy.float32)
        obs = Observations(mask, time, location, measurements, uncorrelatederror, locallycorrelatederror)
        writer.write_day(os.path.join( self.pathname, "source2.20170602.bin"), uuid.UUID('00000000-0000-0000-0000-000000244c02'), obs, time_index)
        

    def setup_input_descriptors(self):
        
        import datetime
        from collections import OrderedDict
        import dateutil.parser
        
        input_list = [ ["20160612", []],
                       ["20170601", [
                                    {  "local_correlation_ranges_filenames":  {"Tmean": os.path.join(self.pathname, "fixed_location_correlation_ranges.bin")},
                                        "observable_filenames":               {"Tmean": os.path.join(self.pathname, "source1.20170601.bin")},
                                        "mobile_location_lookup_filenames":   {"Tmean": os.path.join(self.pathname, "fixed_location_lookup_filename.bin")},
                                    },
                                    {  "local_correlation_ranges_filenames":  {"Tmean": os.path.join(self.pathname, "mobile_location_correlation_ranges.bin")},
                                        "observable_filenames":               {"Tmean": os.path.join(self.pathname, "source2.20170601.bin")},
                                        "fixed_location_lookup_filename":     os.path.join(self.pathname, "mobile_location_lookup_filename.20170601.bin")
                                    }
                                    ]],
                       ["20170602", [
                                    {  "local_correlation_ranges_filenames":  {"Tmean": os.path.join(self.pathname, "fixed_location_correlation_ranges.bin")},
                                        "observable_filenames":               {"Tmean": os.path.join(self.pathname, "source1.20170602.bin")},
                                        "mobile_location_lookup_filenames":   {"Tmean": os.path.join(self.pathname, "fixed_location_lookup_filename.bin")},
                                    },
                                    {  "local_correlation_ranges_filenames":  {"Tmean": os.path.join(self.pathname, "mobile_location_correlation_ranges.bin")},
                                        "observable_filenames":               {"Tmean": os.path.join(self.pathname, "source2.20170602.bin")},
                                        "fixed_location_lookup_filename":     os.path.join(self.pathname, "mobile_location_lookup_filename.20170602.bin")
                                    }
                                    ]],
                       ["20170603", []],
                       ["20170604", []], ]

        self.input_descriptors = OrderedDict( input_list )
