"""ModelSetup classes for storing model configurations"""

from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.analysissystem import AnalysisSystem

from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.elements.spatial_bias import SpatialBiasElement
from eustace.analysis.advanced_standard.elements.bias_insitu_land import InsituLandBiasElement
from eustace.analysis.advanced_standard.elements.combination import CombinationElement
from eustace.analysis.advanced_standard.elements.combination import CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters
from eustace.analysis.advanced_standard.elements.covariate import LoadCovariateElement
from eustace.analysis.advanced_standard.elements.kronecker import SpaceTimeKroneckerElement
from eustace.analysis.advanced_standard.elements.kronecker_annual import AnnualKroneckerElement
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeSPDEHyperparameters
from eustace.analysis.advanced_standard.elements.grandmean import GrandMeanElement
from eustace.analysis.advanced_standard.elements.latitudeharmonics import LatitudeHarmonicsElement
from eustace.analysis.advanced_standard.elements.latitudespline import LatitudeSplineElement
from eustace.analysis.advanced_standard.elements.local import LocalElement
from eustace.analysis.advanced_standard.elements.local import LocalHyperparameters
from eustace.analysis.advanced_standard.elements.local_view import NonStationaryLocal, ExpandedLocalHyperparameters, ExtendedCombinationHyperparameters
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalElement
from eustace.analysis.advanced_standard.elements.seasonal import SeasonalHyperparameters

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent

from eustace.timeutils.decimaltime import datetime_to_decimal_year
from datetime import datetime
from dateutil import relativedelta
from dateutil.rrule import rrule, MONTHLY
import eustace.timeutils.epoch

import numpy

class ModelSetup(object):
    
    def build_report(self):
        output = []
        
        for name, item in self.__dict__.iteritems():
            if isinstance(item, ModelSetup):
                output.append("".join( (name, ":", str(item)) ))
                output.extend(["".join( (name, ".", value) ) for value in item.build_report()])
            else:
                output.append("".join( (name, ":", str(item)) ))
            
        return output
                
    def report(self):
        print str(self)
        for item in self.build_report():
            print item
        
class AnalysisModelSetup(ModelSetup):
    def __init__(self):
        self.climatology_settings = MovingClimatologySetup()
        self.large_scale_settings = MidScaleSetup()
        self.local_settings = ShortScaleSetup()

"""

ModelSetup classes for individual model components

"""

class MovingClimatologySetup(ModelSetup):
    """Define the setup for the climatology component"""
    
    def __init__(self):
        
        self.seasonal_settings  = self.SeasonalSetup()
        self.seasonal_spline_settings = self.SeasonalSpline()
        self.covariate_settings = self.CovariateSetup()
        self.bias_settings      = self.BiasSetup()        
        #self.slow_settings      = self.SlowSetup()
        self.latitude_settings  = self.LatitudeSpline()
    
    class SeasonalSetup(ModelSetup):
        def __init__(self):
            self.n_triangulation_divisions=6
            self.n_harmonics=3 # 4
            self.include_local_mean = True
            self.n_spatial_components=self.n_harmonics +self.include_local_mean
            
            # placeholders for initial setup - will be overwritten below
            self.amplitude=2.
            self.space_length_scale=10. #length scale in units of degrees - prev 5.0
            
            # actual values to be used ordered [local_mean, 1st harmonic order, 2nd harmonic order, 3rd harmonic order]
            self.harmonic_amplitudes = [1.0, 8.0, 4.0, 2.0]#[2.0, 20.0, 10.0, 5.0]
            self.harmonic_length_scales = [5.0, 15.0, 10.0, 10.0]
    
    class LatitudeSpline(ModelSetup):
        def __init__(self):
            self.amplitude=12.
            self.length_scale=30.#15.
            self.overlap_factor=2.5
            self.H = 1.
            self.alpha = 2
            self.n_nodes = 361#73
            
    class SeasonalSpline(ModelSetup):
        def __init__(self):
            # Does not include a local offset model - not needed for random field representation? Maybe for higher frequency variation.
            
            self.amplitude=5.
            #self.space_length_scale=10.# 5.#10. #length scale in units of degrees - prev 5.0
            self.space_length_scale=25.# 5.#10. #length scale in units of degrees - prev 5.0
            self.n_triangulation_divisions=6
            self.alpha = 2
            
            self.time_length_scale=2. / 12.0 # lengthscale in years
            self.n_nodes = 12
            
            self.overlap_factor=2.5
            self.H = 1.0
            
            self.starttime = 0.0 # year start in decimal years
            self.endtime = 1.0   # year end in decimal years
            self.wrap_dimensions = [True,] # list of dimensions to wrap. only one dimension in time.
            

    class CovariateSetup(ModelSetup):
        def __init__(self):
            self.grandmean_amplitude=2.0 * 273.15

    
    class BiasSetup(ModelSetup):
        def __init__(self):
            self.bias_amplitude=.9

class MidScaleSetup(ModelSetup):
    
    def __init__(self):
        
        self.midscale_settings = self.MonthlySetup()
        self.bias_settings = self.BiasSetup()
        self.slow_settings      = self.SlowSetup()
        
    class MonthlySetup(ModelSetup):
        def __init__(self):
            # Note: node spacing is constant in seconds, not months.
            # Months are not the same length in seconds so nodes
            # are not exactly placed at month starts.
            time_extension = relativedelta.relativedelta(years = 2)
            
            #model_startdate = datetime(1850, 1, 1)-time_extension#datetime(1995, 1, 1)-time_extension # datetime(1848, 1, 1)-time_extension
            model_startdate = datetime(1880, 1, 1)-time_extension
            model_enddate = datetime(2016, 1, 1)+time_extension # datetime(1862, 1, 1)+time_extension
            avg_months_per_node = 3
            
            self.n_triangulation_divisions=4
            self.alpha=2
            self.starttime=eustace.timeutils.epoch.days_since_epoch( model_startdate )
            self.endtime=eustace.timeutils.epoch.days_since_epoch( model_enddate )
            self.n_nodes = (len( list(rrule(MONTHLY, dtstart=model_startdate, until=model_enddate)) ) - 1 ) / avg_months_per_node +1
            self.overlap_factor=2.5
            self.H = 1
            self.amplitude=1.0
            #self.space_length_scale=10.0 #10.0  # length scale in units of degrees
            self.space_length_scale=15.0 #10.0  # length scale in units of degrees
            self.time_length_scale=135.0   # length scale in units of days

    class SlowSetup(ModelSetup):
        def __init__(self):
            # Note: node spacing is constant in seconds, not months.
            # Months are not the same length in seconds so nodes
            # are not exactly place at month starts.
            
            def round_to_quantisation(value, quantisation):
                return (value+quantisation/2) // quantisation * quantisation
            
            def floor_to_quantisation(value, quantisation):
                return value // quantisation * quantisation
            
            def ceil_to_quantisation(value, quantisation):
                return numpy.ceil( value / node_spacing ) * node_spacing
            
            self.n_triangulation_divisions=4
            self.alpha=2
            
            # Desired minimum model period
            #start_time = datetime_to_decimal_year( datetime(1830, 1, 1) )
            start_time = datetime_to_decimal_year( datetime(1860, 1, 1) )
            end_time   = datetime_to_decimal_year( datetime(2035, 1, 1) )
            
            # Model time resolution in years
            node_spacing = 4.0
            
            # Model start/end shifted to nearst grid times encompasing the desred start and end times
            self.starttime = floor_to_quantisation(start_time, node_spacing)
            self.endtime = ceil_to_quantisation(end_time, node_spacing)
            
            # number of model nodes in period
            self.n_nodes = int( round( (self.endtime - self.starttime) / node_spacing + 1 ) )
            
            self.overlap_factor=2.5
            self.H = 1
            self.amplitude=0.5
            self.space_length_scale=30.#25 #5.0  # length scale in units of degrees
            self.time_length_scale= 10.   # length scale in units of years
            
    

    class BiasSetup(ModelSetup):
        def __init__(self):
            self.bias_amplitude=.9

class ShortScaleSetup(ModelSetup):
    
    def __init__(self, hyperparameter_file = None):
        
        if hyperparameter_file is None:
            self.local_settings = self.StationaryLocalSetup()
        else:
            self.local_settings = self.NonStationaryLocalSetup(hyperparameter_file)
            
        self.bias_settings = self.BiasSetup()
        self.spatial_bias_settings = self.SpatialBiasSetup()
        
    class StationaryLocalSetup(ModelSetup):
        def __init__(self):
            self.n_triangulation_divisions=7
            self.amplitude=2.0# 3.0.
            self.space_length_scale=5.0 #2.0   #length scale in units of degrees
    
    class NonStationaryLocalSetup(ModelSetup):
        def __init__(self, hyperparameter_file):
            self.n_triangulation_divisions=7
            self.amplitude=2.
            self.space_length_scale=5.0 #2.0   #length scale in units of degrees
            
            self.hyperparameters = numpy.load(hyperparameter_file)
    
    class BiasSetup(ModelSetup):
        def __init__(self):
            self.bias_amplitude=2.0#15.0
            
    class SpatialBiasSetup(ModelSetup):
        def __init__(self):
            self.n_triangulation_divisions=6
            self.spatial_bias_amplitutde=1.0#15.0
            self.spatial_bias_length_scale=25.0


"""

Component definitions for each model component

"""

class ClimatologyDefinition(ComponentStorage_InMemory):
    """Define climatology component."""

    def __init__(self, covariates_descriptor):
        if covariates_descriptor is not None:
            loader = LoadCovariateElement(covariates_descriptor)
            loader.check_keys()
            covariate_elements, covariate_hyperparameters = loader.load_covariates_and_hyperparameters()
            print('The following fields have been added as covariates of the climatology model')
            print(loader.data.keys())
        else:
            covariate_elements, covariate_hyperparameters = [], []

        setup = MovingClimatologySetup()
    
        model_elements = CombinationElement( [AnnualKroneckerElement(setup.seasonal_spline_settings.n_triangulation_divisions, 
                                                          setup.seasonal_spline_settings.alpha, 
                                                          setup.seasonal_spline_settings.starttime, 
                                                          setup.seasonal_spline_settings.endtime, 
                                                          setup.seasonal_spline_settings.n_nodes, 
                                                          setup.seasonal_spline_settings.overlap_factor, 
                                                          setup.seasonal_spline_settings.H,
                                                          setup.seasonal_spline_settings.wrap_dimensions),
                                                #SeasonalElement(setup.seasonal_settings.n_triangulation_divisions, 
                                                 #setup.seasonal_settings.n_harmonics, 
                                                 #include_local_mean=setup.seasonal_settings.include_local_mean),
                                 GrandMeanElement(),
                                 LatitudeSplineElement(setup.latitude_settings.alpha,
                                                       setup.latitude_settings.n_nodes,
                                                       setup.latitude_settings.overlap_factor,
                                                       setup.latitude_settings.H,),]
                                 + covariate_elements)
        
        
        
        seasonal_hyperparameters = SeasonalHyperparameters(setup.seasonal_settings.n_spatial_components, 
                                                                                     numpy.log(setup.seasonal_settings.amplitude), 
                                                                                     numpy.log(numpy.radians(setup.seasonal_settings.space_length_scale)))
        
        seasonal_spline_hyperparameters = SpaceTimeSPDEHyperparameters(numpy.log(setup.seasonal_spline_settings.amplitude), 
                                                                                          numpy.log(numpy.radians(setup.seasonal_spline_settings.space_length_scale)), 
                                                                                          numpy.log(setup.seasonal_spline_settings.time_length_scale))
        

        seasonal_params = zip( numpy.log(setup.seasonal_settings.harmonic_amplitudes),
                               numpy.log(numpy.radians(setup.seasonal_settings.harmonic_length_scales)))
        seasonal_hyperparameters.set_array([val for pair in seasonal_params for val in pair])
        
        model_hyperparameters = CombinationHyperparameters( [seasonal_spline_hyperparameters,
                                                             #seasonal_hyperparameters, 
                                                             CovariateHyperparameters(numpy.log(setup.covariate_settings.grandmean_amplitude)),
                                                             #SpaceTimeSPDEHyperparameters(numpy.log(setup.slow_settings.amplitude), 
                                                                                          #numpy.log(numpy.radians(setup.slow_settings.space_length_scale)), 
                                                                                          #numpy.log(setup.slow_settings.time_length_scale)),
                                                             LocalHyperparameters(numpy.log(setup.latitude_settings.amplitude), 
                                                                                  numpy.log(setup.latitude_settings.length_scale))
                                                                                          ]
                                                             + covariate_hyperparameters)
    
        super(ClimatologyDefinition, self).__init__(model_elements, model_hyperparameters)

class LargeScaleDefinition(ComponentStorage_InMemory):
    """Define large scale component."""

    def __init__(self, bias_terms = False, breakpoints_file = None):
      
        setup = MidScaleSetup()
        
        if bias_terms:
            bias_element = [InsituLandBiasElement(breakpoints_file, apply_policy=True, cut_value=3)]
            bias_hyperparameters = [CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude))]
        else:
            bias_element, bias_hyperparameters = [], []

        super(LargeScaleDefinition, self).__init__(
            CombinationElement([SpaceTimeKroneckerElement(setup.midscale_settings.n_triangulation_divisions, 
                                                        setup.midscale_settings.alpha, 
                                                        setup.midscale_settings.starttime, 
                                                        setup.midscale_settings.endtime, 
                                                        setup.midscale_settings.n_nodes, 
                                                        setup.midscale_settings.overlap_factor, 
                                                        setup.midscale_settings.H),
                                 AnnualKroneckerElement(setup.slow_settings.n_triangulation_divisions, 
                                                          setup.slow_settings.alpha, 
                                                          setup.slow_settings.starttime, 
                                                          setup.slow_settings.endtime, 
                                                          setup.slow_settings.n_nodes, 
                                                          setup.slow_settings.overlap_factor, 
                                                          setup.slow_settings.H),]
                                                          + bias_element),
            CombinationHyperparameters([SpaceTimeSPDEHyperparameters(numpy.log(setup.midscale_settings.amplitude), 
                                                                    numpy.log(numpy.radians(setup.midscale_settings.space_length_scale)), 
                                                                    numpy.log(setup.midscale_settings.time_length_scale)),
                                        SpaceTimeSPDEHyperparameters(numpy.log(setup.slow_settings.amplitude), 
                                                                    numpy.log(numpy.radians(setup.slow_settings.space_length_scale)), 
                                                                    numpy.log(setup.slow_settings.time_length_scale)),]
                                                                    + bias_hyperparameters))
            
class LocalDefinition(ComponentStorage_InMemory):
    """Define local component."""
    
    def __init__(self, bias_terms = False, global_biases_group_list = []):         
    
        setup = ShortScaleSetup()
        
        bias_elements, bias_hyperparameters = [], []
        if bias_terms:
            for groupname in global_biases_group_list:
                #if (groupname == 'surfaceairmodel_land_global') or (groupname == 'surfaceairmodel_ocean_global'):
                    #bias_elements.append( BiasElement(groupname, 1) )
                    #bias_hyperparameters.append( CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude)) )
                if (groupname == 'surfaceairmodel_land_global'):
                    # global mean term
                    bias_elements.append( BiasElement(groupname, 1) )
                    bias_hyperparameters.append( CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude)) )
                    # spatial bias term
                    bias_elements.append( SpatialBiasElement(groupname, setup.spatial_bias_settings.n_triangulation_divisions) )
                    bias_hyperparameters.append( LocalHyperparameters(numpy.log(setup.spatial_bias_settings.spatial_bias_amplitutde), 
                                                            numpy.log(numpy.radians(setup.spatial_bias_settings.spatial_bias_length_scale))) )
                elif (groupname == 'surfaceairmodel_ocean_global'):
                    # global mean term
                    bias_elements.append( BiasElement(groupname, 1) )
                    bias_hyperparameters.append( CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude)) )
                elif (groupname == 'surfaceairmodel_ice_global'):
                    # hemispheric term
                    bias_elements.append( BiasElement(groupname, 2) )
                    bias_hyperparameters.append( CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude)) )
                    # spatial bias term
                    bias_elements.append( SpatialBiasElement(groupname, setup.spatial_bias_settings.n_triangulation_divisions) )
                    bias_hyperparameters.append( LocalHyperparameters(numpy.log(setup.spatial_bias_settings.spatial_bias_amplitutde), 
                                                            numpy.log(numpy.radians(setup.spatial_bias_settings.spatial_bias_length_scale))) )

        super(LocalDefinition, self).__init__(
            CombinationElement([LocalElement(setup.local_settings.n_triangulation_divisions)] + bias_elements),         
            CombinationHyperparameters([LocalHyperparameters(numpy.log(setup.local_settings.amplitude), 
                                                            numpy.log(numpy.radians(setup.local_settings.space_length_scale)))] + bias_hyperparameters))

class NonStationaryLocalDefinition(ComponentStorage_InMemory):
    """Define non stational with stored hyperparameters local component."""
    
    def __init__(self, bias_terms = False, global_biases_group_list = [], local_hyperparameter_file = None):         
    
        setup = ShortScaleSetup(local_hyperparameter_file)
        
        bias_elements, bias_hyperparameters = [], []
        if bias_terms:
            for groupname in global_biases_group_list:
                if (groupname == 'surfaceairmodel_land_global') or (groupname == 'surfaceairmodel_ocean_global'):
                    bias_elements.append( BiasElement(groupname, 1) )
                    bias_hyperparameters.append( CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude)) )
                elif (groupname == 'surfaceairmodel_ice_global'):
                    bias_elements.append( BiasElement(groupname, 2) )
                    bias_hyperparameters.append( CovariateHyperparameters(numpy.log(setup.bias_settings.bias_amplitude)) )

        local_hyperparameters = ExpandedLocalHyperparameters(log_sigma = None, log_rho = None)
        local_hyperparameters.values_from_npy_savefile(local_hyperparameter_file)
        
        super(NonStationaryLocalDefinition, self).__init__(
            CombinationElement([NonStationaryLocal(setup.local_settings.n_triangulation_divisions)] + bias_elements),         
            CombinationHyperparameters([local_hyperparameters]+bias_hyperparameters))



