"""Testing model resolutions

For out of core operation Pardiso needs export MKL_PARDISO_OOC_PATH=/work/scratch/$USER/
and the Pardiso wrapper needs to recieve out_of_core=True

"""

import numpy
from datetime import datetime

from eustace.timeutils.decimaltime import datetime_to_decimal_year
from moving_climatology import ModelSetup
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.components.storage_inmemory import SpaceTimeComponentSolutionStorage_InMemory
from eustace.analysis.advanced_standard.elements.combination import CombinationElement, CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.kronecker_annual import AnnualKroneckerElement
from eustace.analysis.advanced_standard.elements.factor import SpaceTimeSPDEHyperparameters

def round_to_quantisation(value, quantisation):
    return (value+quantisation/2) // quantisation * quantisation

def floor_to_quantisation(value, quantisation):
    return value // quantisation * quantisation

def ceil_to_quantisation(value, quantisation):
    return numpy.ceil( value / quantisation ) * quantisation

class SlowSetupT2M1(ModelSetup):
    
    def __init__(self):
        # Note: node spacing is constant in seconds, not months.
        # Months are not the same length in seconds so nodes
        # are not exactly place at month starts.
        
        self.n_triangulation_divisions=2
        self.alpha=2
        
        # Desired minimum model period
        start_time = datetime_to_decimal_year( datetime(1825, 1, 1) )
        end_time   = datetime_to_decimal_year( datetime(2035, 1, 1) )
        
        # Model time resolution in years
        node_spacing = 1./12.
        
        # Model start/end shifted to nearst grid times encompasing the desred start and end times
        self.starttime = floor_to_quantisation(start_time, node_spacing)
        self.endtime = ceil_to_quantisation(end_time, node_spacing)
        
        # number of model nodes in period
        self.n_nodes = int( round( (self.endtime - self.starttime) / node_spacing + 1 ) )
        
        self.overlap_factor=2.5
        self.H = 1
        self.amplitude=2.
        self.space_length_scale=25 #5.0  # length scale in units of degrees
        self.time_length_scale= 10   # length scale in units of years
        
class SlowSetupT5Y5(ModelSetup):
    
    def __init__(self):
        # Note: node spacing is constant in seconds, not months.
        # Months are not the same length in seconds so nodes
        # are not exactly place at month starts.
        
        self.n_triangulation_divisions=5
        self.alpha=2
        
        # Desired minimum model period
        start_time = datetime_to_decimal_year( datetime(1840, 1, 1) )
        end_time   = datetime_to_decimal_year( datetime(2025, 1, 1) )
        
        # Model time resolution in years
        node_spacing = 5.
        
        # Model start/end shifted to nearst grid times encompasing the desred start and end times
        self.starttime = floor_to_quantisation(start_time, node_spacing)
        self.endtime = ceil_to_quantisation(end_time, node_spacing)
        
        # number of model nodes in period
        self.n_nodes = int( round( (self.endtime - self.starttime) / node_spacing + 1 ) )
        
        self.overlap_factor=2.5
        self.H = 1
        self.amplitude=2.
        self.space_length_scale=25 #5.0  # length scale in units of degrees
        self.time_length_scale= 10   # length scale in units of years
        

class TestComponentDefinition(ComponentStorage_InMemory):
    """Define climatology component."""

    def __init__(self):
        setup = SlowSetupT2M1()

        model_elements = CombinationElement( [AnnualKroneckerElement( setup.n_triangulation_divisions, 
                                                                        setup.alpha, 
                                                                        setup.starttime, 
                                                                        setup.endtime, 
                                                                        setup.n_nodes, 
                                                                        setup.overlap_factor, 
                                                                        setup.H ),] )

        model_hyperparameters = CombinationHyperparameters( [SpaceTimeSPDEHyperparameters(numpy.log(setup.amplitude), 
                                                                                            numpy.log(numpy.radians(setup.space_length_scale)), 
                                                                                            numpy.log(setup.time_length_scale)),] )

        super(TestComponentDefinition, self).__init__(model_elements, model_hyperparameters)

def run_prior_solve():

    componentstorage = TestComponentDefinition()    
    solutionstorage = SpaceTimeComponentSolutionStorage_InMemory()
    printstats=True
    compute_uncertainties = False
    method='EXACT'
    compute_sample=False
    sample_size=1
    compute_prior_sample=True
    
    component = SpaceTimeComponent(componentstorage, solutionstorage, printstats, compute_uncertainties, method, compute_sample, sample_size, compute_prior_sample)  
    component_solution = component.component_solution()
    component_solution.update()
    
if __name__ == '__main__':
    
    run_prior_solve()