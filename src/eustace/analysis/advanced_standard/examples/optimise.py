import moving_climatology
from eustace.analysis.advanced_standard.optimizationsystem import OptimizationSystem
from eustace.analysis.advanced_standard.elements.local_view import LocalSubRegion, LocalSuperTriangle, ExtendedCombinationHyperparameters, SphereMeshViewGlobal
from eustace.analysis.advanced_standard.elements.local_view import NonStationaryLocal, ExpandedLocalHyperparameters
from eustace.analysis.advanced_standard.components.storage_inmemory import ComponentStorage_InMemory
from eustace.analysis.advanced_standard.elements.combination import CombinationElement, CombinationHyperparameters
from eustace.analysis.advanced_standard.elements.local import LocalElement, LocalHyperparameters
from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.elements.covariate import CovariateHyperparameters

from eustace.analysis.advanced_standard.components.spatial import SpatialComponent
from eustace.analysis.advanced_standard.components.spacetime import SpaceTimeComponent

from eustace.analysis.observationsource import ObservationSource

from example_optimization import RegionOptimizationSystem_EUSTACE
from example_optimization import optimize_region, merge_local_parameterisations

import numpy
import argparse
import os

from eumopps.jsonobjects import jsonobjects
from eustace.outputformats.ensuredirectory import ensuredirectory

def read_inputdescriptor(filename):
    
    import json
    from collections import OrderedDict
    
    with open(filename) as f:
        inputdescriptor = OrderedDict( json.load(f) )
    
    return inputdescriptor

def local_optimisation(neighbourhood_level = 0, region_index = 0, regionspec = 'LocalSubRegion'):
    """Run the optimisation outside of EUMOPPS with manually setup analysis system"""
    
    from eustace.analysis.advanced_standard.components.storage_files import SpaceTimeComponentSolutionStorage_Files, DelayedSpatialComponentSolutionStorageFlexible_Files
    from eustace.analysis.advanced_standard.components.storage_inmemory import SpatialComponentSolutionStorage_InMemory
    
    # climatology storage
    climatology_dir = "/work/scratch/cmorice/advanced_standard/solution_climatology/"
    climatology_state_file = os.path.join(climatology_dir, "solution_climatology.pickle")
    #climatology_state_file = None
    storage_climatology = SpaceTimeComponentSolutionStorage_Files( statefilename_read = climatology_state_file )
    
    # large scale storage
    large_scale_dir = "/work/scratch/cmorice/advanced_standard/solution_large_scale/"
    large_scale_state_file = os.path.join(large_scale_dir, "solution_large_scale.pickle")
    #large_scale_state_file = None
    storage_large_scale = SpaceTimeComponentSolutionStorage_Files( statefilename_read = large_scale_state_file )
    
    # local bias storage
    local_bias_read_pattern = ["/work/scratch/cmorice/advanced_standard/optimisation_local_bias/%Y/", "solution_bias_%Y%m%d.pickle"]
    storage_local_bias  = DelayedSpatialComponentSolutionStorageFlexible_Files( state_read_pattern = local_bias_read_pattern )
    
    # local spde storage
    storage_region_spde = SpatialComponentSolutionStorage_InMemory()
    
    # climatology covariates
    covariates_descriptor = "/gws/nopw/j04/eustace/data/internal/climatology_covariates/covariates.json"
    
    # bias terms
    insitu_biases         = 1
    breakpoints_file      = "/gws/nopw/j04/eustace/data/internal/D1.7/daily/eustace_stations_global_R001127_daily_status.nc"
    global_biases         = 1
    global_biases_group_list = ["surfaceairmodel_ice_global", "surfaceairmodel_land_global", "surfaceairmodel_ocean_global"]
    
    # uncertainty calculation
    compute_uncertainties    = False
    method = 'APPROXIMATED'
    
    # observation descriptor json file
    input_dir = "/work/scratch/cmorice/advanced_standard/local_hyperparameters/"
    inputdescriptor_file = os.path.join(input_dir, "merged_input_summary.json")
    inputdescriptor = read_inputdescriptor(inputdescriptor_file)
    #print inputdescriptor
    
    # list of time steps used in optimisation
    time_keys = ["2009-12-25T00:00:00","2009-12-26T00:00:00", "2009-12-27T00:00:00",]# "2009-12-28T00:00:00", "2009-12-29T00:00:00", "2009-12-30T00:00:00", "2009-12-31T00:00:00"]
    
    
    time_keys = ["2003-01-01T00:00:00", "2003-02-01T00:00:00", "2003-03-01T00:00:00", "2003-04-01T00:00:00","2003-05-01T00:00:00", "2003-06-01T00:00:00", 
                 "2003-07-01T00:00:00","2003-08-01T00:00:00", "2003-09-01T00:00:00", "2003-10-01T00:00:00","2003-11-01T00:00:00", "2003-12-01T00:00:00",
                 "2004-01-01T00:00:00", "2004-02-01T00:00:00", "2004-03-01T00:00:00", "2004-04-01T00:00:00","2004-05-01T00:00:00", "2004-06-01T00:00:00", 
                 "2004-07-01T00:00:00","2004-08-01T00:00:00", "2004-09-01T00:00:00", "2004-10-01T00:00:00","2004-11-01T00:00:00", "2004-12-01T00:00:00",]
    
    # output hyperparameter file
    output_dir = "/work/scratch/cmorice/advanced_standard/local_optimisation/"
    ensuredirectory(output_dir)
    local_setup = moving_climatology.ShortScaleSetup()
    hyperparameter_storage_file = os.path.join(output_dir, "regional_hyperparameters.%i.%i.%i.npy" % ( local_setup.local_settings.n_triangulation_divisions, neighbourhood_level, region_index ) )



    compute_sample = False
    sample_size = 0

    # run the optimisation for this region - model configuration is a combination of storage and ancilliary files specified here and that imported by optimize_region
    optimize_region(storage_climatology, storage_large_scale, storage_local_bias, storage_region_spde,
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method,
                        compute_sample, sample_size,
                        neighbourhood_level, region_index, regionspec,
                        inputdescriptor, time_keys,
                        hyperparameter_storage_file)

def merge_hyperparameters(neighbourhood_level = 0, regionspec = 'LocalSubRegion'):
    
    from submit_optimisation import n_operations
    
    local_setup = moving_climatology.ShortScaleSetup()
    full_resolution_level = local_setup.local_settings.n_triangulation_divisions
    neighbourhood_level = 3
    regionspec = 'LocalSubRegion'
    
    number_of_regions = n_operations(neighbourhood_level)
    #number_of_regions = 20
    
    # output hyperparameter file
    output_dir = "/work/scratch/cmorice/advanced_standard/local_optimisation/"
    
    output_filename = os.path.join(output_dir, "regional_hyperparameters.%i.%i.npy" % ( local_setup.local_settings.n_triangulation_divisions, neighbourhood_level ) )
    
    hyperparameter_filenames = []
    for region_index in range(number_of_regions):
        hyperparameter_file = os.path.join(output_dir, "regional_hyperparameters.%i.%i.%i.npy" % ( local_setup.local_settings.n_triangulation_divisions, neighbourhood_level, region_index ) )
        hyperparameter_filenames.append( hyperparameter_file )
    
    merge_local_parameterisations(full_resolution_level, neighbourhood_level, hyperparameter_filenames, output_filename)



def parse_and_run(argvalues=None):

    parser = argparse.ArgumentParser()
    parser.add_argument('--neighbourhood_level', required=False, default=3, type=int, help='triangulation subdivision level at which the local region is defined')
    parser.add_argument('--region_index', required=False, default=0, type=int, help='index of region centre triangle at the subdivision level specified in neighbourhood_level')
    parser.add_argument('--regionspec', required=False, default='LocalSubRegion', choices = ['LocalSubRegion', 'LocalSuperTriangle'], help='type of subdivision of "LocalSubRegion" or "LocalSuperTriangle"')
    parser.add_argument('--merge', required=False, default=False, help='run merge of local hyperparameter files')
    args = parser.parse_args(argvalues)

    # call build method using command line args
    #build_from_json_file(**vars(args))
    
    #print vars(args)
    #local_optimisation(**vars(args))
    if args.merge:
        merge_hyperparameters(neighbourhood_level=args.neighbourhood_level, regionspec=args.regionspec)
    else:
        local_optimisation(neighbourhood_level=args.neighbourhood_level, region_index=args.region_index, regionspec=args.regionspec)
        

    
if __name__ == '__main__':
    
    parse_and_run(argvalues=None)
    