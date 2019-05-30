from eustace.analysis.advanced_standard.examples.moving_climatology import AnalysisSystem_EUSTACE
from eustace.timeutils.epoch import days_since_epoch
import pickle
import numpy
import datetime
from eustace.outputformats.ensuredirectory import ensuredirectory
from eustace.analysis.advanced_standard.elements.bias import BiasElement
from eustace.analysis.advanced_standard.fileio.observation_source_surface_effects import insitu_land_covariate_effect

BIAS_COMPONENT_INDEX=1

def save_model(outdict, modelfile):
    # load a bias model pickle file
    ensuredirectory(modelfile)
    with open(modelfile, 'w') as f:
        pickle.dump(outdict, f)


def load_model(modelfile):
    # load a bias model pickle file
    
    return pickle.load(open(modelfile, 'r'))

#def extract_station_biases(storage_climatology, storage_large_scale, storage_local, 
            #outputfile, climatologyfile, largescalefile, localfile,
            #processdate, time_index, 
            #covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
            #compute_uncertainties, method,
            #compute_sample, sample_size, compute_prior_sample):
                
def extract_station_biases(storage_climatology, storage_large_scale, storage_local, 
            outputfile, first_year, last_year,
            component_index,
            covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
            compute_uncertainties, method,
            compute_sample, sample_size, compute_prior_sample):
                
    from eustace.preprocess.fileio.insitu_land_breakpoints import ObservationBreakPointSourceInsituLand
    from eustace.analysis.advanced_standard.fileio.output_projector import Projector

    #print 'VERSION: {0}'.format(get_revision_id_for_module(eustace))

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list,
                        compute_uncertainties, method)


    component = analysissystem.components[BIAS_COMPONENT_INDEX]
    element = component.storage.element_read()
    hyperparameters = component.storage.hyperparameters_read()
    
    currentstate = component.solutionstorage.state_read()
    element_states = component.storage.element_read().element_prior(hyperparameters).element_states(currentstate)
    
    breakpoints_reader = ObservationBreakPointSourceInsituLand(breakpoints_file)
    station_locations = breakpoints_reader.observation_location_lookup()
    breakpoints_reader.dataset.close()
    
    #evaluation_dates = [datetime.datetime(year, 1, 1) for year in range(1850, 2016)]
    import pickle
    for i, subelement in enumerate(element.combination):

        if isinstance(subelement, BiasElement):
            mystate = element_states[i]
            
            if hasattr(subelement, 'observed_breakpoints'):
                
                
                # station_indices for each homogenisation segment period 
                expanded_station_indices = subelement.observed_breakpoints.break_station-1 # with station indices shift to start at zero
                
                print station_locations
                print expanded_station_indices
                # station spatial locations expanded for each break point
                expanded_locations = station_locations[:,expanded_station_indices]
                    
                outdict = {'breakpoints':        subelement.observed_breakpoints,                 # the breakpoint object
                           'station_locations':  station_locations.T,                             # one location per station index
                           'expanded_locations': expanded_locations.T,                            # one location per bias estimate
                           'expanded_station_indices': expanded_station_indices,                  # the station index for each bias state
                           'biases':             mystate.ravel(),   # the bias estimates
                           }                                         
               
            
            build_bias_matrix(outdict, int(first_year), int(last_year))
            
            save_model(outdict, outputfile)



def build_bias_matrix(model, first_year, last_year):
    """Evaluate the bias model at the first day of each year for each station and store"""
    
    #model = load_fitted_breakmodel(modelfile)
    
    evaluation_dates = [datetime.datetime(year, 1, 1) for year in range(first_year, last_year+1)]
    
    n_stations = max(model['expanded_station_indices']) + 1
    station_indices = range(n_stations)
    
    time_indices = [days_since_epoch(date) for date in evaluation_dates]
                    
    bias_index = numpy.zeros((n_stations, len(time_indices))) + numpy.nan
    bias_grid = numpy.zeros((n_stations, len(time_indices))) + numpy.nan
    
    for ind, time_index in enumerate(time_indices):
        effect = insitu_land_covariate_effect(time_index, station_indices, model['breakpoints'])
        
        if effect is not None:
            bias_index[effect[:,0],ind] = effect[:,1]                    
            bias_grid[effect[:,0],ind] = model['biases'][effect[:,1]]
    
    model['time_indices'] = time_indices
    model['evaluation_dates']= evaluation_dates    # dates represented by time_indices
    model['bias_index']= bias_index                # inidcies to the biases mapping into the adjustment matrix
    model['bias_grid']=   bias_grid                # the bias matrix


def store_analysissystem(storage_climatology, storage_large_scale, storage_local, 
      component_index, covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
      compute_uncertainties, method,
          compute_sample, sample_size, compute_prior_sample):

    # Build analysis system
    analysissystem = AnalysisSystem_EUSTACE(storage_climatology, storage_large_scale, storage_local, 
                        covariates_descriptor, insitu_biases, breakpoints_file, global_biases, global_biases_group_list, 
                        compute_uncertainties, method,
                        compute_sample, sample_size, compute_prior_sample)

    # Store
    save_model(analysissystem, '/work/scratch/cmorice/advanced_standard/analysis_system.pickle')

def break_series( observed_breakpoints, station_index ):
    # return a time series of biases for a give station
    pass