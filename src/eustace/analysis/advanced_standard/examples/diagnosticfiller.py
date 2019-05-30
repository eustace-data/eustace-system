"""Utilities for automatic generation of operation descriptors in Fullstace"""
import argparse
import collections
import copy
import json
import numpy
import os.path
from eustace.outputformats import definitions

from jsonfiller import FullstaceJsonFiller

class DiagnosticJsonFiller(FullstaceJsonFiller):
    """Fill up operations json descriptors starting from a given number of templates"""
       
    STORAGE_LIST = ['storage_climatology', 'storage_large_scale', 'storage_local']
    UNCERTAINTIES_LIST = ['uncertainties_climatology', 'uncertainties_large_scale', 'uncertainties_local']
    SAMPLES_LIST = ['sample_climatology', 'sample_large_scale', 'sample_local']
    SAMPLES_LIST_PRIOR = ['sample_climatology_prior', 'sample_large_scale_prior', 'sample_local_prior']
    
    INPUTSOURCE_LIST = ['insitu_land', 'insitu_ocean', 'surfaceairmodel_land', 'surfaceairmodel_ocean', 'surfaceairmodel_ice']
    #INPUTSOURCE_LIST = ['insitu_land', 'insitu_ocean']
    
    STATE_SAMPLE_SIZE  = 30
    OUTPUT_SAMPLE_SIZE = definitions.GLOBAL_SAMPLE_SHAPE[3]
    
    def __init__(self, n_iterations, 
                       inpath, outpath,
                       operation_template, batch_operation_template, output_grid_template,
                       start_date, end_date,
                       land_biases, global_biases, 
                       measurement_climatology_command, measurement_climatology_name, solution_climatology_command, solution_climatology_name,
                       measurement_large_scale_command, measurement_large_scale_name, solution_large_scale_command, solution_large_scale_name,
                       local_scale_command, solution_local_scale_name,
                       output_grid_command, output_grid_name,
                       output_climatology_name, output_large_scale_name, output_local_name):
   
        self.n_iterations = n_iterations
        self.inpath = inpath
        self.outpath = outpath
        self.operation_template = operation_template
        self.batch_operation_template = batch_operation_template
        self.output_grid_template = output_grid_template
        self.start_date = start_date+'000000'
        self.end_date = end_date+'000000'
 
        self.land_biases = land_biases
        self.global_biases = global_biases
        self.measurement_climatology_command = measurement_climatology_command
        self.measurement_climatology_name = measurement_climatology_name
        self.solution_climatology_command = solution_climatology_command
        self.solution_climatology_name = solution_climatology_name
        self.measurement_large_scale_command = measurement_large_scale_command
        self.measurement_large_scale_name = measurement_large_scale_name
        self.solution_large_scale_command = solution_large_scale_command
        self.solution_large_scale_name = solution_large_scale_name
        self.local_scale_command = local_scale_command
        self.solution_local_scale_name = solution_local_scale_name
               
        self.output_grid_command = output_grid_command
        self.output_grid_name = output_grid_name
        
        self.output_climatology_name = output_climatology_name
        self.output_large_scale_name = output_large_scale_name
        self.output_local_name = output_local_name
    
    def fill_analysissystem_save_operation(self, iteration_index, component_index, operation_command, dataset_name):
        
        new_dict = copy.deepcopy(self.batch_operation_dict)
        command = operation_command+'_analysissytem_'+str(iteration_index)
        name = dataset_name+'_'+str(iteration_index)
       
        #self.fill_primary_keys(new_dict, command, name, 'pickle', 'spacetime_solution')
        
        new_dict['runmodule']['component_index'] = component_index
        #new_dict['step']={ 'python_class' : 'eumopps.catalogue.step.StepOnce'}
        
        new_dict['step']={ 'python_class' : 'eustace.analysis.advanced_standard.examples.StepAnnual'}
        
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.station_bias_model.store_analysissystem'

        new_dict['name'] = command
        new_dict['runmodule']['insitu_biases'] = self.land_biases
        new_dict['runmodule']['global_biases'] = self.global_biases
        new_dict['runmodule']['compute_uncertainties'] = 0
        new_dict['runmodule']['method'] = 'EXACT'
        new_dict['runmodule']['compute_sample'] = 0
        new_dict['runmodule']['sample_size'] = DiagnosticJsonFiller.STATE_SAMPLE_SIZE
        new_dict['runmodule']['compute_prior_sample'] = 0
        
        for point, value in zip(['start', 'end'], [self.start_date, self.end_date]):
            new_dict['step'][point] = value
        
        datasets = new_dict['newdatasets'][0]
        datasets['name'] = name
        subsets = datasets['subsets'][0]
        subsets['layout']['patterns'] = [ name, name+'_%Y.pickle' ]
        
        # add instructions for state file reading
        self.fill_state_inputs(new_dict, iteration_index)
        
        for to_remove in ['inputsources', 'time_index', 'time_indices']:                                                                                                   
            if to_remove in new_dict['runmodule']:
                del new_dict['runmodule'][to_remove]
                
        return new_dict

    def fill_state_inputs(self, new_dict, iteration, component_filter_list = None):
        """Fill secondary, less common keys
        
        Component reading for other components
        
        TODO check the storage classes for each component
        
        """
        
        # loop through components filling in EUMOPPS instructions to retrieve state, uncertainty and sample filenames
        for storage, solution_datasetname, uncertainty_datasetname, sample_datasetname, sample_prior_datasetname in zip(FullstaceJsonFiller.STORAGE_LIST, 
                                                                                              [self.solution_climatology_name, self.solution_large_scale_name, self.solution_local_scale_name], 
                                                                                              FullstaceJsonFiller.UNCERTAINTIES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST_PRIOR):

            if component_filter_list is not None:
                if solution_datasetname not in component_filter_list:
                    continue

            print new_dict['runmodule'][storage].keys()
            
            if solution_datasetname != self.solution_local_scale_name:
                # space time components
                new_dict['runmodule'][storage]['statefilename_read'] = {'python_class' : 'eumopps.catalogue.placeholder.InputFile',
                                                                        'datasetname': solution_datasetname+'_'+str(iteration)}
                if iteration == self.n_iterations - 1:
                    # add samples for last iteration
                    new_dict['runmodule'][storage]['sample_filename_read'] = {'python_class' : "eumopps.catalogue.placeholder.InputFile",
                                                                            'datasetname': solution_datasetname+'_sample_'+str(iteration)}
                    new_dict['runmodule'][storage]['prior_sample_filename_read'] = {'python_class': 'eumopps.catalogue.placeholder.InputFile',
                                                                                    'datasetname':   solution_datasetname+'_prior_sample_'+str(iteration)}
            else:
                # daily local compontent
                new_dict['runmodule'][storage]['statefilelist_read'] =         {"python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays",
                                                                                                 "datasetname" : solution_datasetname+'_'+str(iteration),
                                                                                                 "missing_data":"allowed"
                                                                                                }
                new_dict['runmodule'][storage]['time_list'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }
                
                if iteration == self.n_iterations - 1:
                    # add samples for last iteration
                    new_dict['runmodule'][storage]['sample_filelist_read'] =       {"python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays",
                                                                                                    "datasetname" : solution_datasetname+'_sample_'+str(iteration),
                                                                                                    "missing_data":"allowed"
                                                                                                    }
                    new_dict['runmodule'][storage]['prior_sample_filelist_read'] = {"python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays",
                                                                                                    "datasetname" : solution_datasetname+'_prior_sample_'+str(iteration),
                                                                                                    "missing_data":"allowed"
                                                                                                    }
                

    def fill_extract_station_biases_operation(self, iteration_index, component_index, operation_command, dataset_name):
        
        new_dict = copy.deepcopy(self.operation_dict)
        command = operation_command+'_'+str(iteration_index)
        name = dataset_name+'_'+str(iteration_index)
       
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'spacetime_solution')
        
        new_dict['runmodule']['component_index'] = component_index
        new_dict['step']={ 'python_class' : 'eumopps.catalogue.step.StepOnce'}
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.station_bias_model.extract_station_biases'

        new_dict['name'] = command
        new_dict['runmodule']['insitu_biases'] = self.land_biases
        new_dict['runmodule']['global_biases'] = self.global_biases
        new_dict['runmodule']['compute_uncertainties'] = 0
        new_dict['runmodule']['method'] = 'EXACT'
        new_dict['runmodule']['compute_sample'] = 0
        new_dict['runmodule']['sample_size'] = DiagnosticJsonFiller.STATE_SAMPLE_SIZE
        new_dict['runmodule']['compute_prior_sample'] = 0
        
        new_dict['runmodule']['first_year'] = self.start_date[:4]
        new_dict['runmodule']['last_year'] = self.end_date[:4]
        new_dict['runmodule']['outputfile'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : dataset_name+'_'+str(iteration_index) }
        
        #for point, value in zip(['start', 'end'], [self.start_date, self.end_date]):
            #new_dict['step'][point] = value
        
        datasets = new_dict['newdatasets'][0]
        datasets['name'] = name
        subsets = datasets['subsets'][0]
        subsets['layout']['patterns'] = [ name, name+'.pickle' ]
        
        self.fill_state_inputs(new_dict, iteration_index, component_filter_list = [self.solution_large_scale_name])
        
        for to_remove in ['inputsources', 'time_index', 'time_indices']:                                                                                                   
            if to_remove in new_dict['runmodule']:
                del new_dict['runmodule'][to_remove]
                
        return new_dict

    """
    
    Early look gridding
    
    """

    def fill_early_look_batch_output_gridding(self, iteration_index):
        """Fill early look output gridding operation without uncertianty information"""
       
        new_dict = copy.deepcopy(self.batch_operation_dict)
       
        self.fill_primary_keys(new_dict, self.output_grid_command+'_preview_'+str(iteration_index), self.output_grid_name+'_preview_'+str(iteration_index), 'nc', 'measurement')
        new_dict['runmodule']['sample_size'] = FullstaceJsonFiller.OUTPUT_SAMPLE_SIZE
        
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.early_look_grid_batch'
        new_dict['runmodule']['outputfiles'] = { 'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput', 'datasetname' : self.output_grid_name+'_preview_'+str(iteration_index)  }
        
        for storage, solution_datasetname, uncertainty_datasetname, sample_datasetname, sample_prior_datasetname in zip(FullstaceJsonFiller.STORAGE_LIST, 
                                                                                              [self.solution_climatology_name, self.solution_large_scale_name, self.solution_local_scale_name], 
                                                                                              FullstaceJsonFiller.UNCERTAINTIES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST_PRIOR):
                                                                                                  
            print new_dict['runmodule'][storage].keys()
            
            if solution_datasetname != self.solution_local_scale_name:
                new_dict['runmodule'][storage]['statefilename_read'] = {'python_class' : 'eumopps.catalogue.placeholder.InputFile',
                                                                        'datasetname': solution_datasetname+'_'+str(iteration_index)}
            else:
                new_dict['runmodule'][storage]['statefilelist_read'] =         {"python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays",
                                                                                                 "datasetname" : solution_datasetname+'_'+str(iteration_index),
                                                                                                 "missing_data":"allowed"
                                                                                                }
                new_dict['runmodule'][storage]['time_list'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }

        
        # remove unnecessary keys
        for to_remove in ['inputsources', 'time_index']:                                                                                                   
            if to_remove in new_dict['runmodule']:
                del new_dict['runmodule'][to_remove]
                 
        new_dict['runmodule']['time_indices'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }
        
        return new_dict

    """
    
    Diagnostic descriptor creation
    
    """

    def fill_diagnostic_descriptor(self):
        self.load_templates()
        
        if not self.n_iterations:
            message = 'Number  of iterations must be a positive integer, at least 1.'
            raise ValueError(message)
        operations = []

        
        for iteration_index in range(self.n_iterations):
            #system_save = self.fill_analysissystem_save_operation(iteration_index, 0, 'store_analysis', 'analysis_storage')
            #operations += [system_save]
            
            bias_save = self.fill_extract_station_biases_operation(iteration_index, 1, 'store_station_biases', 'station_biases')
            operations += [bias_save]
            
            #early_grid = self.fill_early_look_batch_output_gridding(iteration_index)
            #operations += [early_grid]
        
        final_dictionary = {'python_class' : 'eumopps.catalogue.catalogue.Catalogue', 'operations':operations }
        return final_dictionary


    def create_diagnostic_descriptor(self, name, check_none=False):
        """Crate a descriptor json file, name it using user defined name"""
             
        output_filename = os.path.join(self.outpath, name)
        print(output_filename)
        final_dictionary = self.fill_diagnostic_descriptor()
        
        with open(output_filename, 'w') as opt:
            json.dump(final_dictionary, opt, indent=4, separators=(',', ': '), sort_keys=True)
        if check_none:
            print('Checking that all relevant fields have been filled')
            self.check_none(final_dictionary)
        message = 'Creating json descriptor for Fullstace sub-system. Descriptor name is {}'.format(output_filename)
        print(message)
        
