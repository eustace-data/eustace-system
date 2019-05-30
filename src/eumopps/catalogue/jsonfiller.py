"""Utilities for automatic generation of operation descriptors in Fullstace"""
import argparse
import collections
import copy
import json
import numpy
import os.path

from eustace.outputformats import definitions

class FullstaceJsonFiller(object):
    """Fill up operations json descriptors starting from a given number of templates"""
       
    STORAGE_LIST = ['storage_climatology', 'storage_large_scale', 'storage_local']
    UNCERTAINTIES_LIST = ['uncertainties_climatology', 'uncertainties_large_scale', 'uncertainties_local']
    SAMPLES_LIST = ['sample_climatology', 'sample_large_scale', 'sample_local']
    SAMPLES_LIST_PRIOR = ['sample_climatology_prior', 'sample_large_scale_prior', 'sample_local_prior']
    
    # Input sources are filtered to those included in this list
    INPUTSOURCE_LIST = ['insitu_land', 'insitu_ocean', 'surfaceairmodel_land', 'surfaceairmodel_ocean', 'surfaceairmodel_ice']
    
    STATE_SAMPLE_SIZE  = 10
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
               
    def load_templates(self):
        """Load json templates from given filepaths"""
       
        with open(os.path.join(self.inpath, self.operation_template)) as opt:
            print('XXXXXXXXXXXXXXXXXXXX')
            print(opt)
            print('XXXXXXXXXXXXXXXXXXXX')
            self.operation_dict = json.load(opt)
            self.restrict_inputs(self.operation_dict)
            
        with open(os.path.join(self.inpath, self.batch_operation_template)) as opt:
            print('YYYYYYYYYYYYYYYYYYYYYY')
            print(opt)
            print('YYYYYYYYYYYYYYYYYYYYYY')

            self.batch_operation_dict = json.load(opt)
            self.restrict_inputs(self.batch_operation_dict)
       
        with open(os.path.join(self.inpath, self.output_grid_template)) as opt:
            self.output_gridding_dict = json.load(opt)
            self.restrict_inputs(self.output_gridding_dict)
    
    def restrict_inputs(self, operation_dict):
        """Restrict the list of inputsources to only those specified in FullstaceJsonFiller.INPUTSOURCE_LIST"""
        
        if 'inputsources' in operation_dict['runmodule']:
            restricted_sources = [inputsource for inputsource in operation_dict['runmodule']['inputsources'] if inputsource['name'] in FullstaceJsonFiller.INPUTSOURCE_LIST]
            
            operation_dict['runmodule']['inputsources'] = restricted_sources
    
    def create_descriptor(self, name, check_none=False):
        """Crate a descriptor json file, name it using user defined name"""
             
        output_filename = os.path.join(self.outpath, name)
        print(output_filename)
        final_dictionary = self.fill_descriptor()
        
        with open(output_filename, 'w') as opt:
            json.dump(final_dictionary, opt, indent=4, separators=(',', ': '), sort_keys=True)
        if check_none:
            print('Checking that all relevant fields have been filled')
            self.check_none(final_dictionary)
        message = 'Creating json descriptor for Fullstace sub-system. Descriptor name is {}'.format(output_filename)
        print(message)
        
    def create_batch_descriptor(self, name, check_none=False):
        """Create a descriptor json file, name it using user defined name"""
             
        output_filename = os.path.join(self.outpath, name)
        print(output_filename)
        final_dictionary = self.fill_batch_descriptor()
        
        with open(output_filename, 'w') as opt:
            json.dump(final_dictionary, opt, indent=4, separators=(',', ': '), sort_keys=True)
        if check_none:
            print('Checking that all relevant fields have been filled')
            self.check_none(final_dictionary)
        message = 'Creating json descriptor for Fullstace sub-system. Descriptor name is {}'.format(output_filename)
        print(message)
        

    def fill_descriptor(self):
        """Fill the descriptor"""
        
        self.load_templates()
        if not self.n_iterations:
            message = 'Number  of iterations must be a positive integer, at least 1.'
            raise ValueError(message)
        operations = []
        for iteration_index in range(self.n_iterations):
            # Filling climatology
            measurement_climatology = self.fill_single_measurement_operation_space_time(iteration_index, 0, self.measurement_climatology_command, self.measurement_climatology_name)
            solution_climatology =self.fill_single_solution_operation_space_time(iteration_index, 0, self.solution_climatology_command, self.solution_climatology_name, self.measurement_climatology_name)
            # Filling large scale
            measurement_large_scale = self.fill_single_measurement_operation_space_time(iteration_index, 1, self.measurement_large_scale_command, self.measurement_large_scale_name)
            solution_large_scale = self.fill_single_solution_operation_space_time(iteration_index, 1, self.solution_large_scale_command, self.solution_large_scale_name, self.measurement_large_scale_name)
            # Filling local scale
            meas_sol_local = self.fill_single_operation_space(iteration_index, 2)
            operations += [measurement_climatology, solution_climatology, measurement_large_scale, solution_large_scale, meas_sol_local]
        self.fill_output_gridding()
        operations += [self.output_gridding_dict]
        
        final_dictionary = {'python_class' : 'eumopps.catalogue.catalogue.Catalogue', 'operations':operations }
        return final_dictionary
    
    
    def fill_batch_descriptor(self):
        """Fill the descriptor"""
        
        self.load_templates()
        if not self.n_iterations:
            message = 'Number  of iterations must be a positive integer, at least 1.'
            raise ValueError(message)
        operations = []
        
        for iteration_index in range(self.n_iterations):
            
            # Filling climatology
            measurement_climatology = self.fill_batch_measurement_operation_space_time(iteration_index, 0, self.measurement_climatology_command, self.measurement_climatology_name)
            solution_climatology =self.fill_single_posterior_solution_operation_space_time(iteration_index, 0, self.solution_climatology_command, self.solution_climatology_name, self.measurement_climatology_name)
            operations += [measurement_climatology, solution_climatology]
            if iteration_index==(self.n_iterations-1):
                solution_prior_climatology =self.fill_single_prior_solution_operation_space_time(iteration_index, 0, self.solution_climatology_command, self.solution_climatology_name, self.measurement_climatology_name)
                operations += [solution_prior_climatology]
                
            # Filling large scale
            measurement_large_scale = self.fill_batch_measurement_operation_space_time(iteration_index, 1, self.measurement_large_scale_command, self.measurement_large_scale_name)
            solution_large_scale = self.fill_single_posterior_solution_operation_space_time(iteration_index, 1, self.solution_large_scale_command, self.solution_large_scale_name, self.measurement_large_scale_name)
            operations += [measurement_large_scale, solution_large_scale]
            if iteration_index==(self.n_iterations-1):
                solution_prior_large_scale = self.fill_single_prior_solution_operation_space_time(iteration_index, 1, self.solution_large_scale_command, self.solution_large_scale_name, self.measurement_large_scale_name)
                operations += [solution_prior_large_scale]
            
            # Filling local scale
            # meas_sol_local = self.fill_batch_measurement_operation_space(iteration_index, 2) # Batch operation problematic here. Prohibited by apparent memory leak in calls to Pardiso.
            meas_sol_local = self.fill_single_operation_space(iteration_index, 2)
            operations += [meas_sol_local]
        
        grid_operation = self.fill_batch_output_gridding()
        operations += [grid_operation]
        
        final_dictionary = {'python_class' : 'eumopps.catalogue.catalogue.Catalogue', 'operations':operations }
        return final_dictionary
    
    """
    
    Single output operation fillers
    
    """
    
    
    def fill_single_operation_space(self, iteration_index, component_index):
        """Fill an operation for local space model at a given iteration"""
                   
        new_dict = copy.deepcopy(self.operation_dict)
       
        command = self.local_scale_command+'_'+str(iteration_index)
        name = self.solution_local_scale_name+'_'+str(iteration_index)
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'space_solution')
        self.fill_secondary_keys(new_dict, iteration_index, component_index, name)
        if iteration_index==(self.n_iterations-1):
            #self.add_uncertainties_information(new_dict, component_index, 'pickle', 'space_uncertainty')
            sample_name = self.solution_local_scale_name+'_sample_'+str(iteration_index)
            self.add_sample_information(new_dict, component_index, sample_name, 'pickle', 'space_sample', 'single')
            
            prior_sample_name = self.solution_local_scale_name+'_prior_sample_'+str(iteration_index)
            self.add_prior_sample_information(new_dict, component_index, prior_sample_name, 'pickle', 'space_sample', 'single')
         
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.process_inputs'
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['time_index'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisStepIndex" }
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['statefilename_write'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 
                                                                                                            'datasetname' : name}
                                                                                                             
        return new_dict
    
    def fill_single_measurement_operation_space_time(self, iteration_index, component_index, measurement_command, measurement_name):
        """Fill a measurement operation for local spacetime model at a given iteration"""
       
        new_dict = copy.deepcopy(self.operation_dict)
        command = measurement_command+'_'+str(iteration_index)
        name = measurement_name+'_'+str(iteration_index)
       
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'measurement')
        self.fill_secondary_keys(new_dict, iteration_index, component_index, name)
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.process_inputs'
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['measurement_time_index_write'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisStepIndex" }
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['measurementfilename_write'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 
                                                                                                                  'datasetname': name }
        return new_dict                   
               
    def fill_single_solution_operation_space_time(self, iteration_index, component_index, solution_command, solution_name, measurement_name):
        """Fill a solution operation for local spacetime model at a given iteration"""
        new_dict = copy.deepcopy(self.operation_dict)
        command = solution_command+'_'+str(iteration_index)
        name = solution_name+'_'+str(iteration_index)
        meas_name = measurement_name+'_'+str(iteration_index)
       
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'spacetime_solution')
        
        if iteration_index==(self.n_iterations-1):
            #self.add_uncertainties_information(new_dict, component_index, 'pickle', 'spacetime_uncertainty')
            sample_name = solution_name+'_sample_'+str(iteration_index)
            self.add_sample_information(new_dict, component_index, sample_name, 'pickle', 'spacetime_sample', 'single')
            
            prior_sample_name = solution_name+'_prior_sample_'+str(iteration_index)
            self.add_prior_sample_information(new_dict, component_index, prior_sample_name, 'pickle', 'spacetime_sample', 'single')
            
        new_dict['runmodule']['component_index'] = component_index
        new_dict['step']={ 'python_class' : 'eumopps.catalogue.step.StepOnce'}
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.solve'
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['measurementfilelist_read'] = { 'python_class' : 'eumopps.catalogue.placeholder.InputFileList', 
                                                                                                                 'datasetname' : meas_name }
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['statefilename_write'] = { 'python_class': 'eumopps.catalogue.placeholder.OutputFile', 
                                                                                                            'datasetname' : name }
        for to_remove in ['inputsources', 'time_index', 'time_indices']:                                                                                                   
            if to_remove in new_dict['runmodule']:
                del new_dict['runmodule'][to_remove]
                 
        return new_dict
    
    def fill_single_posterior_solution_operation_space_time(self, iteration_index, component_index, solution_command, solution_name, measurement_name):
        """Fill a solution operation for local spacetime model at a given iteration"""
        new_dict = copy.deepcopy(self.operation_dict)
        command = solution_command+'_'+str(iteration_index)
        name = solution_name+'_'+str(iteration_index)
        meas_name = measurement_name+'_'+str(iteration_index)
       
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'spacetime_solution')
        
        if iteration_index==(self.n_iterations-1):
            #self.add_uncertainties_information(new_dict, component_index, 'pickle', 'spacetime_uncertainty')
            sample_name = solution_name+'_sample_'+str(iteration_index)
            self.add_sample_information(new_dict, component_index, sample_name, 'pickle', 'spacetime_sample', 'single')
            
        new_dict['runmodule']['component_index'] = component_index
        new_dict['step']={ 'python_class' : 'eumopps.catalogue.step.StepOnce'}
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.solve'
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['measurementfilelist_read'] = { 'python_class' : 'eumopps.catalogue.placeholder.InputFileList', 
                                                                                                                 'datasetname' : meas_name }
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['statefilename_write'] = { 'python_class': 'eumopps.catalogue.placeholder.OutputFile', 
                                                                                                            'datasetname' : name }
        for to_remove in ['inputsources', 'time_index', 'time_indices']:                                                                                                   
            if to_remove in new_dict['runmodule']:
                del new_dict['runmodule'][to_remove]
                 
        return new_dict
    
    def fill_single_prior_solution_operation_space_time(self, iteration_index, component_index, solution_command, solution_name, measurement_name):
        """Fill a solution operation for local spacetime model at a given iteration"""
        new_dict = copy.deepcopy(self.operation_dict)
        command = solution_command+'_prior_sample_'+str(iteration_index)
        name = solution_name+'_prior_sample_'+str(iteration_index)
       
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'spacetime_solution')
        
        if iteration_index==(self.n_iterations-1):
            #self.add_uncertainties_information(new_dict, component_index, 'pickle', 'spacetime_uncertainty')
            
            prior_sample_name = solution_name+'_prior_sample_'+str(iteration_index)
            self.add_prior_sample_information(new_dict, component_index, prior_sample_name, 'pickle', 'spacetime_sample', 'single')
            
        new_dict['runmodule']['component_index'] = component_index
        new_dict['step']={ 'python_class' : 'eumopps.catalogue.step.StepOnce'}
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.solve'
        
        for to_remove in ['inputsources', 'time_index', 'time_indices']:                                                                                                   
            if to_remove in new_dict['runmodule']:
                del new_dict['runmodule'][to_remove]
                 
        return new_dict
    
    """
    
    Multi-sub-set/mutli-output operations
    
    """
    
        
    def fill_batch_measurement_operation_space(self, iteration_index, component_index):
        """Fill an operation with multiple substeps in time for local space model at a given iteration"""
                   
        new_dict = copy.deepcopy(self.batch_operation_dict)
       
        command = self.local_scale_command+'_'+str(iteration_index)
        name = self.solution_local_scale_name+'_'+str(iteration_index)
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'space_solution')
        self.fill_secondary_keys_batch(new_dict, iteration_index, component_index, name)
        if iteration_index==(self.n_iterations-1):
            #self.add_uncertainties_information(new_dict, component_index, 'pickle', 'space_uncertainty')
            sample_name = self.solution_local_scale_name+'_sample_'+str(iteration_index)
            self.add_sample_information(new_dict, component_index, sample_name, 'pickle', 'space_sample', 'batch')
            
            prior_sample_name = self.solution_local_scale_name+'_prior_sample_'+str(iteration_index)
            print component_index
            self.add_prior_sample_information(new_dict, component_index, prior_sample_name, 'pickle', 'space_sample', 'batch')
         
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.process_inputs_batch'
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['time_list'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['statefilelist_write'] = { 'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput', 
                                                                                                                  'datasetname' : name}
        
        return new_dict
    
    def fill_batch_measurement_operation_space_time(self, iteration_index, component_index, measurement_command, measurement_name):
        """Fill a batch measurement operation for spacetime model at a given iteration"""
       
        new_dict = copy.deepcopy(self.batch_operation_dict)
        command = measurement_command+'_'+str(iteration_index)
        name = measurement_name+'_'+str(iteration_index)
       
        self.fill_primary_keys(new_dict, command, name, 'pickle', 'monthly_measurement')
        self.fill_secondary_keys_batch(new_dict, iteration_index, component_index, name)
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.process_inputs_batch'
        new_dict['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['measurementfilename_write'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 
                                                                                                                  'datasetname': name }
        return new_dict
    
    """
    
    Output gridding
    
    """
    
    def fill_output_gridding(self):
        """Fill final output gridding operation"""
       
        self.fill_primary_keys(self.output_gridding_dict, self.output_grid_command, self.output_grid_name, 'nc', 'measurement')
        self.output_gridding_dict['runmodule']['sample_size'] = FullstaceJsonFiller.OUTPUT_SAMPLE_SIZE
        
        self.output_gridding_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.output_grid'
        self.output_gridding_dict['runmodule']['outputfile'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : self.output_grid_name  }
        self.output_gridding_dict['runmodule']['climatologyfile'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : self.output_climatology_name  }
        self.output_gridding_dict['runmodule']['largescalefile'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : self.output_large_scale_name  }
        self.output_gridding_dict['runmodule']['localfile'] = { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : self.output_local_name  }
        
        for storage, solution_datasetname, uncertainty_datasetname, sample_datasetname, sample_prior_datasetname in zip(FullstaceJsonFiller.STORAGE_LIST, 
                                                                                              [self.solution_climatology_name, self.solution_large_scale_name, self.solution_local_scale_name], 
                                                                                              FullstaceJsonFiller.UNCERTAINTIES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST_PRIOR):
                                                                                                  
            print self.output_gridding_dict['runmodule'][storage].keys()
            self.output_gridding_dict['runmodule'][storage]['statefilename_read']['datasetname'] = solution_datasetname+'_'+str(self.n_iterations-1)
            #self.output_gridding_dict['runmodule'][storage]['marginal_std_filename_read']['datasetname'] = uncertainty_datasetname
            self.output_gridding_dict['runmodule'][storage]['statefilename_read']['datasetname'] = solution_datasetname+'_'+str(self.n_iterations-1)
            self.output_gridding_dict['runmodule'][storage]['sample_filename_read']['datasetname'] = solution_datasetname+'_sample_'+str(self.n_iterations-1)#sample_datasetname
            self.output_gridding_dict['runmodule'][storage]['prior_sample_filename_read']['datasetname'] = solution_datasetname+'_prior_sample_'+str(self.n_iterations-1)#sample_prior_datasetname
        
        
        # additional output datasets for components
        climatology_dataset = copy.deepcopy( self.output_gridding_dict['newdatasets'][0] )
        climatology_dataset['name'] = self.output_climatology_name
        subsets = climatology_dataset['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(self.output_climatology_name, 'nc', 'measurement')
        self.output_gridding_dict['newdatasets'].append( climatology_dataset )
        
        large_scale_dataset = copy.deepcopy( self.output_gridding_dict['newdatasets'][0] )
        large_scale_dataset['name'] = self.output_large_scale_name
        subsets = large_scale_dataset['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(self.output_large_scale_name, 'nc', 'measurement')
        self.output_gridding_dict['newdatasets'].append( large_scale_dataset )
        
        local_dataset = copy.deepcopy( self.output_gridding_dict['newdatasets'][0] )
        local_dataset['name'] = self.output_local_name
        subsets = local_dataset['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(self.output_local_name, 'nc', 'measurement')
        self.output_gridding_dict['newdatasets'].append( local_dataset )

    def fill_batch_output_gridding(self):
        """Fill final output gridding operation"""
       
        new_dict = copy.deepcopy(self.batch_operation_dict)
       
        self.fill_primary_keys(new_dict, self.output_grid_command, self.output_grid_name, 'nc', 'measurement')
        new_dict['runmodule']['sample_size'] = FullstaceJsonFiller.OUTPUT_SAMPLE_SIZE
        
        new_dict['runmodule']['python_function'] = 'eustace.analysis.advanced_standard.examples.moving_climatology.output_grid_batch'
        new_dict['runmodule']['outputfiles'] = { 'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput', 'datasetname' : self.output_grid_name  }
        new_dict['runmodule']['climatologyfiles'] = { 'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput', 'datasetname' : self.output_climatology_name  }
        new_dict['runmodule']['largescalefiles'] = { 'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput', 'datasetname' : self.output_large_scale_name  }
        new_dict['runmodule']['localfiles'] = { 'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput', 'datasetname' : self.output_local_name  }
        
        for storage, solution_datasetname, uncertainty_datasetname, sample_datasetname, sample_prior_datasetname in zip(FullstaceJsonFiller.STORAGE_LIST, 
                                                                                              [self.solution_climatology_name, self.solution_large_scale_name, self.solution_local_scale_name], 
                                                                                              FullstaceJsonFiller.UNCERTAINTIES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST,
                                                                                              FullstaceJsonFiller.SAMPLES_LIST_PRIOR):
                                                                                                  
            print new_dict['runmodule'][storage].keys()
            
            if solution_datasetname != self.solution_local_scale_name:
                new_dict['runmodule'][storage]['statefilename_read'] = {'python_class' : 'eumopps.catalogue.placeholder.InputFile',
                                                                        'datasetname': solution_datasetname+'_'+str(self.n_iterations-1)}
                new_dict['runmodule'][storage]['sample_filename_read'] = {'python_class' : "eumopps.catalogue.placeholder.InputFile",
                                                                          'datasetname': solution_datasetname+'_sample_'+str(self.n_iterations-1)}
                new_dict['runmodule'][storage]['prior_sample_filename_read'] = {'python_class': 'eumopps.catalogue.placeholder.InputFile',
                                                                                'datasetname':   solution_datasetname+'_prior_sample_'+str(self.n_iterations-1)}
            else:
                new_dict['runmodule'][storage]['statefilelist_read'] =         {"python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays",
                                                                                                 "datasetname" : solution_datasetname+'_'+str(self.n_iterations-1),
                                                                                                 "missing_data":"allowed"
                                                                                                }
                new_dict['runmodule'][storage]['sample_filelist_read'] =       {"python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays",
                                                                                                 "datasetname" : solution_datasetname+'_sample_'+str(self.n_iterations-1),
                                                                                                 "missing_data":"allowed"
                                                                                                }
                new_dict['runmodule'][storage]['prior_sample_filelist_read'] = {"python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays",
                                                                                                 "datasetname" : solution_datasetname+'_prior_sample_'+str(self.n_iterations-1),
                                                                                                 "missing_data":"allowed"
                                                                                                }
                new_dict['runmodule'][storage]['time_list'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }


        
        # additional output datasets for components
        climatology_dataset = copy.deepcopy( new_dict['newdatasets'][0] )
        climatology_dataset['name'] = self.output_climatology_name
        subsets = climatology_dataset['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(self.output_climatology_name, 'nc', 'measurement')
        new_dict['newdatasets'].append( climatology_dataset )
        
        large_scale_dataset = copy.deepcopy( new_dict['newdatasets'][0] )
        large_scale_dataset['name'] = self.output_large_scale_name
        subsets = large_scale_dataset['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(self.output_large_scale_name, 'nc', 'measurement')
        new_dict['newdatasets'].append( large_scale_dataset )
        
        local_dataset = copy.deepcopy( new_dict['newdatasets'][0] )
        local_dataset['name'] = self.output_local_name
        subsets = local_dataset['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(self.output_local_name, 'nc', 'measurement')
        new_dict['newdatasets'].append( local_dataset )
        
        # remove unnecessary keys
        for to_remove in ['inputsources', 'time_index']:                                                                                                   
            if to_remove in new_dict['runmodule']:
                del new_dict['runmodule'][to_remove]
                 
        new_dict['runmodule']['time_indices'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }
        
        # the string form of the date will be computed within the called function instead
        #new_dict['runmodule']['processdates'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayTimes" }
        
        return new_dict

    """
    
    Uncertainty information
    
    """

    def add_uncertainties_information(self, dictionary, component_index, file_format, name_format, mode):
        """Add information about computation of uncertainties"""
       
        dataset_name = FullstaceJsonFiller.UNCERTAINTIES_LIST[component_index]
        
        if mode == 'single':
            dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['marginal_std_filename_write'] = { "python_class": "eumopps.catalogue.placeholder.OutputFile", 
                                                                                                                                "datasetname" : dataset_name }
        elif mode == 'batch':            
            dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['marginal_std_filelist_write'] = { "python_class": "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput", 
                                                                                                                                "datasetname" : dataset_name }
        else:
            raise ValueError('Unrecognised processing mode:'+mode)
        
        dictionary['runmodule']['compute_uncertainties'] = 1
        dictionary['runmodule']['method'] = 'APPROXIMATED'
        datasets = copy.deepcopy(dictionary['newdatasets'][0])
        datasets['name'] = dataset_name
        subsets = datasets['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(dataset_name, file_format, name_format)
        dictionary['newdatasets'].append(datasets)
        
        
    def add_sample_information(self, dictionary, component_index, dataset_name, file_format, name_format, mode):
        """Add information about computation of samples"""
        
        #dataset_name = FullstaceJsonFiller.SAMPLES_LIST[component_index]
        if mode == 'single':
            dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['sample_filename_write'] = { "python_class": "eumopps.catalogue.placeholder.OutputFile", 
                                                                                                                                "datasetname" : dataset_name }
        elif mode == 'batch':
            dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['sample_filelist_write'] = { "python_class": "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput", 
                                                                                                                                "datasetname" : dataset_name }

        dictionary['runmodule']['compute_sample'] = 1
        dictionary['runmodule']['sample_size'] = FullstaceJsonFiller.STATE_SAMPLE_SIZE
        
        datasets = copy.deepcopy(dictionary['newdatasets'][0])
        datasets['name'] = dataset_name
        subsets = datasets['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(dataset_name, file_format, name_format)
        dictionary['newdatasets'].append(datasets)

    def add_prior_sample_information(self, dictionary, component_index, dataset_name, file_format, name_format, mode):
        """Add information about computation of samples"""
        print dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]
        #dataset_name = FullstaceJsonFiller.SAMPLES_LIST_PRIOR[component_index]
        if mode == 'single':
            dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['prior_sample_filename_write'] = { "python_class": "eumopps.catalogue.placeholder.OutputFile", 
                                                                                                                                "datasetname" : dataset_name }
        elif mode == 'batch':            
            dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[component_index]]['prior_sample_filelist_write'] = { "python_class": "eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDaysOutput", 
                                                                                                                                "datasetname" : dataset_name }

        dictionary['runmodule']['compute_prior_sample'] = 1
        dictionary['runmodule']['sample_size'] = FullstaceJsonFiller.STATE_SAMPLE_SIZE
        
        # check whether the dataset has already been added - work around for fill primary keys adding datasets for dedicated prior sampling operation for space-time components
        if not any([dataset['name'] == dataset_name for dataset in dictionary['newdatasets']]):
            datasets = copy.deepcopy(dictionary['newdatasets'][0])
            datasets['name'] = dataset_name
            subsets = datasets['subsets'][0]
            subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(dataset_name, file_format, name_format)
            dictionary['newdatasets'].append(datasets)

    """
    
    Additional keyword filling
    
    """

    def fill_primary_keys(self, dictionary, command_name, dataset_name, file_format, name_format):
        """Fill primary, common keys
        
        Main analysis settings, date range and output file patterns for new data sets
        
        """
        
        dictionary['name'] = command_name
        dictionary['runmodule']['insitu_biases'] = self.land_biases
        dictionary['runmodule']['global_biases'] = self.global_biases
        dictionary['runmodule']['compute_uncertainties'] = 0
        dictionary['runmodule']['method'] = 'EXACT'
        dictionary['runmodule']['compute_sample'] = 0
        dictionary['runmodule']['sample_size'] = 0
        dictionary['runmodule']['compute_prior_sample'] = 0
        
        for point, value in zip(['start', 'end'], [self.start_date, self.end_date]):
            dictionary['step'][point] = value
        
        datasets = dictionary['newdatasets'][0]
        datasets['name'] = dataset_name
        subsets = datasets['subsets'][0]
        subsets['layout']['patterns'] = FullstaceJsonFiller.return_output_patterns(dataset_name, file_format, name_format)

    def fill_secondary_keys(self, dictionary, iteration_index, component_index, name):
        """Fill secondary, less common keys
        
        Component reading for other components
        
        """
       
        dictionary['runmodule']['component_index'] = component_index                                                                       
        datasets = [self.solution_climatology_name, self.solution_large_scale_name, self.solution_local_scale_name ]
       
        # local scale solution is produced on daily basis
        if component_index != 2:
            dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[2]]['time_index'] = {'python_class': 'eustace.analysis.advanced_standard.examples.placeholder.AnalysisStepIndex'}
        for shifter in range(1,3):
            index = component_index-shifter
           
            if index >= 0 and iteration_index>=0:
                dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[index]]['statefilename_read'] = {'python_class' : 'eumopps.catalogue.placeholder.InputFile',
                                                                             'datasetname':datasets[index]+'_'+str(iteration_index)}
            elif index < 0 and iteration_index>0:
                dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[index]]['statefilename_read'] = {'python_class' : 'eumopps.catalogue.placeholder.InputFile',
                                                                                                             'datasetname':datasets[index]+'_'+str(iteration_index-1)}

    def fill_secondary_keys_batch(self, dictionary, iteration_index, component_index, name):
        """Fill secondary, less common keys
        
        Component reading for other components
        
        """
       
        dictionary['runmodule']['component_index'] = component_index                                                                       
        datasets = [self.solution_climatology_name, self.solution_large_scale_name, self.solution_local_scale_name ]
       
        for storage_index in range(3):
            
            if storage_index < component_index and iteration_index>=0:
                if storage_index != 2:
                    dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[storage_index]]['statefilename_read'] = {'python_class' : 'eumopps.catalogue.placeholder.InputFile',
                                                                                                              'datasetname':datasets[storage_index]+'_'+str(iteration_index)}
                else:
                    dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[storage_index]]['statefilelist_read'] = {'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays',
                                                                                                              'datasetname':datasets[storage_index]+'_'+str(iteration_index), }
                                                                                                              #'time_list': { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" } }
                    dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[storage_index]]['time_list'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }
            elif storage_index >= component_index and iteration_index>0:
                if storage_index != 2:
                    dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[storage_index]]['statefilename_read'] = {'python_class' : 'eumopps.catalogue.placeholder.InputFile',
                                                                                                             'datasetname':datasets[storage_index]+'_'+str(iteration_index-1)}
                else:
                    dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[storage_index]]['statefilelist_read'] = {'python_class' : 'eustace.analysis.advanced_standard.examples.placeholder.AnnualBatchDays',
                                                                                                              'datasetname':datasets[storage_index]+'_'+str(iteration_index-1), }
                                                                                                              #'time_list': { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" } }
                    dictionary['runmodule'][FullstaceJsonFiller.STORAGE_LIST[storage_index]]['time_list'] = { "python_class" : "eustace.analysis.advanced_standard.examples.placeholder.AnalysisAnnualBatchDayIndices" }
                    

    @staticmethod
    def return_output_patterns(datasetname, file_format, name_format):
        if (name_format=='spacetime_solution') or (name_format=='spacetime_uncertainty') or (name_format=='spacetime_sample'):
            return [ datasetname+'/', datasetname+'.'+file_format ]
        elif (name_format=='measurement') or (name_format=='space_solution') or (name_format=='space_uncertainty') or (name_format=='space_sample'):
            return [ datasetname+'/%Y/', datasetname+'_%Y%m%d.'+file_format ]
        elif (name_format=='annual_measurement'):
            return [ datasetname+'/%Y/', datasetname+'_%Y.'+file_format ]
        elif (name_format=='monthly_measurement'):
            return [ datasetname+'/%Y/', datasetname+'_%Y%m.'+file_format ]
        elif name_format == 'region_optimization':
            return [ datasetname+'/optimization/', datasetname+'_%i.%i.%i.'+file_format ]
        else:
            raise ValueError('Invalid name format flag')
       
    @staticmethod   
    def check_none(descriptor):
        """Ensures that all relevant fields in the descriptor have been filled"""
        if isinstance(descriptor, dict):
            for key, value in descriptor.items():
                FullstaceJsonFiller.check_none(value)
        elif isinstance(descriptor, list):
            return [ FullstaceJsonFiller.check_none(item) for item in descriptor ]
        else:
            if descriptor==None:
                message = 'Field {} has not been filled! Cannot create json descriptor'.format(descriptor)
                raise ValueError(message)
