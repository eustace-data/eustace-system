"""Test the implementation of classes devised for filling json descriptors"""

import json
import numpy
import os
import tempfile
import unittest

from eustace.analysis.advanced_standard.examples.jsonfiller import FullstaceJsonFiller
from eustace.outputformats import definitions

class TestFullstaceJsonFiller(unittest.TestCase):

    def setUp(self):
      
      self.test_dictionaries = [{'A': 1434}, {'B': 'cat', 'please':456}]
      
      test_template_files = [tempfile.NamedTemporaryFile(prefix=str(index)+'.json', delete=False) for index in range(len(self.test_dictionaries))]
      self.test_template_filenames = [test_template_file.name for test_template_file in test_template_files]
      
      for dictionary, test_template_file in zip(self.test_dictionaries, test_template_files):
          json.dump(dictionary, test_template_file)
      
      self.dictionaries_with_None = [{'A': None}, 
                     {'B': [{'cat':3, 'please':None}, {'cato':2}]}, 
                     [{'B': {'cat':3, 'please':4}}, {'B': [{'cat':3, 'please':None}, {'dogo':1}]}]]
                     
      test_template_file = tempfile.NamedTemporaryFile(prefix='output_grid.json', delete=False)
      
      template_dictionary = {"python_class" : "eumopps.catalogue.operation.Operation",
                  "name" : None,
                  "runmodule" : {
                    "storage_climatology":{"statefilename_read" : {"datasetname" : None, "missing_data":"allowed" },
                               "marginal_std_filename_read" : {"datasetname" : None, "missing_data":"allowed" },
                               "sample_filename_read" : {"datasetname" : None, "missing_data":"allowed" }},
                    "storage_large_scale" :{"statefilename_read" : {"datasetname" : None, "missing_data":"allowed" },
                                    "marginal_std_filename_read" : {"datasetname" : None, "missing_data":"allowed" },
                                    "sample_filename_read" : {"datasetname" : None, "missing_data":"allowed" }},
                    "storage_local" :{"statefilename_read" : {"datasetname" : None, "missing_data":"allowed" },
                                          "marginal_std_filename_read" : {"datasetname" : None, "missing_data":"allowed" },
                                          "sample_filename_read" : {"datasetname" : None, "missing_data":"allowed" }},
                    "insitu_biases": None ,
                    "global_biases": None,
                    "compute_uncertainties": None,
                    "method": None,
                    "compute_sample": None,
                    "sample_size": None},
                  "step" : { "python_class" : "eumopps.catalogue.step.StepDaily", "start" : None, "end" : None },
                  "newdatasets": [{"name"  : None,
                          "subsets" : [{"layout" :{"patterns" : None}}]}]}
                                  
      self.test_template_file = test_template_file.name
      json.dump(template_dictionary, test_template_file)
      
      # To uncomment when more detailed tests on descriptor generation are required
      #self.original_descriptor_file = 'test/test_descriptor.json'
      
    def tearDown(self):
    
      for name in self.test_template_filenames:
          os.remove(name)
      
      os.remove(self.test_template_file)

    def test_init(self):
      current_dir = os.getcwd()
      template = 'test.json'
      date = '19870405'
      bias = 1
      filler = FullstaceJsonFiller(3, 
                   current_dir, current_dir, current_dir,
                   template, template, template,
                   date, date,
                   1, 1,
                   'A', 'B', 'C', 'D', 
                   'E', 'F', 'G', 'H',
                   'I', 'J', 
                   'K', 'L',
                   'M', 'N', 'O')
                   
      self.assertEqual(3, filler.n_iterations)

      self.assertEqual(current_dir, filler.inpath)
      self.assertEqual(current_dir, filler.outpath)

      self.assertEqual(template, filler.operation_template)
      self.assertEqual(template, filler.output_grid_template)

      self.assertEqual(date+'000000', filler.start_date)
      self.assertEqual(date+'000000', filler.end_date)

      self.assertEqual(bias, filler.land_biases)
      self.assertEqual(bias, filler.global_biases)
      
      self.assertEqual('A', filler.measurement_climatology_command)
      self.assertEqual('B', filler.measurement_climatology_name)
      self.assertEqual('C', filler.solution_climatology_command)
      self.assertEqual('D', filler.solution_climatology_name)

      self.assertEqual('E', filler.measurement_large_scale_command)
      self.assertEqual('F', filler.measurement_large_scale_name)
      self.assertEqual('G', filler.solution_large_scale_command)
      self.assertEqual('H', filler.solution_large_scale_name)
                       
      self.assertEqual('I', filler.local_scale_command)
      self.assertEqual('J', filler.solution_local_scale_name)

      self.assertEqual('K', filler.output_grid_command)
      self.assertEqual('L', filler.output_grid_name)
      
    def test_load_templates(self):
      
      current_dir=''
      template_operation = self.test_template_filenames[0]
      template_grid = self.test_template_filenames[1]
      date = '19870405'
      bias = 1
      filler = FullstaceJsonFiller(3, 
                   current_dir, current_dir, current_dir,
                   template_operation, template_grid,
                   date, date,
                   1, 1,
                   'A', 'B', 'C', 'D', 
                   'E', 'F', 'G', 'H',
                   'I', 'J', 
                   'K', 'L',
                   'M', 'N', 'O')


      filler.load_templates()
      self.assertDictEqual(self.test_dictionaries[0], filler.operation_dict)
      self.assertDictEqual(self.test_dictionaries[1], filler.output_gridding_dict)
      
    def test_fill_output_gridding(self):
      
      current_dir=''
      template_operation = self.test_template_filenames[0]
      template_grid = self.test_template_file
      start_date = '19870405'
      end_date = '19880405'
      bias = 1
      filler = FullstaceJsonFiller(3, 
                   current_dir, current_dir, current_dir,
                   template_operation, template_grid,
                   start_date, end_date,
                   1, 1,
                   'A', 'B', 'C', 'D', 
                   'E', 'F', 'G', 'H',
                   'I', 'J', 
                   'output_grid', 'eustace_infilled',
                   'M', 'N', 'O')
                   
      filler.load_templates()
      filler.fill_output_gridding()
      
      self.assertEqual(filler.output_gridding_dict['name'], 'output_grid')
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'D_2')
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], 'H_2')
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_local']['statefilename_read']['datasetname'], 'J_2')    
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_climatology']['marginal_std_filename_read']['datasetname'], 'uncertainties_climatology')
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_large_scale']['marginal_std_filename_read']['datasetname'], 'uncertainties_large_scale')
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_local']['marginal_std_filename_read']['datasetname'], 'uncertainties_local')    
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_climatology']['sample_filename_read']['datasetname'], 'sample_climatology')
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_large_scale']['sample_filename_read']['datasetname'], 'sample_large_scale')
      self.assertEqual(filler.output_gridding_dict['runmodule']['storage_local']['sample_filename_read']['datasetname'], 'sample_local')    


      self.assertEqual(filler.output_gridding_dict['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.output_grid')
      self.assertDictEqual(filler.output_gridding_dict['runmodule']['outputfile'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'eustace_infilled'})
      self.assertEqual(filler.output_gridding_dict['runmodule']['insitu_biases'], 1)
      self.assertEqual(filler.output_gridding_dict['runmodule']['global_biases'], 1)
      self.assertEqual(filler.output_gridding_dict['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(filler.output_gridding_dict['runmodule']['compute_sample'], 0)

      self.assertEqual(filler.output_gridding_dict['step']['start'], start_date+'000000')
      self.assertEqual(filler.output_gridding_dict['step']['end'], end_date+'000000')
      self.assertEqual(filler.output_gridding_dict['newdatasets'][0]['name'], 'eustace_infilled')
      self.assertListEqual(filler.output_gridding_dict['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['eustace_infilled/%Y/', 'eustace_infilled_%Y%m%d.nc'])
      
    def test_fill_single_operation_space(self):
      current_dir=''
      template_operation = self.test_template_file
      template_grid = self.test_template_filenames[0]
      start_date = '19870405'
      end_date = '19880405'
      bias = 1
      filler = FullstaceJsonFiller(3, 
                   current_dir, current_dir,
                   template_operation, template_grid,
                   start_date, end_date,
                   1, 1,
                   'A', 'B', 'C', 'D', 
                   'E', 'F', 'G', 'H',
                   'local_input_and_solve', 'solution_local', 
                   'K', 'L',
                   'M', 'N', 'O')
                   
      filler.load_templates()
      result = filler.fill_single_operation_space(0, 2)

      self.assertEqual(result['name'], 'local_input_and_solve_0')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'D_0')
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], 'H_0')
      self.assertDictEqual(result['runmodule']['storage_local']['time_index'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })

      self.assertEqual(result['runmodule']['insitu_biases'], 1)
      self.assertEqual(result['runmodule']['global_biases'], 1)

      self.assertEqual(result['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(result['runmodule']['method'], 'EXACT')
      
      self.assertEqual(result['runmodule']['compute_sample'], 0)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])
      
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      self.assertEqual(result['runmodule']['component_index'], 2)
      self.assertDictEqual(result['runmodule']['storage_local']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'solution_local_0'})

      self.assertEqual(result['step']['start'], start_date+'000000')
      self.assertEqual(result['step']['end'], end_date+'000000')
      self.assertEqual(len(result['newdatasets']),1)
      self.assertEqual(result['newdatasets'][0]['name'], 'solution_local_0')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_local_0/%Y/', 'solution_local_0_%Y%m%d.pickle'])

      result = filler.fill_single_operation_space(2, 2)

      self.assertEqual(result['name'], 'local_input_and_solve_2')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'D_2')
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], 'H_2')
      self.assertDictEqual(result['runmodule']['storage_local']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'solution_local_2'})
      self.assertDictEqual(result['runmodule']['storage_local']['marginal_std_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'uncertainties_local'})
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 1)
      self.assertEqual(result['runmodule']['method'], 'APPROXIMATED')
      
      self.assertEqual(result['runmodule']['compute_sample'], 1)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])
      
      self.assertEqual(len(result['newdatasets']), 3)
      self.assertEqual(result['newdatasets'][0]['name'], 'solution_local_2')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_local_2/%Y/', 'solution_local_2_%Y%m%d.pickle'])
      self.assertEqual(result['newdatasets'][1]['name'], 'uncertainties_local')
      self.assertListEqual(result['newdatasets'][1]['subsets'][0]['layout']['patterns'], ['uncertainties_local/%Y/', 'uncertainties_local_%Y%m%d.pickle'])
      self.assertEqual(result['newdatasets'][2]['name'], 'sample_local')
      self.assertListEqual(result['newdatasets'][2]['subsets'][0]['layout']['patterns'], ['sample_local/%Y/', 'sample_local_%Y%m%d.pickle'])


    def test_fill_single_measurement_operation_space_time(self):
      current_dir=''
      template_operation = self.test_template_file
      template_grid = self.test_template_filenames[0]
      start_date = '19870405'
      end_date = '19880405'
      bias = 1
      filler = FullstaceJsonFiller(3, 
                   current_dir, current_dir,
                   template_operation, template_grid,
                   start_date, end_date,
                   1, 1,
                   'climatology_input', 'measurement_climatology', 'C', 'D', 
                   'E', 'F', 'G', 'H',
                   'I', 'J', 
                   'K', 'L',
                   'M', 'N', 'O')
                   
      filler.load_templates()
      result = filler.fill_single_measurement_operation_space_time(0, 0, filler.measurement_climatology_command, filler.measurement_climatology_name)

      self.assertEqual(result['name'], 'climatology_input_0')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertDictEqual(result['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      
      self.assertEqual(result['runmodule']['insitu_biases'], 1)
      self.assertEqual(result['runmodule']['global_biases'], 1)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(result['runmodule']['method'], 'EXACT')
      
      self.assertEqual(result['runmodule']['compute_sample'], 0)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])
      
      self.assertEqual(result['runmodule']['component_index'], 0)
      self.assertDictEqual(result['runmodule']['storage_climatology']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_climatology_0'})
      self.assertDictEqual(result['runmodule']['storage_climatology']['measurement_time_index_write'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })
      self.assertEqual(result['step']['start'], start_date+'000000')
      self.assertEqual(result['step']['end'], end_date+'000000')
      self.assertEqual(result['newdatasets'][0]['name'], 'measurement_climatology_0')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_climatology_0/%Y/', 'measurement_climatology_0_%Y%m%d.pickle'])

      result = filler.fill_single_measurement_operation_space_time(51, 0, filler.measurement_climatology_command, filler.measurement_climatology_name)

      self.assertEqual(result['name'], 'climatology_input_51')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], 'H_50')
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], 'J_50')
      self.assertDictEqual(result['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      self.assertEqual(result['runmodule']['component_index'], 0)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(result['runmodule']['method'], 'EXACT')
      
      self.assertEqual(result['runmodule']['compute_sample'], 0)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])
            
      self.assertDictEqual(result['runmodule']['storage_climatology']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_climatology_51'})
      self.assertEqual(result['newdatasets'][0]['name'], 'measurement_climatology_51')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_climatology_51/%Y/', 'measurement_climatology_51_%Y%m%d.pickle'])

      filler = FullstaceJsonFiller(3, 
                   current_dir, current_dir,
                   template_operation, template_grid,
                   start_date, end_date,
                   1, 1,
                   'A', 'C', 'C', 'D', 
                   'large_scale_input', 'measurement_large_scale', 'G', 'H',
                   'I', 'J', 
                   'K', 'L',
                   'M', 'N', 'O')
                   
      filler.load_templates()
      result = filler.fill_single_measurement_operation_space_time(0, 1, filler.measurement_large_scale_command, filler.measurement_large_scale_name)

      self.assertEqual(result['name'], 'large_scale_input_0')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'D_0')
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertDictEqual(result['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      
      self.assertEqual(result['runmodule']['insitu_biases'], 1)
      self.assertEqual(result['runmodule']['global_biases'], 1)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(result['runmodule']['method'], 'EXACT')
      
      self.assertEqual(result['runmodule']['compute_sample'], 0)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(result['runmodule']['component_index'], 1)
      self.assertDictEqual(result['runmodule']['storage_large_scale']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_large_scale_0'})
      self.assertDictEqual(result['runmodule']['storage_large_scale']['measurement_time_index_write'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })
      self.assertEqual(result['step']['start'], start_date+'000000')
      self.assertEqual(result['step']['end'], end_date+'000000')
      self.assertEqual(result['newdatasets'][0]['name'], 'measurement_large_scale_0')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_large_scale_0/%Y/', 'measurement_large_scale_0_%Y%m%d.pickle'])

      result = filler.fill_single_measurement_operation_space_time(14, 1, filler.measurement_large_scale_command, filler.measurement_large_scale_name)
      self.assertEqual(result['name'], 'large_scale_input_14')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'D_14')
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], 'J_13')
      self.assertDictEqual(result['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      self.assertEqual(result['runmodule']['component_index'], 1)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(result['runmodule']['method'], 'EXACT')
      
      self.assertEqual(result['runmodule']['compute_sample'], 0)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertDictEqual(result['runmodule']['storage_large_scale']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_large_scale_14'})
      self.assertEqual(result['newdatasets'][0]['name'], 'measurement_large_scale_14')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_large_scale_14/%Y/', 'measurement_large_scale_14_%Y%m%d.pickle'])

    def test_fill_single_solution_operation_space_time(self):
      current_dir=''
      template_operation = self.test_template_file
      template_grid = self.test_template_filenames[0]
      start_date = '19870405'
      end_date = '19880405'
      bias = 1
      filler = FullstaceJsonFiller(52, 
                   current_dir, current_dir,
                   template_operation, template_grid,
                   start_date, end_date,
                   1, 0,
                   'A', 'measurement_climatology', 'climatology_solve', 'solution_climatology', 
                   'E', 'F', 'G', 'H',
                   'I', 'J', 
                   'K', 'L',
                   'M', 'N', 'O')
                   
      filler.load_templates()
      result = filler.fill_single_solution_operation_space_time(0, 0, filler.solution_climatology_command, filler.solution_climatology_name, filler.measurement_climatology_name)

      self.assertEqual(result['name'], 'climatology_solve_0')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      
      self.assertEqual(result['runmodule']['insitu_biases'], 1)
      self.assertEqual(result['runmodule']['global_biases'], 0)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(result['runmodule']['method'], 'EXACT')

      self.assertEqual(result['runmodule']['compute_sample'], 0)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(result['runmodule']['component_index'], 0)
      self.assertDictEqual(result['runmodule']['storage_climatology']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'solution_climatology_0'})
      self.assertDictEqual(result['runmodule']['storage_climatology']['measurementfilelist_read'], { 'python_class' : 'eumopps.catalogue.placeholder.InputFileList','datasetname' : 'measurement_climatology_0'})
      self.assertDictEqual(result['step'], { 'python_class' : 'eumopps.catalogue.step.StepOnce'})
      self.assertEqual(len(result['newdatasets']), 1)
      self.assertEqual(result['newdatasets'][0]['name'], 'solution_climatology_0')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_climatology_0/', 'solution_climatology_0.pickle'])

      result = filler.fill_single_solution_operation_space_time(51, 0, filler.solution_climatology_command, filler.solution_climatology_name, filler.measurement_climatology_name)

      self.assertEqual(result['name'], 'climatology_solve_51')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      self.assertEqual(result['runmodule']['component_index'], 0)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 1)
      self.assertEqual(result['runmodule']['method'], 'APPROXIMATED')

      self.assertEqual(result['runmodule']['compute_sample'], 1)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertDictEqual(result['runmodule']['storage_climatology']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'solution_climatology_51'})
      self.assertDictEqual(result['runmodule']['storage_climatology']['marginal_std_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'uncertainties_climatology'})
      self.assertDictEqual(result['runmodule']['storage_climatology']['sample_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'sample_climatology'})

      self.assertDictEqual(result['runmodule']['storage_climatology']['measurementfilelist_read'], { 'python_class' : 'eumopps.catalogue.placeholder.InputFileList','datasetname' : 'measurement_climatology_51'})
      self.assertEqual(len(result['newdatasets']), 3)
      self.assertEqual(result['newdatasets'][0]['name'], 'solution_climatology_51')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_climatology_51/', 'solution_climatology_51.pickle'])
      self.assertEqual(result['newdatasets'][1]['name'], 'uncertainties_climatology')
      self.assertListEqual(result['newdatasets'][1]['subsets'][0]['layout']['patterns'], ['uncertainties_climatology/', 'uncertainties_climatology.pickle'])
      self.assertEqual(result['newdatasets'][2]['name'], 'sample_climatology')
      self.assertListEqual(result['newdatasets'][2]['subsets'][0]['layout']['patterns'], ['sample_climatology/', 'sample_climatology.pickle'])


      filler = FullstaceJsonFiller(15, 
                   current_dir, current_dir,
                   template_operation, template_grid,
                   start_date, end_date,
                   0, 1,
                   'A', 'C', 'C', 'D', 
                   'E', 'measurement_large_scale', 'large_scale_solve', 'solution_large_scale',
                   'I', 'J', 
                   'K', 'L',
                   'M', 'N', 'O')
                   
      filler.load_templates()
      result = filler.fill_single_solution_operation_space_time(0, 1, filler.solution_large_scale_command, filler.solution_large_scale_name,filler.measurement_large_scale_name)

      self.assertEqual(result['name'], 'large_scale_solve_0')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      
      self.assertEqual(result['runmodule']['insitu_biases'], 0)
      self.assertEqual(result['runmodule']['global_biases'], 1)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(result['runmodule']['method'], 'EXACT')
      
      self.assertEqual(result['runmodule']['compute_sample'], 0)
      self.assertEqual(result['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(result['runmodule']['component_index'], 1)
      self.assertDictEqual(result['runmodule']['storage_large_scale']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'solution_large_scale_0'})
      self.assertDictEqual(result['runmodule']['storage_large_scale']['measurementfilelist_read'], { 'python_class' : 'eumopps.catalogue.placeholder.InputFileList','datasetname' : 'measurement_large_scale_0'})
      self.assertDictEqual(result['step'], { 'python_class' : 'eumopps.catalogue.step.StepOnce'})
      self.assertEqual(len(result['newdatasets']), 1)
      self.assertEqual(result['newdatasets'][0]['name'], 'solution_large_scale_0')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_large_scale_0/', 'solution_large_scale_0.pickle'])

      result = filler.fill_single_solution_operation_space_time(14, 1, filler.solution_large_scale_command, filler.solution_large_scale_name,filler.measurement_large_scale_name)
      self.assertEqual(result['name'], 'large_scale_solve_14')
  
      self.assertEqual(result['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(result['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      self.assertEqual(result['runmodule']['component_index'], 1)
      
      self.assertEqual(result['runmodule']['compute_uncertainties'], 1)
      self.assertEqual(result['runmodule']['method'], 'APPROXIMATED')

      self.assertDictEqual(result['runmodule']['storage_large_scale']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'solution_large_scale_14'})
      self.assertDictEqual(result['runmodule']['storage_large_scale']['marginal_std_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'uncertainties_large_scale'})
      self.assertDictEqual(result['runmodule']['storage_large_scale']['sample_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile','datasetname' : 'sample_large_scale'})
      self.assertDictEqual(result['runmodule']['storage_large_scale']['measurementfilelist_read'], { 'python_class' : 'eumopps.catalogue.placeholder.InputFileList','datasetname' : 'measurement_large_scale_14'})
      self.assertEqual(len(result['newdatasets']), 3)
      self.assertEqual(result['newdatasets'][0]['name'], 'solution_large_scale_14')
      self.assertListEqual(result['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_large_scale_14/', 'solution_large_scale_14.pickle'])
      self.assertEqual(result['newdatasets'][1]['name'], 'uncertainties_large_scale')
      self.assertListEqual(result['newdatasets'][1]['subsets'][0]['layout']['patterns'], ['uncertainties_large_scale/', 'uncertainties_large_scale.pickle'])
      self.assertEqual(result['newdatasets'][2]['name'], 'sample_large_scale')
      self.assertListEqual(result['newdatasets'][2]['subsets'][0]['layout']['patterns'], ['sample_large_scale/', 'sample_large_scale.pickle'])


    def test_fill_descriptor(self):
      current_dir=''
      template_operation = self.test_template_file
      template_grid = self.test_template_file
      start_date = '19870405'
      end_date = '19880405'
      bias = 1
      filler = FullstaceJsonFiller(3, 
                   current_dir, current_dir,
                   template_operation, template_grid,
                   start_date, end_date,
                   1, 1,
                   'climatology_input', 'measurement_climatology', 'climatology_solve', 'solution_climatology', 
                   'large_scale_input', 'measurement_large_scale', 'large_scale_solve', 'solution_large_scale',
                   'local_input_and_solve', 'solution_local', 
                   'output_grid', 'eustace_infilled',
                   'output_climatology_name', 'output_large_scale_name', 'output_local_name')
                   
      result = filler.fill_descriptor()
      
      self.assertEqual(result['python_class'], 'eumopps.catalogue.catalogue.Catalogue')
      self.assertTrue(isinstance(result['operations'], list))
      self.assertEqual(len(result['operations']), 16)
      
      # Analyzing more in detail

      # Climatology
      climatology_input_0 = result['operations'][0]
      self.assertEqual(climatology_input_0['name'], 'climatology_input_0')
  
      self.assertEqual(climatology_input_0['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_input_0['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_input_0['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertDictEqual(climatology_input_0['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(climatology_input_0['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      self.assertEqual(climatology_input_0['runmodule']['insitu_biases'], 1)
      self.assertEqual(climatology_input_0['runmodule']['global_biases'], 1)
      self.assertEqual(climatology_input_0['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(climatology_input_0['runmodule']['method'], 'EXACT')
      self.assertEqual(climatology_input_0['runmodule']['compute_sample'], 0)
      self.assertEqual(climatology_input_0['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(climatology_input_0['runmodule']['component_index'], 0)
      self.assertDictEqual(climatology_input_0['runmodule']['storage_climatology']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_climatology_0'})
      self.assertDictEqual(climatology_input_0['runmodule']['storage_climatology']['measurement_time_index_write'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })
      self.assertEqual(climatology_input_0['step']['start'], start_date+'000000')
      self.assertEqual(climatology_input_0['step']['end'], end_date+'000000')
      self.assertEqual(len(climatology_input_0['newdatasets']), 1)
      self.assertEqual(climatology_input_0['newdatasets'][0]['name'], 'measurement_climatology_0')
      self.assertListEqual(climatology_input_0['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_climatology_0/%Y/', 'measurement_climatology_0_%Y%m%d.pickle'])

      climatology_solve_0 = result['operations'][1]
      self.assertEqual(climatology_solve_0['name'], 'climatology_solve_0')
  
      self.assertEqual(climatology_solve_0['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_solve_0['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_solve_0['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_solve_0['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      self.assertEqual(climatology_solve_0['runmodule']['insitu_biases'], 1)
      self.assertEqual(climatology_solve_0['runmodule']['global_biases'], 1)
      self.assertEqual(climatology_solve_0['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(climatology_solve_0['runmodule']['method'], 'EXACT')
      self.assertEqual(climatology_solve_0['runmodule']['compute_sample'], 0)
      self.assertEqual(climatology_solve_0['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(climatology_solve_0['runmodule']['component_index'], 0)
      self.assertDictEqual(climatology_solve_0['runmodule']['storage_climatology']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'solution_climatology_0'})
      self.assertDictEqual(climatology_solve_0['step'], { 'python_class' : 'eumopps.catalogue.step.StepOnce'})
      self.assertEqual(climatology_solve_0['newdatasets'][0]['name'], 'solution_climatology_0')
      self.assertListEqual(climatology_solve_0['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_climatology_0/', 'solution_climatology_0.pickle'])


      climatology_input_2 = result['operations'][10]
      self.assertEqual(climatology_input_2['name'], 'climatology_input_2')
  
      self.assertEqual(climatology_input_2['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_input_2['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], 'solution_large_scale_1')
      self.assertEqual(climatology_input_2['runmodule']['storage_local']['statefilename_read']['datasetname'], 'solution_local_1')
      self.assertDictEqual(climatology_input_2['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(climatology_input_2['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      self.assertEqual(climatology_input_2['runmodule']['insitu_biases'], 1)
      self.assertEqual(climatology_input_2['runmodule']['global_biases'], 1)
      self.assertEqual(climatology_input_2['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(climatology_input_2['runmodule']['method'], 'EXACT')
      self.assertEqual(climatology_input_2['runmodule']['compute_sample'], 0)
      self.assertEqual(climatology_input_2['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(climatology_input_2['runmodule']['component_index'], 0)
      self.assertDictEqual(climatology_input_2['runmodule']['storage_climatology']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_climatology_2'})
      self.assertDictEqual(climatology_input_2['runmodule']['storage_climatology']['measurement_time_index_write'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })
      self.assertEqual(climatology_input_2['step']['start'], start_date+'000000')
      self.assertEqual(climatology_input_2['step']['end'], end_date+'000000')
      self.assertEqual(climatology_input_2['newdatasets'][0]['name'], 'measurement_climatology_2')
      self.assertListEqual(climatology_input_2['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_climatology_2/%Y/', 'measurement_climatology_2_%Y%m%d.pickle'])

      climatology_solve_2 = result['operations'][11]
      self.assertEqual(climatology_solve_2['name'], 'climatology_solve_2')
  
      self.assertEqual(climatology_solve_2['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_solve_2['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_solve_2['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(climatology_solve_2['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      self.assertEqual(climatology_solve_2['runmodule']['insitu_biases'], 1)
      self.assertEqual(climatology_solve_2['runmodule']['global_biases'], 1)
      self.assertEqual(climatology_solve_2['runmodule']['compute_uncertainties'], 1)
      self.assertEqual(climatology_solve_2['runmodule']['method'], 'APPROXIMATED')
      self.assertEqual(climatology_solve_2['runmodule']['compute_sample'], 1)
      self.assertEqual(climatology_solve_2['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(climatology_solve_2['runmodule']['component_index'], 0)
      self.assertDictEqual(climatology_solve_2['runmodule']['storage_climatology']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'solution_climatology_2'})
      self.assertDictEqual(climatology_solve_2['runmodule']['storage_climatology']['marginal_std_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'uncertainties_climatology'})
      self.assertDictEqual(climatology_solve_2['runmodule']['storage_climatology']['sample_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'sample_climatology'})
      self.assertDictEqual(climatology_solve_2['step'], { 'python_class' : 'eumopps.catalogue.step.StepOnce'})
      self.assertEqual(climatology_solve_2['newdatasets'][0]['name'], 'solution_climatology_2')
      self.assertListEqual(climatology_solve_2['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_climatology_2/', 'solution_climatology_2.pickle'])
      self.assertEqual(climatology_solve_2['newdatasets'][1]['name'], 'uncertainties_climatology')
      self.assertListEqual(climatology_solve_2['newdatasets'][1]['subsets'][0]['layout']['patterns'], ['uncertainties_climatology/', 'uncertainties_climatology.pickle'])
      self.assertEqual(climatology_solve_2['newdatasets'][2]['name'], 'sample_climatology')
      self.assertListEqual(climatology_solve_2['newdatasets'][2]['subsets'][0]['layout']['patterns'], ['sample_climatology/', 'sample_climatology.pickle'])

      # Large scale
      large_scale_input_0 = result['operations'][2]
      self.assertEqual(large_scale_input_0['name'], 'large_scale_input_0')
  
      self.assertEqual(large_scale_input_0['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'solution_climatology_0')
      self.assertEqual(large_scale_input_0['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_input_0['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertDictEqual(large_scale_input_0['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(large_scale_input_0['runmodule']['insitu_biases'], 1)
      self.assertEqual(large_scale_input_0['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      self.assertEqual(large_scale_input_0['runmodule']['global_biases'], 1)
      self.assertEqual(large_scale_input_0['runmodule']['component_index'], 1)
      self.assertEqual(large_scale_input_0['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(large_scale_input_0['runmodule']['method'], 'EXACT')
      self.assertEqual(large_scale_input_0['runmodule']['compute_sample'], 0)
      self.assertEqual(large_scale_input_0['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])


      self.assertDictEqual(large_scale_input_0['runmodule']['storage_large_scale']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_large_scale_0'})
      self.assertDictEqual(large_scale_input_0['runmodule']['storage_large_scale']['measurement_time_index_write'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })
      self.assertEqual(large_scale_input_0['step']['start'], start_date+'000000')
      self.assertEqual(large_scale_input_0['step']['end'], end_date+'000000')
      self.assertEqual(len(large_scale_input_0['newdatasets']), 1)
      self.assertEqual(large_scale_input_0['newdatasets'][0]['name'], 'measurement_large_scale_0')
      self.assertListEqual(large_scale_input_0['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_large_scale_0/%Y/', 'measurement_large_scale_0_%Y%m%d.pickle'])

      large_scale_solve_0 = result['operations'][3]
      self.assertEqual(large_scale_solve_0['name'], 'large_scale_solve_0')
  
      self.assertEqual(large_scale_solve_0['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_solve_0['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_solve_0['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_solve_0['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      self.assertEqual(large_scale_solve_0['runmodule']['insitu_biases'], 1)
      self.assertEqual(large_scale_solve_0['runmodule']['global_biases'], 1)
      self.assertEqual(large_scale_solve_0['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(large_scale_solve_0['runmodule']['method'], 'EXACT')
      self.assertEqual(large_scale_solve_0['runmodule']['compute_sample'], 0)
      self.assertEqual(large_scale_solve_0['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(large_scale_solve_0['runmodule']['component_index'], 1)
      self.assertDictEqual(large_scale_solve_0['runmodule']['storage_large_scale']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'solution_large_scale_0'})
      self.assertDictEqual(large_scale_solve_0['step'], { 'python_class' : 'eumopps.catalogue.step.StepOnce'})
      self.assertEqual(len(large_scale_solve_0['newdatasets']), 1)
      self.assertEqual(large_scale_solve_0['newdatasets'][0]['name'], 'solution_large_scale_0')
      self.assertListEqual(large_scale_solve_0['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_large_scale_0/', 'solution_large_scale_0.pickle'])


      large_scale_input_2 = result['operations'][12]
      self.assertEqual(large_scale_input_2['name'], 'large_scale_input_2')
  
      self.assertEqual(large_scale_input_2['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'solution_climatology_2')
      self.assertEqual(large_scale_input_2['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_input_2['runmodule']['storage_local']['statefilename_read']['datasetname'], 'solution_local_1')
      self.assertDictEqual(large_scale_input_2['runmodule']['storage_local']['time_index'], {'python_class': 'eumopps.catalogue.placeholder.StepIndex'})
      self.assertEqual(large_scale_input_2['runmodule']['insitu_biases'], 1)
      self.assertEqual(large_scale_input_2['runmodule']['global_biases'], 1)
      self.assertEqual(large_scale_input_2['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(large_scale_input_2['runmodule']['method'], 'EXACT')
      self.assertEqual(large_scale_input_2['runmodule']['compute_sample'], 0)
      self.assertEqual(large_scale_input_2['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(large_scale_input_2['runmodule']['component_index'], 1)
      self.assertDictEqual(large_scale_input_2['runmodule']['storage_large_scale']['measurementfilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'measurement_large_scale_2'})
      self.assertDictEqual(large_scale_input_2['runmodule']['storage_large_scale']['measurement_time_index_write'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })
      self.assertEqual(large_scale_input_2['step']['start'], start_date+'000000')
      self.assertEqual(large_scale_input_2['step']['end'], end_date+'000000')
      self.assertEqual(len(large_scale_input_2['newdatasets']), 1)
      self.assertEqual(large_scale_input_2['newdatasets'][0]['name'], 'measurement_large_scale_2')
      self.assertListEqual(large_scale_input_2['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['measurement_large_scale_2/%Y/', 'measurement_large_scale_2_%Y%m%d.pickle'])

      large_scale_solve_2 = result['operations'][13]
      self.assertEqual(large_scale_solve_2['name'], 'large_scale_solve_2')
  
      self.assertEqual(large_scale_solve_2['runmodule']['storage_climatology']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_solve_2['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_solve_2['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(large_scale_solve_2['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.solve')
      self.assertEqual(large_scale_solve_2['runmodule']['insitu_biases'], 1)
      self.assertEqual(large_scale_solve_2['runmodule']['global_biases'], 1)
      self.assertEqual(large_scale_solve_2['runmodule']['compute_uncertainties'], 1)
      self.assertEqual(large_scale_solve_2['runmodule']['method'], 'APPROXIMATED')
      self.assertEqual(large_scale_solve_2['runmodule']['component_index'], 1)
      self.assertEqual(large_scale_solve_2['runmodule']['compute_sample'], 1)
      self.assertEqual(large_scale_solve_2['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertDictEqual(large_scale_solve_2['runmodule']['storage_large_scale']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'solution_large_scale_2'})
      self.assertDictEqual(large_scale_solve_2['runmodule']['storage_large_scale']['marginal_std_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'uncertainties_large_scale'})
      self.assertDictEqual(large_scale_solve_2['runmodule']['storage_large_scale']['sample_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'sample_large_scale'})
      self.assertEqual(large_scale_solve_2['step'], { 'python_class' : 'eumopps.catalogue.step.StepOnce'})
      self.assertEqual(large_scale_solve_2['newdatasets'][0]['name'], 'solution_large_scale_2')
      self.assertEqual(len(large_scale_solve_2['newdatasets']), 3)
      self.assertListEqual(large_scale_solve_2['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_large_scale_2/', 'solution_large_scale_2.pickle'])
      self.assertEqual(large_scale_solve_2['newdatasets'][1]['name'], 'uncertainties_large_scale')
      self.assertListEqual(large_scale_solve_2['newdatasets'][1]['subsets'][0]['layout']['patterns'], ['uncertainties_large_scale/', 'uncertainties_large_scale.pickle'])
      self.assertEqual(large_scale_solve_2['newdatasets'][2]['name'], 'sample_large_scale')
      self.assertListEqual(large_scale_solve_2['newdatasets'][2]['subsets'][0]['layout']['patterns'], ['sample_large_scale/', 'sample_large_scale.pickle'])


      # Local scale
      local_input_and_solve_2 = result['operations'][14]
      self.assertEqual(local_input_and_solve_2['name'], 'local_input_and_solve_2')
  
      self.assertEqual(local_input_and_solve_2['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'solution_climatology_2')
      self.assertEqual(local_input_and_solve_2['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], 'solution_large_scale_2')
      self.assertEqual(local_input_and_solve_2['runmodule']['storage_local']['statefilename_read']['datasetname'], None)
      self.assertEqual(local_input_and_solve_2['runmodule']['python_function'], 'eustace.analysis.advanced_standard.examples.example_eustace.process_inputs')
      self.assertEqual(local_input_and_solve_2['runmodule']['insitu_biases'], 1)
      self.assertEqual(local_input_and_solve_2['runmodule']['global_biases'], 1)
      self.assertEqual(local_input_and_solve_2['runmodule']['compute_uncertainties'], 1)
      self.assertEqual(local_input_and_solve_2['runmodule']['method'], 'APPROXIMATED')
      self.assertEqual(local_input_and_solve_2['runmodule']['compute_sample'], 1)
      self.assertEqual(local_input_and_solve_2['runmodule']['sample_size'], definitions.GLOBAL_SAMPLE_SHAPE[3])

      self.assertEqual(local_input_and_solve_2['runmodule']['component_index'], 2)
      self.assertDictEqual(local_input_and_solve_2['runmodule']['storage_local']['statefilename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'solution_local_2'})
      self.assertDictEqual(local_input_and_solve_2['runmodule']['storage_local']['marginal_std_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'uncertainties_local'})
      self.assertDictEqual(local_input_and_solve_2['runmodule']['storage_local']['sample_filename_write'], { 'python_class' : 'eumopps.catalogue.placeholder.OutputFile', 'datasetname' : 'sample_local'})
      self.assertDictEqual(local_input_and_solve_2['runmodule']['storage_local']['time_index'], { "python_class" : "eumopps.catalogue.placeholder.StepIndex" })
      self.assertEqual(local_input_and_solve_2['step']['start'], start_date+'000000')
      self.assertEqual(local_input_and_solve_2['step']['end'], end_date+'000000')
      self.assertEqual(len(local_input_and_solve_2['newdatasets']), 3)
      self.assertEqual(local_input_and_solve_2['newdatasets'][0]['name'], 'solution_local_2')
      self.assertListEqual(local_input_and_solve_2['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['solution_local_2/%Y/', 'solution_local_2_%Y%m%d.pickle'])
      self.assertEqual(local_input_and_solve_2['newdatasets'][1]['name'], 'uncertainties_local')
      self.assertListEqual(local_input_and_solve_2['newdatasets'][1]['subsets'][0]['layout']['patterns'], ['uncertainties_local/%Y/', 'uncertainties_local_%Y%m%d.pickle'])
      self.assertEqual(local_input_and_solve_2['newdatasets'][2]['name'], 'sample_local')
      self.assertListEqual(local_input_and_solve_2['newdatasets'][2]['subsets'][0]['layout']['patterns'], ['sample_local/%Y/', 'sample_local_%Y%m%d.pickle'])

      # Output gridding
      output_grid = result['operations'][15]
      self.assertEqual(output_grid['name'], 'output_grid')
  
      self.assertEqual(output_grid['runmodule']['storage_climatology']['statefilename_read']['datasetname'], 'solution_climatology_2')
      self.assertEqual(output_grid['runmodule']['storage_climatology']['marginal_std_filename_read']['datasetname'], 'uncertainties_climatology')
      self.assertEqual(output_grid['runmodule']['storage_climatology']['sample_filename_read']['datasetname'], 'sample_climatology')
      self.assertEqual(output_grid['runmodule']['storage_large_scale']['statefilename_read']['datasetname'], 'solution_large_scale_2')
      self.assertEqual(output_grid['runmodule']['storage_large_scale']['marginal_std_filename_read']['datasetname'], 'uncertainties_large_scale')
      self.assertEqual(output_grid['runmodule']['storage_large_scale']['sample_filename_read']['datasetname'], 'sample_large_scale')
      self.assertEqual(output_grid['runmodule']['storage_local']['statefilename_read']['datasetname'], 'solution_local_2')
      self.assertEqual(output_grid['runmodule']['storage_local']['marginal_std_filename_read']['datasetname'], 'uncertainties_local')
      self.assertEqual(output_grid['runmodule']['storage_local']['sample_filename_read']['datasetname'], 'sample_local')
      self.assertEqual(output_grid['runmodule']['outputfile']['datasetname'], 'eustace_infilled')
      self.assertEqual(output_grid['runmodule']['insitu_biases'], 1)
      self.assertEqual(output_grid['runmodule']['global_biases'], 1)
      self.assertEqual(output_grid['step']['start'], start_date+'000000')
      self.assertEqual(output_grid['step']['end'], end_date+'000000')
      self.assertEqual(output_grid['runmodule']['compute_uncertainties'], 0)
      self.assertEqual(output_grid['runmodule']['method'], 'EXACT')
      self.assertEqual(output_grid['runmodule']['compute_sample'], 0)

      self.assertEqual(output_grid['newdatasets'][0]['name'], 'eustace_infilled')
      self.assertListEqual(output_grid['newdatasets'][0]['subsets'][0]['layout']['patterns'], ['eustace_infilled/%Y/', 'eustace_infilled_%Y%m%d.nc'])

    def test_create_descriptor(self):
      current_dir=''
      template_operation = self.test_template_file
      template_grid = self.test_template_file
      start_date = '19870405'
      end_date = '19880405'
      bias = 1
      filler = FullstaceJsonFiller(3, 
                   current_dir, './',
                   template_operation, '', template_grid,
                   start_date, end_date,
                   1, 1,
                   'climatology_input', 'measurement_climatology', 'climatology_solve', 'solution_climatology', 
                   'large_scale_input', 'measurement_large_scale', 'large_scale_solve', 'solution_large_scale',
                   'local_input_and_solve', 'solution_local', 
                   'output_grid', 'eustace_infilled',
                   'output_climatology_name', 'output_large_scale_name', 'output_local_name')
                   
      filler.create_descriptor('baba.json')
      self.assertTrue(os.path.exists('./baba.json'))
      os.remove('./baba.json')
    
    def test_return_output_patterns(self):
      
      for input_name, output_list in zip(['A', 'abgs'], [['A/%Y/', 'A_%Y%m%d.nc'], ['abgs/%Y/', 'abgs_%Y%m%d.nc']]):
          self.assertListEqual(FullstaceJsonFiller.return_output_patterns(input_name, 'nc', 'measurement'), output_list)

      for input_name, output_list in zip(['A', 'abgs'], [['A/%Y/', 'A_%Y%m%d.nc'], ['abgs/%Y/', 'abgs_%Y%m%d.nc']]):
          self.assertListEqual(FullstaceJsonFiller.return_output_patterns(input_name, 'nc', 'space_solution'), output_list)

      for input_name, output_list in zip(['A', 'abgs'], [['A/', 'A.nc'], ['abgs/', 'abgs.nc']]):
          self.assertListEqual(FullstaceJsonFiller.return_output_patterns(input_name, 'nc', 'spacetime_solution'), output_list)

      for input_name, output_list in zip(['A', 'abgs'], [['A/', 'A.nc'], ['abgs/', 'abgs.nc']]):
          self.assertRaises(ValueError, FullstaceJsonFiller.return_output_patterns, input_name, 'nc', 'spjkbchw_time_solution')


    def test_check_none(self):
      
      for test_dictionary in self.dictionaries_with_None:
          self.assertRaises(ValueError, FullstaceJsonFiller.check_none, test_dictionary)

# To uncomment when more detailed tests on descriptor generation are required  
"""
    def test_against_original_descriptor(self):
      
      current_dir='test/'
      template_operation = 'test_advstd_operation_template.json'
      template_grid = 'test_advstd_outputgrid_template.json'
      start_date = '20060101'
      end_date = '20101231'
      bias = 1
      filler = FullstaceJsonFiller(1, 
                   current_dir, './',
                   template_operation, template_grid,
                   start_date, end_date,
                   1, 1,
                   'climatology_input', 'measurement_climatology', 'climatology_solve', 'solution_climatology', 
                   'large_scale_input', 'measurement_large_scale', 'large_scale_solve', 'solution_large_scale',
                   'local_input_and_solve', 'solution_local', 
                   'output_grid', 'eustace_example_infilled')
                   
      result = filler.fill_descriptor()['operations']
      with open(self.original_descriptor_file) as fp:
    comparison = json.load(fp)
      comparison = comparison['operations']
      self.assertEqual(len(result), len(comparison))
      
      self.check_all_keys(result, comparison)
      self.check_all_items(result, comparison)
      
    def check_all_keys(self, descriptor, comparison):

    if isinstance(descriptor, dict) and isinstance(comparison, dict):
        if len(descriptor) != len(comparison):
          message = 'Sub dictionaries differ!\nlen(descriptor) = {}, keys = {}\nlen(comparison)={}, keys={}'.format(len(descriptor), descriptor.keys(), len(comparison), comparison.keys())
          raise ValueError(message)
        else:
          left = descriptor.keys()
          left.sort()
          right = comparison.keys()
          right.sort()
          self.assertListEqual(left, right)
          for key in descriptor.keys():
        return self.check_all_keys(descriptor[key], comparison[key])
    elif isinstance(descriptor, list) and isinstance(comparison, list):
        return [ self.check_all_keys(item_1, item_2) for item_1, item_2 in zip(descriptor, comparison) ]
        
    def check_all_items(self, descriptor, comparison):

    if isinstance(descriptor, dict) and isinstance(comparison, dict):
       for key in descriptor.keys():
          return self.check_all_items(descriptor[key], comparison[key])
    elif isinstance(descriptor, list) and isinstance(comparison, list):
       return [ self.check_all_items(item_1, item_2) for item_1, item_2 in zip(descriptor, comparison) ]
    else:
       self.assertEqual(descriptor, comparison)
"""
