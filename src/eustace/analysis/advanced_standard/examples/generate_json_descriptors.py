"""Generate json descriptor describing operations that have to be iterated."""

import argparse
import json
import numpy
import os.path

from eustaceconfig import WORKSPACE_PATH
from jsonfiller import FullstaceJsonFiller
      
def main():
    
    print 'Generation of json descriptor for iterating operations'
    
    parser = argparse.ArgumentParser(description='Generation of json descriptor for iterating operations')
    parser.add_argument('inpath', help='directory containing input template descriptors.\n Two templates are expected: "advstd_operation_template.json" and "advstd_outputgrid_template.json"')
    parser.add_argument('outpath', help='directory where the output should be redirected')
    
    parser.add_argument('first_iteration', type=int, default=0, help='first iteration to add to the descriptor')
    parser.add_argument('last_iteration', type=int, default=2, help='last iteration to add to the descriptor')
    parser.add_argument('grid_last_iteration', type=int, default=1, help='should samples be drawn and output be gridded on last iteration?')
    
    parser.add_argument('start_date', help='starting date (YYYYMMDD)')
    parser.add_argument('end_date', help='final date (YYYYMMDD)')
    parser.add_argument('descriptor_name', help='name of the file descriptor')

    parser.add_argument('--operation_template', default = 'advstd_operation_template.json', help='template operation descriptor name')
    parser.add_argument('--batch_operation_template', default = 'advstd_batch_operation_template.json', help='template operation descriptor name for processing time blocks')
    parser.add_argument('--output_grid_template', default = 'advstd_outputgrid_template.json', help='template gridding descriptor name')

    parser.add_argument('--land_biases', default=1, help='include insitu land homogenization bias terms')
    parser.add_argument('--global_biases', default=1, help='include global satellite bias terms')

    parser.add_argument('--measurement_climatology_command', default = 'climatology_input', help='template name for climatology measuremement command')
    parser.add_argument('--measurement_climatology_name', default = 'measurement_climatology', help='template name for climatology measuremement dataset')
    parser.add_argument('--solution_climatology_command', default = 'climatology_solve', help='template name for climatology solution command')
    parser.add_argument('--solution_climatology_name', default = 'solution_climatology', help='template name for climatology solution dataset')

    parser.add_argument('--measurement_large_scale_command', default = 'large_scale_input', help='template name for large scale measuremement command')
    parser.add_argument('--measurement_large_scale_name', default = 'measurement_large_scale', help='template name for large scale measuremement dataset')
    parser.add_argument('--solution_large_scale_command', default = 'large_scale_solve', help='template name for large scale solution command')
    parser.add_argument('--solution_large_scale_name', default = 'solution_large_scale', help='template name for large scale solution dataset')

    parser.add_argument('--local_scale_command', default = 'local_input_and_solve', help='template name for local scale measuremement and solution command')
    parser.add_argument('--solution_local_scale_name', default = 'solution_local', help='template name for local scale solution dataset')

    parser.add_argument('--output_grid_command', default = 'output_grid', help='template name for output grid command')
    parser.add_argument('--output_grid_name', default = 'eustace_infilled', help='template name for output grid dataset')
    
    parser.add_argument('--output_climatology_name', default = 'eustace_infilled', help='template name for output grid dataset')
    parser.add_argument('--output_large_scale_name', default = 'eustace_infilled', help='template name for output grid dataset')
    parser.add_argument('--output_local_name', default = 'eustace_infilled', help='template name for output grid dataset')

    args = parser.parse_args()
  	
    filler = FullstaceJsonFiller(args.first_iteration, args.last_iteration, args.grid_last_iteration,
				 args.inpath, args.outpath,
                                 args.operation_template, args.batch_operation_template, args.output_grid_template,
		                 args.start_date, args.end_date,
                                 args.land_biases, args.global_biases, 
                                 args.measurement_climatology_command, args.measurement_climatology_name, args.solution_climatology_command, args.solution_climatology_name,
                                 args.measurement_large_scale_command, args.measurement_large_scale_name, args.solution_large_scale_command, args.solution_large_scale_name,
                                 args.local_scale_command, args.solution_local_scale_name,
                                 args.output_grid_command, args.output_grid_name,
                                 args.output_climatology_name, args.output_large_scale_name, args.output_local_name)
                                 
    filler.create_batch_descriptor(args.descriptor_name, check_none=False)

    
    
if __name__=='__main__':

    # Call main method
    main()
