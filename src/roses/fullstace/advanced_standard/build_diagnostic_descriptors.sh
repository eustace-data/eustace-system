#!/bin/bash
#
# Build a json descriptor for a set of analysis iterations
#

# Expect output directory as parameter
if [ ! -d "$1" ] ; then
  echo "ERROR: expected output directory as parameter"
  echo "Usage: $0 <outputdirectory>"
  exit 1
fi 

# find full path of templates directory
# assuming it's in the same directory as this script
SCRIPT_PATH=`dirname ${BASH_SOURCE[0]}`
INPATH=${SCRIPT_PATH}'/templates'
OUTPATH="$1"
N_ITERATIONS="5"
START_DATE="18800101"
END_DATE="20151231"
DESCRIPTOR_NAME="advstd_diagnostics.json"

OPERATION_TEMPLATE="advstd_operation_template.json"
BATCH_OPERATION_TEMPLATE="advstd_batch_operation_template.json"
OUTPUT_GRID_TEMPLATE="advstd_outputgrid_template.json"
LAND_BIASES="1"
GLOBAL_BIASES="1"
MEASUREMENT_CLIMATOLOGY_COMMAND="climatology_input" #"eustace.analysis.advanced_standard.examples.moving_climatology.process_inputs"
MEASUREMENT_CLIMATOLOGY_NAME="climatology_increment" #"climatology_input"
SOLUTION_CLIMATOLOGY_COMMAND="climatology_solve" #"eustace.analysis.advanced_standard.examples.moving_climatology.solve"
SOLUTION_CLIMATOLOGY_NAME="climatology_solution"
MEASUREMENT_LARGE_SCALE_COMMAND="large_scale_input" #"eustace.analysis.advanced_standard.examples.moving_climatology.process_inputs"
MEASUREMENT_LARGE_SCALE_NAME="large_scale_increment"
SOLUTION_LARGE_SCALE_COMMAND="large_scale_solve" #"eustace.analysis.advanced_standard.examples.moving_climatology.solve"
SOLUTION_LARGE_SCALE_NAME="large_scale_solution"
LOCAL_SCALE_COMMAND="local_input_and_solve" #"eustace.analysis.advanced_standard.examples.moving_climatology.process_inputs"
SOLUTION_LOCAL_SCALE_NAME="local_solution"
OUTPUT_GRID_COMMAND="output_grid"
OUTPUT_GRID_NAME="eustace_analysis"
OUTPUT_CLIMATOLOGY_NAME="eustace_climatology_component"
OUTPUT_LARGE_SCALE_NAME="eustace_large_scale_component"
OUTPUT_LOCAL_NAME="eustace_local_component"

python -m eustace.analysis.advanced_standard.examples.generate_diagnostic_descriptors $INPATH $OUTPATH $N_ITERATIONS $START_DATE $END_DATE $DESCRIPTOR_NAME \
--operation_template $OPERATION_TEMPLATE \
--batch_operation_template $BATCH_OPERATION_TEMPLATE \
--output_grid_template $OUTPUT_GRID_TEMPLATE \
--land_biases $LAND_BIASES \
--global_biases $GLOBAL_BIASES \
--measurement_climatology_command $MEASUREMENT_CLIMATOLOGY_COMMAND \
--measurement_climatology_name $MEASUREMENT_CLIMATOLOGY_NAME \
--solution_climatology_command $SOLUTION_CLIMATOLOGY_COMMAND \
--solution_climatology_name $SOLUTION_CLIMATOLOGY_NAME \
--measurement_large_scale_command $MEASUREMENT_LARGE_SCALE_COMMAND \
--measurement_large_scale_name $MEASUREMENT_LARGE_SCALE_NAME \
--solution_large_scale_command $SOLUTION_LARGE_SCALE_COMMAND \
--solution_large_scale_name $SOLUTION_LARGE_SCALE_NAME \
--local_scale_command $LOCAL_SCALE_COMMAND \
--solution_local_scale_name $SOLUTION_LOCAL_SCALE_NAME \
--output_grid_command $OUTPUT_GRID_COMMAND \
--output_grid_name $OUTPUT_GRID_NAME \
--output_climatology_name $OUTPUT_CLIMATOLOGY_NAME \
--output_large_scale_name $OUTPUT_LARGE_SCALE_NAME \
--output_local_name $OUTPUT_LOCAL_NAME