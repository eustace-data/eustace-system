"""Submit EUMOPPS catalogue tasks using LSF"""

import sys
import subprocess
import re
import os

from eumopps.catalogue.fileio.formatnetcdf import CatalogueReaderNetCDF

from StringIO import StringIO

LSF_BATCH_SCRIPT_GROUP = "lsf_commandrun_group.sh"

LSF_BATCH_SCRIPT = "lsf_commandrun.sh"
LSF_BATCH_SCRIPT_WIDE = "lsf_commandrun.sh"
LSF_BATCH_SCRIPT_TALL = "lsf_commandrun_tall.sh"
LSF_SOLVE_SCRIPT = "lsf_commandrun_solver_testing.sh" # "lsf_commandrun_solver.sh"

LSF_SOLVE_SCRIPT_LARGESCALE = "lsf_commandrun_solver_large_scale_high_memory.sh" # "lsf_commandrun_solver.sh"
LSF_SOLVE_SCRIPT_CLIMATOLOGY = "lsf_commandrun_solver_climatology.sh" # "lsf_commandrun_solver.sh"

LSF_BATCH_SCRIPT_GROUPED="lsf_commandrun_group.sh"

OUTDIR = "/work/scratch/cmorice/advanced_standard/"

def done_conditon(job_list):
    """Combine a list of lsf job name strings into a single string of lsf done conditions"""
    
    print job_list
    
    if len(job_list) == 0:
        wait_condition = ""
    else:
        wait_condition = " && ".join(["done("+job_id+")" for job_id in job_list])
    
    return wait_condition

def run_cmd( arguments, environmentvars, submission_script ):
    """Run a command using subprocess.Popen with specified environment stdin redirected from submission_script"""
    
    with open(submission_script) as myinput:
        
        process = subprocess.Popen( arguments,
                                    env=environmentvars,
                                    stdin=myinput,
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE)
           
        procout, procerr = process.communicate()
        
        procstatus = process.returncode   
    
    return procstatus, procout, procerr


def submit_eumopps_jobs(cataloguefilename, module_name, wait_options="", lsf_script = LSF_BATCH_SCRIPT, eumopps_batchsize = 366, output_directory=OUTDIR):
    """Submit EUMOPPS catalogue tasks using LSF"""
    
    # Get the job specification from the EUMOPPS catalogue
    eumopps_catalogue = os.path.abspath( cataloguefilename )
    
    commandcount = CatalogueReaderNetCDF().operationcount(eumopps_catalogue, module_name, 1)    
    batchcount = CatalogueReaderNetCDF().operationcount(eumopps_catalogue, module_name, eumopps_batchsize)
    if commandcount == 0:
        raise ValueError("Task "+module_name+" not found in catalogue")
    
    print( "Running: "+lsf_script+" for EUMOPPS module "+ module_name +" with "+str(batchcount)+" batches of "+str(eumopps_batchsize) )
    print( "Total tasks: "+str(commandcount) )

    # Get the full name of the scipt to be passed to bsub
    lsf_filename = os.path.abspath( lsf_script )

    # Setup and submit the jobs in batches collecting the job id for each batch
    submitted_tasks = 0
    job_ids = []
    for eumopps_batchnumber in range(batchcount):
        #if eumopps_batchnumber != 156:
            #continue
        #if not eumopps_batchnumber in [76, 118]:
            #continue
        
        if eumopps_batchnumber == batchcount - 1:
            tasks_this_submission = commandcount - submitted_tasks
        else:
            tasks_this_submission = eumopps_batchsize
        
        batchstart=submitted_tasks+1
        batchend=batchstart+tasks_this_submission-1
        
        # Setup environment variables used in the lsf_script
        job_environ = { "EUMOPPS_CATALOGUE":  eumopps_catalogue,
                        "EUMOPPS_MODULENAME": module_name,
                        "EUMOPPS_BATCHSIZE":  str(eumopps_batchsize),
                        "EUMOPPS_BATCHCOUNT": str(commandcount),
                        "EUMOPPS_BATCHNUMBER":str(eumopps_batchnumber),
                        "OUTDIR": output_directory
                       }
        
        submission_environ = dict(os.environ, **job_environ)
        
        # Construct the batch array including specification of batch size for lsf
        batcharray=module_name+"_batch"+str(eumopps_batchnumber)+"["+str(batchstart)+"-"+str(batchend)+"]"
        
        # Let the user know what is going on
        print( "Submitting:   "+lsf_script+" for EUMOPPS module "+ module_name +" for ["+str(batchstart)+"-"+str(batchend)+"] of "+str(commandcount)+" tasks" )
        print( "Dependencies: "+wait_options )
        print( "Batch name:   "+batcharray )
        
        # Construct the bsub command
        if wait_options == "":
            args = ["bsub", "-r", "-J", batcharray]
        else:
            args = ["bsub", "-r", "-w", wait_options, "-J", batcharray]

        # Submit the jobs
        procstatus, procout, procerr = run_cmd(args , submission_environ, lsf_filename) 
        
        if procstatus != 0:
            print procerr
            raise RuntimeError('Submission failed')
        
        print "Submission output:"
        print procout
        
        # Get the job id for the submitted batch
        job_id = re.compile('<([0-9]+)>').search( procout ).group(1)
        
        print( "Submitted as job id: "+job_id )
        
        # Collect job id
        job_ids.append(job_id)
        submitted_tasks += tasks_this_submission 
        
    return job_ids

    
def submit_eumopps_jobs_grouped(cataloguefilename, module_name, wait_options="", lsf_script = LSF_BATCH_SCRIPT, jobs_per_lsfbatch = 100, tasks_per_job  = 20, output_directory=OUTDIR):
    """Submit EUMOPPS catalogue tasks using LSF
    
    This version groups a large number of EUMOPPS tasks into invididual LSF submissions.
    Intended to reduce LSF scheduler problems associated having large numbers of dependencies
    on large numbers of jobs.
    
    Splits large numbers of jobs into groups of n_tasks_per_job with n_jobs_per_batch in each of n_batches.
    
    tasks_per_job  = 30   # number of eumopps operations in a lsf job
    jobs_per_lsfbatch = 100  # number of jobs in an lsf batch submission
    
    """
    
    eumopps_batchsize = tasks_per_job * jobs_per_lsfbatch    # A eumopps batch of operations is not the same as an lsf batch of jobs
                                                             # This is the number of eumopps operations in an lsf batch.
    
    # Get the job specification from the EUMOPPS catalogue
    eumopps_catalogue = os.path.abspath( cataloguefilename )
    
    commandcount = CatalogueReaderNetCDF().operationcount(eumopps_catalogue, module_name, 1)                # total number of operations to be run
    batchcount = CatalogueReaderNetCDF().operationcount(eumopps_catalogue, module_name, eumopps_batchsize)  # number of batches for commandcount jobs with for a specified batch size
    jobcount = CatalogueReaderNetCDF().operationcount(eumopps_catalogue, module_name, tasks_per_job)        # total number of lsf jobs to be submitted in all batches
    
    if commandcount == 0:
        raise ValueError("Task "+module_name+" not found in catalogue")
    
    print( "Running: "+lsf_script+" for EUMOPPS module "+ module_name +" with "+str(batchcount)+" batches of "+str(eumopps_batchsize) )
    print( "Total tasks: "+str(commandcount) )

    # Get the full name of the scipt to be passed to bsub
    lsf_filename = os.path.abspath( lsf_script )

    # Setup and submit the jobs in batches collecting the job id for each batch
    submitted_jobs = 0
    job_ids = []
    for eumopps_batchnumber in range(batchcount):
        
        
        
        if eumopps_batchnumber == batchcount - 1:
            jobs_this_submission = jobcount - submitted_jobs
        else:
            jobs_this_submission = jobs_per_lsfbatch
        
        startjob=submitted_jobs+1
        endjob=startjob+jobs_this_submission-1
        
        #if not eumopps_batchnumber in [3,]:
            #continue
        #if eumopps_batchnumber >= 21:
            #submitted_jobs += jobs_this_submission
            #continue
        
        # Setup environment variables used in the lsf_script
        job_environ = { "EUMOPPS_CATALOGUE":  eumopps_catalogue,
                        "EUMOPPS_MODULENAME": module_name,
                        "EUMOPPS_BATCHSIZE":  str(eumopps_batchsize),   # Total number of eumopps operations considered to be a batch in eumopps. *Not the number of LSF jobs*
                        "EUMOPPS_BATCHCOUNT": str(commandcount),        # Total number of operations in batch
                        "EUMOPPS_BATCHNUMBER":str(eumopps_batchnumber), # Use by EUMOPPS to identify which set of jobs are being run
                        "EUMOPPS_NJOBS": str(jobs_per_lsfbatch),     # 
                        "EUMOPPS_NTASKS":str(tasks_per_job),     # The number of tasks to be run consecutively within each lsf job
                        "OUTDIR": output_directory
                       }
        
        submission_environ = dict(os.environ, **job_environ)
        
        # Construct the batch array including specification of batch size for lsf
        batcharray=module_name+"_batch"+str(eumopps_batchnumber)+"["+str(startjob)+"-"+str(endjob)+"]"
        
        # Let the user know what is going on
        print( "Submitting:   "+lsf_script+" for EUMOPPS module "+ module_name +" for ["+str(startjob)+"-"+str(endjob)+"] of "+str(jobcount)+" jobs" )
        print( "Dependencies: "+wait_options )
        print( "Batch name:   "+batcharray )
        
        # Construct the bsub command
        if wait_options == "":
            args = ["bsub", "-r", "-J", batcharray]
        else:
            args = ["bsub", "-r", "-w", wait_options, "-J", batcharray]

        dryrun = False #eumopps_batchnumber < (batchcount - 1)
        if not dryrun:
            # Submit the jobs
            procstatus, procout, procerr = run_cmd(args , submission_environ, lsf_filename) 
            
            if procstatus != 0:
                print procerr
                raise RuntimeError('Submission failed')
            
            print "Submission output:"
            print procout
            
            # Get the job id for the submitted batch
            job_id = re.compile('<([0-9]+)>').search( procout ).group(1)
            
            print( "Submitted as job id: "+job_id )
            
            # Collect job id
            job_ids.append(job_id)
        
        submitted_jobs += jobs_this_submission
        #stop
        
    return job_ids

def eustace_single_iteration():
    """Runs the single iteration test run of the advanced standard analysis"""
    
    cataloguefilename = '/work/scratch/cmorice/advanced_standard/catalogue.nc'

    wait_options = ""
    job_ids = submit_eumopps_jobs(cataloguefilename, 'climatology_input', lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
    
    wait_options = done_conditon(job_ids)
    job_ids = submit_eumopps_jobs(cataloguefilename, 'climatology_solve', wait_options=wait_options, lsf_script = LSF_SOLVE_SCRIPT, output_directory=OUTDIR)
    
    wait_options = done_conditon(job_ids)
    job_ids = submit_eumopps_jobs(cataloguefilename, 'large_scale_input', wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
    
    wait_options = done_conditon(job_ids)
    job_ids = submit_eumopps_jobs(cataloguefilename, 'large_scale_solve', wait_options=wait_options, lsf_script = LSF_SOLVE_SCRIPT, output_directory=OUTDIR)
    
    wait_options = done_conditon(job_ids)
    job_ids = submit_eumopps_jobs(cataloguefilename, 'local_input_and_solve', wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
    
    wait_options = done_conditon(job_ids)
    job_ids = submit_eumopps_jobs(cataloguefilename, 'output_grid', wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)#, eumopps_batchsize = 10)

def eustace_multi_iteration(first_iteration = 0,
                            last_iteration = 2,
                            final_iteration = 2,
                            climatology_input = True,
                            climatology_solve = True,
                            large_scale_input = True,
                            large_scale_solve = True,
                            local_input_and_solve = True,
                            output_grid = True,
                            cataloguefilename = '/work/scratch/cmorice/advanced_standard/catalogue.nc',
                            wait_options = ""):
    
    for iteration_number in range(first_iteration, last_iteration+1):
        
        grid_output = False
        if iteration_number == final_iteration:
            if output_grid:
                grid_output = True
            run_prior_sample = True
        else:
            run_prior_sample = False

        #run_prior_sample = True

        print "iteration starts after :"+ wait_options
        wait_options = eustace_iteration(iteration_number,
                                         climatology_input = climatology_input,
                                         climatology_solve = climatology_solve,
                                         large_scale_input = large_scale_input,
                                         large_scale_solve = large_scale_solve,
                                         local_input_and_solve = local_input_and_solve,
                                         grid_output = grid_output,
                                         run_prior_sample = run_prior_sample,
                                         wait_options = wait_options,
                                         cataloguefilename = cataloguefilename)

def eustace_iteration(interation_number,
                      cataloguefilename = '/work/scratch/cmorice/advanced_standard/catalogue.nc',
                      wait_options = "",
                      climatology_input = True,
                      climatology_solve = True,
                      large_scale_input = True,
                      large_scale_solve = True,
                      local_input_and_solve = True,
                      grid_output = False,
                      run_prior_sample = False,
                      local_hyperparameters = False):
                      #grid_components = False):
    """Runs a specified iteration of the advanced standard analysis"""
    
    

    # climatology processing
    if climatology_input:
        job_ids = submit_eumopps_jobs(cataloguefilename, 'climatology_input_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT_TALL, output_directory=OUTDIR)
        wait_options = done_conditon(job_ids)
    if climatology_solve:
        job_ids = submit_eumopps_jobs(cataloguefilename, 'climatology_solve_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_SOLVE_SCRIPT_CLIMATOLOGY, output_directory=OUTDIR)
        wait_options = done_conditon(job_ids)
        if run_prior_sample:
            job_ids = submit_eumopps_jobs(cataloguefilename, 'climatology_solve_prior_sample_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_SOLVE_SCRIPT_CLIMATOLOGY, output_directory=OUTDIR)
            wait_options = done_conditon(job_ids)
    
    # largescale processing
    if large_scale_input:
        job_ids = submit_eumopps_jobs(cataloguefilename, 'large_scale_input_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT_TALL, output_directory=OUTDIR)
        wait_options = done_conditon(job_ids)
    if large_scale_solve:
        job_ids = submit_eumopps_jobs(cataloguefilename, 'large_scale_solve_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_SOLVE_SCRIPT_LARGESCALE, output_directory=OUTDIR)
        wait_options = done_conditon(job_ids)
        if run_prior_sample:
            job_ids = submit_eumopps_jobs(cataloguefilename, 'large_scale_solve_prior_sample_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_SOLVE_SCRIPT_LARGESCALE, output_directory=OUTDIR)
            wait_options = done_conditon(job_ids)
    
    # local processing
    if local_input_and_solve:
        if not local_hyperparameters:
            job_ids = submit_eumopps_jobs_grouped(cataloguefilename, 'local_input_and_solve_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT_GROUPED, output_directory=OUTDIR)
        else:
            job_ids = submit_eumopps_jobs_grouped(cataloguefilename, 'nonstationary_local_input_and_solve_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT_GROUPED, output_directory=OUTDIR)
        wait_options = done_conditon(job_ids)
    
    # output grid
    if grid_output:
        job_ids = submit_eumopps_jobs(cataloguefilename, 'output_grid_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT_TALL, output_directory=OUTDIR)
        wait_options = done_conditon(job_ids)
        
        ## extra early period output local solve and grid
        #job_ids = submit_eumopps_jobs(cataloguefilename, 'local_input_and_solve_early_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT_WIDE, output_directory=OUTDIR)
        #wait_options = done_conditon(job_ids)
        #job_ids = submit_eumopps_jobs(cataloguefilename, 'output_grid_early', wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT_TALL, output_directory=OUTDIR)
        #wait_options = done_conditon(job_ids)
        
    #if grid_components:
        #job_ids1 = submit_eumopps_jobs(cataloguefilename, 'climatology_grid_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
        #job_ids2 = submit_eumopps_jobs(cataloguefilename, 'large_scale_grid_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
        #job_ids3 = submit_eumopps_jobs(cataloguefilename, 'local_grid_'+str(interation_number), wait_options=wait_options, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
        #wait_options = done_conditon(job_ids1+job_ids2+job_ids3)
        
    return wait_options

def eustace_optimisation(iteration_number=None,
                         cataloguefilename = '/work/scratch/cmorice/advanced_standard/catalogue.nc',
                         wait_options = ""):

    if iteration_number is None:
        iter_string = ''
    else:
        iter_string = '_'+str(iteration_number)

    # split the local component into spatial model and bias components so that the bias sub-component can be used later
    job_ids = submit_eumopps_jobs(cataloguefilename, 'split_local_component'+iter_string, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
    wait_options = done_conditon(job_ids)
    
    # generate json files listing inputloader information for each time step for use outside of EUMOPPS
    job_ids = submit_eumopps_jobs(cataloguefilename, 'input_summaries_create'+iter_string, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
    wait_options = done_conditon(job_ids)
    
    # consolidate the inputloader information in json files into a single master json file
    job_ids = submit_eumopps_jobs(cataloguefilename, 'input_summaries_merge'+iter_string, lsf_script = LSF_BATCH_SCRIPT, output_directory=OUTDIR)
    wait_options = done_conditon(job_ids)
    
    # submit the optimisation tasks outside of EUMOPPS
    
    return wait_options

import argparse

def main():

    print 'Submission of advanced standard analysis jobs'
    
    parser = argparse.ArgumentParser(description='Generation of json descriptor for iterating operations')
    parser.add_argument('cataloguefilename', default='/work/scratch/cmorice/advanced_standard/catalogue.nc', help='catalogue file location')
    
    parser.add_argument('--first_iteration', type=int, default=0, help='first iteration index to run')
    parser.add_argument('--last_iteration',  type=int, default=2, help='last iteration index to run')
    parser.add_argument('--final_iteration', type=int, default=2, help='final iteration of the analysis at which gridding happens')
    
    parser.add_argument('--climatology_input',     type=int, default=1, help='run measurement input processing for climatology')
    parser.add_argument('--climatology_solve',     type=int, default=1, help='run solver for climatology')
    parser.add_argument('--large_scale_input',     type=int, default=1, help='run measurement input processing for climatology')
    parser.add_argument('--large_scale_solve',     type=int, default=1, help='run solver for climatology')
    parser.add_argument('--local_input_and_solve', type=int, default=1, help='run measurement input and solve for daily local')
    
    parser.add_argument('--output_grid',     type=int, default=1, help='run gridder on final iteration')
    
    parser.add_argument('--wait_options', default = '', help='prerequisite conditions for lsf jobs to start')

    args = parser.parse_args()

    eustace_multi_iteration(first_iteration = args.first_iteration,
                            last_iteration  = args.last_iteration,
                            final_iteration = args.final_iteration,
                            climatology_input = args.climatology_input,
                            climatology_solve = args.climatology_solve,
                            large_scale_input = args.large_scale_input,
                            large_scale_solve = args.large_scale_solve,
                            local_input_and_solve = args.local_input_and_solve,
                            output_grid = args.output_grid,
                            cataloguefilename = args.cataloguefilename,
                            wait_options = args.wait_options)

if __name__ == '__main__':
    
    main()
