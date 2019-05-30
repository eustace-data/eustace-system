import argparse
import lsf_submission_group

def eustace_early_grid(interation_number,
                      cataloguefilename = '/work/scratch/cmorice/advanced_standard/catalogue.nc',
                      wait_options = ""):
                          
    job_ids = lsf_submission_group.submit_eumopps_jobs(cataloguefilename, 'output_grid_preview_'+str(interation_number), wait_options=wait_options, lsf_script = lsf_submission_group.LSF_BATCH_SCRIPT_TALL, output_directory=lsf_submission_group.OUTDIR)
    wait_options = lsf_submission_group.done_conditon(job_ids)
    return wait_options

def main():

    print 'Submit early look gridding jobs prior to final analysis iteration. Output contains no uncertainty information.'
    
    parser = argparse.ArgumentParser(description='Generation of json descriptor for iterating operations')
    parser.add_argument('cataloguefilename', help='EUMOPPS catalogue file location')    
    parser.add_argument('--interation_number', type=int, default=0, help='interation_number index to run')
    parser.add_argument('--wait_options', default = '', help='prerequisite conditions for lsf jobs to start')

    args = parser.parse_args()

    eustace_early_grid(interation_number = args.interation_number,
                            cataloguefilename = args.cataloguefilename,
                            wait_options = args.wait_options)

if __name__ == '__main__':
    
    main()
