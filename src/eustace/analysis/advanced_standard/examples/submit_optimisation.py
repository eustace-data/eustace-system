
import os
from eustace.outputformats.ensuredirectory import ensuredirectory

BASE_DIVISIONS = 20


def n_operations(subdivision_level):
    """Total operations defined by triangulations."""

    number_of_operations = BASE_DIVISIONS
    
    for level in range(subdivision_level):
        number_of_operations = number_of_operations * 4
        
    return number_of_operations

def submit_jobs():

    #eustace, par-single, high-mem

    level = 7
    neighbourhood_level = 3
    regionspec = 'LocalSubRegion'
    #number_of_regions = 1
    number_of_regions = n_operations(neighbourhood_level)
    #number_of_regions = 20
    
    outdir = '/work/scratch/cmorice/advanced_standard/optimise/'
    ensuredirectory(outdir)
    
    for region_index in range(number_of_regions):
        logfile = os.path.join(outdir, 'optimize.'+str(level)+'.'+str(neighbourhood_level)+'.'+str(region_index)+'.out')
        errfile = os.path.join(outdir, 'optimize.'+str(level)+'.'+str(neighbourhood_level)+'.'+str(region_index)+'.err')        
        jobfile = os.path.join(outdir, 'optimize.'+str(level)+'.'+str(neighbourhood_level)+'.'+str(region_index)+'.sh')
        
        pycommand = ' '.join([  'python2.7',
                                #'-m eustace.analysis.advanced_standard.examples.example_eustace_few_days_optimization',
                                '-m eustace.analysis.advanced_standard.examples.optimise',
                                '--neighbourhood_level', str(neighbourhood_level),
                                '--region_index', str(region_index),
                                '--regionspec', regionspec])

        job_spec = ''.join([ '#!/bin/bash\n',
                                '#BSUB -c 2:00\n',
                                '#BSUB -W 2:00\n',
                                #'#BSUB -R "rusage[mem=10]"',
                                '#BSUB -n 1\n',
                                '#BSUB -q short-serial\n',
                                '#BSUB -oo', logfile+'\n',
                                '#BSUB -eo', errfile+'\n',
                                'module load intel/cce/15.0.090\n',
                                'module load intel/mkl/11.3.1.150\n',
                                pycommand+'\n',])

        with open(jobfile, 'w') as f:
            f.write(job_spec)
        
        run_string = 'bsub < '+jobfile
        
        print "runing:"
        print run_string
        
        os.system(run_string)

if __name__ == '__main__':

    submit_jobs()