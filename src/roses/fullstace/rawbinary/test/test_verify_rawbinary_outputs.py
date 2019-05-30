""" Verify outputs from the insitu_land test run of the
Satstace raw binary generation code.
"""

from __future__ import print_function

__version__ = '0.1'
__author__ = 'J Winn'

import os
import fnmatch
import json
from pprint import pprint
import time

import eustaceconfig

EUSTPATH = eustaceconfig.SYSTEM_PATH
USER = os.environ.get('USER')
USERDIR = os.path.join('/work/scratch/', USER)
RAW = os.path.join(EUSTPATH, 'src/roses/fullstace/rawbinary')
RAWTEST = os.path.join(EUSTPATH, 'src/roses/fullstace/rawbinary/test')
QC_HO = os.path.join(USERDIR, 'rawbinary_test_land_qc_ho')
NOQC_HO = os.path.join(USERDIR, 'rawbinary_test_land_noqc_ho')


def print_inputs(jsonfile):
    """ Access and list the dataset name, path location, subsets and search patterns
    for json parameter input files
    """
    inputs = os.path.join(RAWTEST, jsonfile)
    with open(inputs) as filex:
        data = json.load(filex)
    inps = [0, 1, 2, 3, 4]
    for inp in inps:
        pprint(data['datasets'][inp]['name'])
        pprint(data['datasets'][inp]['path'])
        pprint(data['datasets'][inp]['subsets'][0]['layout']['patterns'])
        if inp != 1:
            pprint(data['datasets'][inp]['subsets'][1]['layout']['patterns'])
        print('---')

def print_generation(jsonfile):
    """ Access and list the parameters within the json param input files
    """
    inputs = os.path.join(RAWTEST, jsonfile)
    with open(inputs) as filex:
        data = json.load(filex)
    num_operations = len(data['operations'])
    inps = range(0, num_operations)
    for inp in inps:
        pprint(data['operations'][inp]['name']+ ' Years:  '
               +data['operations'][inp]['step']['start']+ '  to  '
               +data['operations'][inp]['step']['end'])
        print('---')

def get_files(directory, verbose):
    """ Access and lists the files within supplied directory in order to
    report the min and max and total count of years for each set of variables
    """
    dummy_path, dirs, dummy_files = next(os.walk(directory))
    all_subs = []
    for dirx in dirs:
        subs = os.path.join(directory, dirx)
        all_subs.append(subs)
        dummy_path, years, dummy_files = next(os.walk(subs))
        print('')
        print ('Data: '+str(dirx)+'  '+ str(len(years))+' years '
               +' from '+min(years)+' to '+max(years))
        print('')
        for year in years:
            sub = os.path.join(subs, year)
            dummy_subpath, dummy_subdir, subfiles = next(os.walk(sub))
            count = len(subfiles)
            #print(count)
            #print(year)
            test_min = len(fnmatch.filter(subfiles, '*Tmin_*'))
            test_max = len(fnmatch.filter(subfiles, '*Tmax_*'))
            test_mean = len(fnmatch.filter(subfiles, '*Tmean_*'))
            test_mobile = len(fnmatch.filter(subfiles, '*Tmean_mobile*'))
            if verbose:
                print(' Directory:  '+str(dirx)+'  Year:  '+str(year)
                      +' Files: '+str(count)+' Tmin: '+str(test_min)
                      +' Tmax '+str(test_max)+'  Tmean '+str(test_mean)
                      +' Tmean mobile '+str(test_mobile))
    return all_subs

# run checks --

# no tests were carried out for options without the holdout file due to difficulty
# in specifying a lack of holdout file in the json file

print('')
print('rawbinary input test files')
print_inputs('test_rawbinary_inputs.json')

print('---')
print('rawbinary generation test files')
print('-- no qc flags, with holdout file --')
print_generation('test_rawbinary_generation_land_noqc_with_holdout.json')
print('-- qc flags, with holdout file --')
print_generation('test_rawbinary_generation_land_with_qc_with_holdout.json')

print('---')
print(' checking test run file for run: land, with qc, with holdout')
print('Directory: %s ' %(QC_HO) )
get_files(QC_HO, True)
print('---')

print('---')
print(' checking test run file for run: land, with no qc, with holdout')
print('Directory: %s ' %(NOQC_HO) )
get_files(NOQC_HO, True)
print('---')



