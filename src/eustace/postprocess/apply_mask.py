
import argparse
import os.path
import copy
import numpy
from dateutil.relativedelta import relativedelta
import datetime
from dateutil import parser

import iris

from eustace.outputformats.outputvariable import OutputVariable
from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eumopps.version.svn import get_revision_id_for_module
from netCDF4 import default_fillvals
from eustace.outputformats import definitions
import eustace.timeutils.epoch

from operation_count import operation_dates
import flags

import copy

# bit mask flags used in output file
FLAG_TYPE = flags.FLAG_TYPE
TYPE_NAME = flags.TYPE_NAME
FLAG_MAX_N = flags.FLAG_MAX_N
FLAG_MAX_USED = flags.FLAG_MAX_USED
NULL_FLAG = flags.NULL_FLAG

# time window checks
DAY_FLAG = flags.DAY_FLAG
PRIOR_WINDOW_FLAG = flags.PRIOR_WINDOW_FLAG
POST_WINDOW_FLAG = flags.POST_WINDOW_FLAG

# calendar day checks
PRIOR_CALENDAR_FLAG  = flags.PRIOR_CALENDAR_FLAG
POST_CALENDAR_FLAG   = flags.POST_CALENDAR_FLAG

# slow component uncertatinty thresholds
CLIMATOLOGY_UNC_FLAG  = flags.CLIMATOLOGY_UNC_FLAG
LARGE_SCALE_UNC_FLAG  = flags.LARGE_SCALE_UNC_FLAG

# extreme value checks
AREAL_LOW_FLAG      = flags.AREAL_LOW_FLAG
AREAL_HIGH_FLAG     = flags.AREAL_HIGH_FLAG
EXTREME_LOW_FLAG    = flags.EXTREME_LOW_FLAG
EXTREME_HIGH_FLAG   = flags.EXTREME_HIGH_FLAG

# omitted data source flags
MISSING_MARINE_FLAG  = flags.MISSING_MARINE_FLAG

LOCATION_LOW_FLAG = flags.LOCATION_LOW_FLAG
LOCATION_HIGH_FLAG = flags.LOCATION_HIGH_FLAG

# missing data indicator
MISSING_FLAG_FLAG    = flags.MISSING_FLAG_FLAG

FLAG_MEANINGS = flags.FLAG_MEANINGS

# useful derived combinations
CALENDAR_DAY_FLAG    = flags.CALENDAR_DAY_FLAG # helper combination that say there is no observation nearby in time


# Fist items are checked in sequence as a logical_or. Fields are mask where any of the following conditions are true.
# For each item, the first index is a set of flags that together result in masking. The second index is a list of 
# flags that aways override and negate the mask, if they are not set, i.e. don't mask if the local constraint flag
# is not set.

FLAGS_TO_MASK = [ [MISSING_FLAG_FLAG,                       MISSING_MARINE_FLAG],   # Mask if calendar day has nearby constraint less than three times in full period.  Ignore if marine.
                  [LARGE_SCALE_UNC_FLAG,                    NULL_FLAG],             # Mask if large scale is not constrained. Removed day constraint because subject to noise.
                  [CLIMATOLOGY_UNC_FLAG,                    NULL_FLAG],             # Mask if climatology is not constrained. Removed day constraint because subject to noise.
                  [AREAL_LOW_FLAG,                          NULL_FLAG],             # Mask if failed low bound areal QC. Never ignore.
                  [AREAL_HIGH_FLAG,                         NULL_FLAG],             # Mask if failed high bound areal QC. Never ignore.
                  [EXTREME_LOW_FLAG,                        NULL_FLAG],             # Mask if failed extreme low bound threshold check. Never ignore.
                  [EXTREME_HIGH_FLAG,                       NULL_FLAG],             # Mask if failed extreme high bound threshold check. Never ignore.
                  [LOCATION_LOW_FLAG,                       NULL_FLAG],             # Extremely cold in comparison to location based bounds
                  [LOCATION_HIGH_FLAG,                      NULL_FLAG],             # Extremely hot in comparison to location based bounds
                 ]

print FLAGS_TO_MASK


def reset_marine_flag(flag_values):
    print flag_values.shape
    # clear the marine flag
    flag_values = numpy.bitwise_and( flag_values, ~MISSING_MARINE_FLAG )
    
    # set new marine flag
    LAND_VALUE = 0.0
    SEA_VALUE  = 100.0
    MARINE_CLASSIFICATION_THRESHOLD = 95.0
    COAST_FILE='/gws/nopw/j04/eustace/data/internal/climatology_covariates/coastal_influence.test.0.25_0.25.nc'
    COAST_FILE_VAR = 'coastal_influence'
    variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == COAST_FILE_VAR))
    coastal_influence = iris.load_cube( COAST_FILE, variable_constraint ).data[:,:]

    coastal_influence[0,:] = LAND_VALUE
    coastal_influence[-1,:] = SEA_VALUE
    coastal_influence[:,0] = coastal_influence[:,1] 
    coastal_influence[:,-1] = coastal_influence[:,-2]
    
    coastal_influence = coastal_influence.reshape(flag_values.shape)
    
    marine_locations = coastal_influence > MARINE_CLASSIFICATION_THRESHOLD

    flag_values[marine_locations] = flag_values[marine_locations] | MISSING_MARINE_FLAG
    print flag_values.shape
    return flag_values

"""

File I/O

"""
    
def save_flag_file(flag_values, processdate, outputfile):
    """Save the windowing flag information"""
    
    print "Saving: ", outputfile
    
    # Create directory structure if required
    if not os.path.exists(os.path.dirname(outputfile)):
        os.makedirs(os.path.dirname(outputfile))
    
    # Make the QC flag file
    filebuilder = FileBuilderGlobalField(
        outputfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis Flags',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional output',
        __name__, 
        '')
    
    # flag flag definition
    tas_qc_definition = OutputVariable( name='tas_qc',
                                dtype=FLAG_TYPE,
                                fill_value=default_fillvals[TYPE_NAME],
                                standard_name='air_temperature status_flag',
                                long_name='Quality control flags',
                                #valid_range=numpy.array([FLAG_TYPE(0), FLAG_TYPE(1)<<FLAG_MAX_N], FLAG_TYPE),
                                flag_masks=numpy.array([ FLAG_TYPE(1)<<n for n in range(FLAG_MAX_USED)], FLAG_TYPE),
                                flag_meanings=' '.join(FLAG_MEANINGS))      


    # Combine flag information and populate the output file 
    filebuilder.add_global_field( tas_qc_definition, flag_values.reshape(definitions.GLOBAL_FIELD_SHAPE) )
    
    # Close the netCDF4 dataset
    filebuilder.save_and_close()

def load_analysis(directory, iteration, datestring):
    """Load EUSTACE advanced standard analysis data"""
    
    analysis_cubes = iris.load(directory+'/eustace_analysis_'+str(iteration)+'/'+datestring[:4]+'/eustace_analysis_'+str(iteration)+'_'+datestring+'.nc')

    cube_names = [analysis_cubes[i].name() for i in range(len(analysis_cubes))]
    
    ensemble_members = []
    for i in range(len(analysis_cubes)):
        
        if cube_names[i] == u'air_temperature':
            analysis = analysis_cubes[i]
        elif cube_names[i] == u'Total uncertainty in average daily surface air temperature':
            uncertainty = analysis_cubes[i]
        elif cube_names[i] == u'Observation influence for average daily surface air temperature':
            influence = analysis_cubes[i]
        elif cube_names[i] == u'Average daily surface air temperature ensemble':
            ensemble_members.append(analysis_cubes[i])
        else:
            continue
        
    return analysis, uncertainty, influence, ensemble_members

def load_climatology(directory, iteration, datestring):
    """Load EUSTACE advanced standard analysis climatology component data"""
        
    analysis_cubes = iris.load(directory+'/eustace_climatology_component_'+str(iteration)+'/'+datestring[:4]+'/eustace_climatology_component_'+str(iteration)+'_'+datestring+'.nc')

    cube_names = [analysis_cubes[i].name() for i in range(len(analysis_cubes))]
    
    for i in range(len(analysis_cubes)):
        
        if cube_names[i] == u'air_temperature':
            analysis = analysis_cubes[i]
        elif cube_names[i] == u'Total uncertainty in average daily surface air temperature':
            uncertainty = analysis_cubes[i]
        elif cube_names[i] == u'Observation influence for average daily surface air temperature':
            influence = analysis_cubes[i]
        else:
            continue
        
    return analysis, uncertainty, influence

def load_large_scale(directory, iteration, datestring):
    """Load EUSTACE advanced standard analysis large scale component data"""
        
    analysis_cubes = iris.load(directory+'/eustace_large_scale_component_'+str(iteration)+'/'+datestring[:4]+'/eustace_large_scale_component_'+str(iteration)+'_'+datestring+'.nc')

    cube_names = [analysis_cubes[i].name() for i in range(len(analysis_cubes))]
    
    for i in range(len(analysis_cubes)):
        
        if cube_names[i] == u'air_temperature':
            analysis = analysis_cubes[i]
        elif cube_names[i] == u'Total uncertainty in average daily surface air temperature':
            uncertainty = analysis_cubes[i]
        elif cube_names[i] == u'Observation influence for average daily surface air temperature':
            influence = analysis_cubes[i]
        else:
            continue
        
    return analysis, uncertainty, influence

def load_local(directory, iteration, datestring):
    """Load EUSTACE advanced standard analysis local component data"""
    
    analysis_cubes = iris.load(directory+'/eustace_local_component_'+str(iteration)+'/'+datestring[:4]+'/eustace_local_component_'+str(iteration)+'_'+datestring+'.nc')

    cube_names = [analysis_cubes[i].name() for i in range(len(analysis_cubes))]
    
    for i in range(len(analysis_cubes)):
        
        if cube_names[i] == u'air_temperature':
            analysis = analysis_cubes[i]
        elif cube_names[i] == u'Total uncertainty in average daily surface air temperature':
            uncertainty = analysis_cubes[i]
        elif cube_names[i] == u'Observation influence for average daily surface air temperature':
            influence = analysis_cubes[i]
        else:
            continue
        
    return analysis, uncertainty, influence

def load_flags(directory, iteration, datestring):
    """Load the a mask defined by the list of flags in flags_to_mask"""
    
    analysis_cubes = iris.load(directory+'/'+datestring[:4]+'/eustace_analysis_'+str(iteration)+'_qc_flags_'+datestring+'.nc')

    cube_names = [analysis_cubes[i].name() for i in range(len(analysis_cubes))]
    print cube_names
    for i in range(len(analysis_cubes)):
        
        if cube_names[i] == u'Quality control flags':
            qc_flags = analysis_cubes[i]
        else:
            continue
    
    # Apply each of the masking conditions        
    #mask = numpy.logical_or.reduce([numpy.bitwise_and(qc_flags.data[:,:,:], flagspec[0] | flagspec[1]) == flagspec[0] for flagspec in flags_to_mask])
    #mask = numpy.logical_or.reduce([numpy.bitwise_and(qc_flags.data[:,:,:], flagspec[0]) == flagspec[0] - numpy.bitwise_and(qc_flags.data[:,:,:], flagspec[1]) == flagspec[1] for flagspec in flags_to_mask])
    
    #mask = compute_mask(qc_flags.data[:,:,:], flags_to_mask)
    
    return qc_flags

"""

Apply mask

"""

    
def compute_mask(qc_flags, flags_to_mask):
    """Derive a boolean mask from a set of qc_flags for a set of masking conditions specified in flags_to_mask"""
    
    mask = numpy.logical_or.reduce([numpy.bitwise_and(qc_flags, flagspec[0] | flagspec[1]) == flagspec[0] for flagspec in flags_to_mask])

    return mask

def apply_mask(analysis_directory, iteration, mask_directory, mask_directory2, flags_to_mask, processdate, output_directory):
    """Apply QC mask to analysis"""

    datestring = '%04i%02i%02i' % (processdate.year, processdate.month, processdate.day)
    
    # get mask defined by provided options
    qc_flags = load_flags(mask_directory, iteration, datestring).data[:,:,:]
    
    if mask_directory2 is not None:
        qc_flags = qc_flags | load_flags(mask_directory2, iteration, datestring).data[:,:,:]
    
    # set a new marine flag with a more flexible definition of marine (95% water rather than 100%)
    qc_flags = reset_marine_flag(qc_flags)
    
    combined_mask = compute_mask(qc_flags, flags_to_mask)
    
    # get analysis
    analysis, uncertainty, influence, ensemble_members = load_analysis(analysis_directory, iteration, datestring)
    
    # build masked output product file
    #outputfile = output_directory+'/eustace_analysis_'+str(iteration)+'/'+datestring[:4]+'/eustace_analysis_'+str(iteration)+'_masked_'+datestring+'.nc'
    
    #format as variable_collection_framework_realization_YYYYmmdd.nc so tas_global_eustace_0_YYYYmmdd.nc 
    if iteration == 4:
        version_string = "R001400"
        rundate = '20190326'
    elif iteration == 9:
        version_string = "R001413"
        rundate = '20190405'
    else:
        print "unrecognised version"
        
    outputfile = output_directory+'/'+version_string+'/'+rundate+'/global/'+datestring[:4]+'/tas_global_eustace_0_'+datestring+'.nc'
    print outputfile
    filebuilder = FileBuilderGlobalField(
        outputfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        'Met Office',
        '',
        '', 
        '')
    numpy.squeeze(analysis.data[:,:,:])
    numpy.squeeze(uncertainty.data[:,:,:])
    numpy.squeeze(influence.data[:,:,:])
    filebuilder.add_global_field(definitions.TAS, numpy.ma.masked_where( combined_mask, analysis.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, numpy.ma.masked_where( combined_mask, uncertainty.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, numpy.ma.masked_where( combined_mask, influence .data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    # sort ensemble members and apply mask
    ensemble_members = sorted(ensemble_members, key=lambda item: int(item.var_name.split('_')[-1]))
    
    for index, ensemble_member in enumerate(ensemble_members):
        variable = copy.deepcopy(definitions.TASENSEMBLE)
        variable.name = variable.name + '_' + str(index)
        filebuilder.add_global_field(variable, numpy.ma.masked_where( combined_mask, ensemble_member.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    filebuilder.save_and_close()
    
    # climatology only output
    analysis, uncertainty, influence = load_climatology(analysis_directory, iteration, datestring)
    climatologyfile = output_directory+'/eustace_climatology_component_'+str(iteration)+'/'+datestring[:4]+'/eustace_climatology_component_'+str(iteration)+'_masked_'+datestring+'.nc'
    
    filebuilder = FileBuilderGlobalField(
        climatologyfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional component output - climatology',
        __name__, 
        '')
    
    numpy.squeeze(analysis.data[:,:,:])
    numpy.squeeze(uncertainty.data[:,:,:])
    numpy.squeeze(influence.data[:,:,:])
    filebuilder.add_global_field(definitions.TAS, numpy.ma.masked_where( combined_mask, analysis.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, numpy.ma.masked_where( combined_mask, uncertainty.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, numpy.ma.masked_where( combined_mask, influence .data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    filebuilder.save_and_close()
    
    # large scale only output
    analysis, uncertainty, influence = load_large_scale(analysis_directory, iteration, datestring)
    largescalefile = output_directory+'/eustace_large_scale_component_'+str(iteration)+'/'+datestring[:4]+'/eustace_large_scale_component_'+str(iteration)+'_masked_'+datestring+'.nc'
    filebuilder = FileBuilderGlobalField(
        largescalefile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional component output - large scale',
        __name__, 
        '')
    
    numpy.squeeze(analysis.data[:,:,:])
    numpy.squeeze(uncertainty.data[:,:,:])
    numpy.squeeze(influence.data[:,:,:])
    filebuilder.add_global_field(definitions.TASPERTURBATION, numpy.ma.masked_where( combined_mask, analysis.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, numpy.ma.masked_where( combined_mask, uncertainty.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, numpy.ma.masked_where( combined_mask, influence .data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    
    filebuilder.save_and_close()
    
    # local only output
    analysis, uncertainty, influence = load_local(analysis_directory, iteration, datestring)
    localfile = output_directory+'/eustace_local_component_'+str(iteration)+'/'+datestring[:4]+'/eustace_local_component_'+str(iteration)+'_masked_'+datestring+'.nc'
    filebuilder = FileBuilderGlobalField(
        localfile, 
        eustace.timeutils.epoch.days_since_epoch(processdate),
        'EUSTACE Analysis',
        get_revision_id_for_module(eustace),
        definitions.TAS.name,
        '',
        'Provisional component output - local',
        __name__, 
        '')

    numpy.squeeze(analysis.data[:,:,:])
    numpy.squeeze(uncertainty.data[:,:,:])
    numpy.squeeze(influence.data[:,:,:])
    filebuilder.add_global_field(definitions.TASPERTURBATION, numpy.ma.masked_where( combined_mask, analysis.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TASUNCERTAINTY, numpy.ma.masked_where( combined_mask, uncertainty.data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))
    filebuilder.add_global_field(definitions.TAS_OBSERVATION_INFLUENCE, numpy.ma.masked_where( combined_mask, influence .data[:,:,:] ).reshape(definitions.GLOBAL_FIELD_SHAPE))


    filebuilder.save_and_close()

"""

Calls setup to run as batches in lsf

"""

def masking_operation(  reference_time_string,
                        operation_index,
                        analysis_directory,
                        iteration,
                        flag_directory,
                        flag_directory2,
                        output_directory ):
    
    # get dates in the month to be processed in this operation
    reference_time = parser.parse(reference_time_string)
    processing_dates = operation_dates(reference_time, operation_index)
    print processing_dates
    # derive flags for each day
    for date in processing_dates:
        
        year = date.year
        month = date.month
        day = date.day
        
        datestring = '%04i%02i%02i' % (year, month, day)
        
        flags_to_mask = copy.deepcopy( FLAGS_TO_MASK )
        # mask marine regions post 2012
        if year > 2012:
            flags_to_mask.append( [MISSING_MARINE_FLAG, NULL_FLAG] )
        
        apply_mask(analysis_directory, iteration, flag_directory, flag_directory2, flags_to_mask, date, output_directory)



def main():

    print 'Submission of advanced standard analysis masking jobs'
    
    parser = argparse.ArgumentParser(description='Monthly batches of masking operations')
    
    parser.add_argument('--reference_time_string', type=str, default="1880-01-01", help='reference first day to run formatted YYY-MM-DD. Should be first day of month.')
    parser.add_argument('--operation_index',  type=int, default=0, help='the number of months since the reference_time_string that this masking operation corresponds to')
    parser.add_argument('--analysis_directory', type=str, default="/work/scratch/cmorice/advanced_standard/", help='directory in which the input analysis files can be found')
    parser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    parser.add_argument('--flag_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='input location for the combined qc flag file')
    parser.add_argument('--flag_directory2',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='input location for the combined qc flag file')
    parser.add_argument('--output_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='output location for the masked analysis grids')

    args = parser.parse_args()

    masking_operation(  args.reference_time_string,
                        args.operation_index,
                        args.analysis_directory,
                        args.iteration,
                        args.flag_directory,
                        args.flag_directory2,
                        args.output_directory)

if __name__ == '__main__':
    
    main()
