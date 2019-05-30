"""Extreme and areal value QC checks on EUSTACE analysis products"""

import iris
import iris.plot as iplt

import numpy
import copy
from eustace.outputformats.outputvariable import OutputVariable
from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eumopps.version.svn import get_revision_id_for_module
from netCDF4 import default_fillvals
from eustace.outputformats import definitions
import eustace.timeutils.epoch

from dateutil import parser
import datetime
import os.path

import argparse

from apply_mask import load_analysis, save_flag_file
from operation_count import operation_dates
import flags

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
CALENDAR_DAY_FLAG    = flags.CALENDAR_DAY_FLAG
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

# missing data indicator
MISSING_FLAG_FLAG    = flags.MISSING_FLAG_FLAG

FLAG_MEANINGS = flags.FLAG_MEANINGS

# parameters used in areal stastics
DEFAULT_EARTH_RADIUS_KM = 6371.0
COMPUTATION_BLOCKSIZE = 10          # number of grid cells in a dimension to be computed simultaneously
                                    # larger numbers increase memory requirement
                                    # smaller number increase computation time
                                    # has no effect in output

SEARCH_RANGE_KM = 750.
MAX_LATITUDE_SEARCH_RANGE_GRID_CELLS=50
MAX_LONGITUDE_SEARCH_RANGE_GRID_CELLS = 50
LATITUDE_BLOCK_SIZE = 20

# pre-extreme calculation lapse rate adjustment 
LAPSE_RATE   = -4.5                 # in K km^-1 (taken from analysis fit)
METRES_TO_KM = 0.001                # multiplier to convert m to km
ALTITUDE_FILE='/gws/nopw/j04/eustace/data/internal/climatology_covariates/DEM_global_0.25_0.25.nc'
ALTITUDE_FILE_VAR = 'dem'
COAST_FILE='/gws/nopw/j04/eustace/data/internal/climatology_covariates/coastal_influence.test.0.25_0.25.nc'
COAST_FILE_VAR = 'coastal_influence'

# parameters for value qc
QUARTILE_LIMIT = 4.0
EXTREME_LOWER_LIMIT=193.15          # lowest permited temperature in kelvin
EXTREME_UPPER_LIMIT=333.15          # highest permited temperature in kelvin
CLIMATOLOGY_UNC_SMOOTH_WINDOW = 13  # window for smoothing climatology uncertainties before thresholding in number of grid cells
LARGE_SCALE_UNC_SMOOTH_WINDOW = 17  # window for smoothing large-scale uncertainties before thresholding in number of grid cells
OUTLIER_MASK_WINDOW = 5             # window for extending the region masked in outlier checks in number of grid cells

# omitted data source parameters
FINAL_MARINE_YEAR = 2012

"""

Arc length to angle conversions

"""


def latitude_delta(north_range_km):
    # km to degree conversions    
    return numpy.rad2deg( north_range_km / DEFAULT_EARTH_RADIUS_KM )

def longitude_delta(east_range_km, latitude):
    # km to degree conversions    
    return numpy.rad2deg( east_range_km / ( DEFAULT_EARTH_RADIUS_KM * numpy.sin(numpy.pi/ 2.0 - numpy.deg2rad( latitude ) ) ) )

"""

Efficient areal statistics calculation

"""

def compute_window_median_quartiles(image, window_latitude, window_longitude, padmode_latitude='reflect', padmode_longitude='wrap'):
    """Running median and quartiles on a 2D grid"""
    image_shape = image.shape

    image = numpy.pad(image, ( (window_latitude//2, window_latitude//2), (0, 0) ), padmode_latitude)
    image = numpy.pad(image, ( (0, 0), (window_longitude//2, window_longitude//2) ), padmode_longitude)
    
    n_latitude_cells, n_longitude_cells = image.shape
    strided_image = numpy.lib.stride_tricks.as_strided( image, 
                                                        shape=[n_latitude_cells - window_latitude + 1, n_longitude_cells - window_longitude + 1, window_latitude, window_longitude],
                                                        strides=image.strides + image.strides)
    
    # important: trying to reshape image will create complete 4-dimensional copy
    
    # operate on blocks to prevent excessive memory use
    blocksize_latitude  = COMPUTATION_BLOCKSIZE
    blocksize_longitude = COMPUTATION_BLOCKSIZE

    
    lower_quartile = numpy.zeros(image_shape)
    median = numpy.zeros(image_shape)
    upper_quartile  = numpy.zeros(image_shape)

    latitude_steps = int(numpy.ceil(float(image_shape[0]) / blocksize_latitude))
    longitude_steps = int(numpy.ceil(float(image_shape[1]) / blocksize_longitude))

    for latitude_step in range(latitude_steps):
        for longitude_step in range(longitude_steps):
            if latitude_step == latitude_steps-1:
                latitude_slice=slice(latitude_step * blocksize_latitude, None)
            else:
                latitude_slice=slice(latitude_step * blocksize_latitude, (latitude_step + 1) * blocksize_latitude)
            
            if longitude_step == longitude_steps-1:
                longitude_slice=slice(longitude_step * blocksize_longitude, None)
            else:
                longitude_slice=slice(longitude_step * blocksize_longitude, (longitude_step + 1) * blocksize_longitude)

            image_block = strided_image[latitude_slice,longitude_slice,:,:]
            
            lower_quartile_block, median_block, upper_quartile_block = numpy.percentile(image_block, [25, 50, 75], axis=(2,3))
            
            lower_quartile[latitude_slice,longitude_slice] = lower_quartile_block
            median[latitude_slice,longitude_slice] = median_block
            upper_quartile[latitude_slice,longitude_slice] = upper_quartile_block

    return lower_quartile, median, upper_quartile 


def compute_window_mean(image, window_latitude, window_longitude, padmode_latitude='reflect', padmode_longitude='wrap'):
    """Running mean on a 2D grid"""
    
    image_shape = image.shape
    
    image = numpy.pad(image, ( (window_latitude//2, window_latitude//2), (0, 0) ), padmode_latitude)
    image = numpy.pad(image, ( (0, 0), (window_longitude//2, window_longitude//2) ), padmode_longitude)
    
    n_latitude_cells, n_longitude_cells = image.shape
    strided_image = numpy.lib.stride_tricks.as_strided( image, 
                                                        shape=[n_latitude_cells - window_latitude + 1, n_longitude_cells - window_longitude + 1, window_latitude, window_longitude],
                                                        strides=image.strides + image.strides)
    
    # important: trying to reshape image will create complete 4-dimensional copy
    
    # operate on blocks to prevent excessive memory use
    blocksize_latitude = 10
    blocksize_longitude = 10
    
    rolling_mean = numpy.zeros(image_shape)
    
    latitude_steps = image_shape[0] // blocksize_latitude
    longitude_steps = image_shape[1] // blocksize_longitude
    
    
    for latitude_step in range(latitude_steps):
        for longitude_step in range(longitude_steps):
            if latitude_step == latitude_steps-1:
                latitude_slice=slice(latitude_step * blocksize_latitude, None)
            else:
                latitude_slice=slice(latitude_step * blocksize_latitude, (latitude_step + 1) * blocksize_latitude)
            
            if longitude_step == longitude_steps-1:
                longitude_slice=slice(longitude_step * blocksize_longitude, None)
            else:
                longitude_slice=slice(longitude_step * blocksize_longitude, (longitude_step + 1) * blocksize_longitude)

            image_block = strided_image[latitude_slice,longitude_slice,:,:]
            
            rolling_mean_block = numpy.mean(image_block, axis=(2,3))
            rolling_mean[latitude_slice,longitude_slice] = rolling_mean_block

    return rolling_mean 

def block_iqr_by_range(latitudes, longitudes, data_array, range_km, max_latitude_cell_range=20, max_longitude_cell_range = 20, latitude_block_size = 20):
    """
    
    Compute running median and interquartile range with window size governed by
    approximate range in km up to a maximum extent.    
    
    data_array must have dimensions ordered (latitudes, longitudes)
    
    """
    
    # initialise output arrays
    lower_quartile = numpy.zeros(data_array.shape) + numpy.nan
    median = numpy.zeros(data_array.shape) + numpy.nan
    upper_quartile = numpy.zeros(data_array.shape) + numpy.nan
    
    # get number of blocks
    n_blocks = int( numpy.ceil(len(latitudes) / latitude_block_size) )
    
    # get grid resolutions
    latitude_grid_resolution = latitudes[1] - latitudes[0]
    longitude_grid_resolution = longitudes[1] - longitudes[0]
    
    # iterate over blocks computing filter
    for block_ind in range(n_blocks):
        
        lat_ind_0 = block_ind * latitude_block_size
        lat_ind_1 = min((block_ind+1) * latitude_block_size, len(latitudes) )
        
        centre_latitude = numpy.mean(latitudes[lat_ind_0:lat_ind_1])
        
        # compute window distance for filter for this block
        latitude_range_grid = numpy.int( numpy.ceil( latitude_delta(range_km) / latitude_grid_resolution ) )
        latitude_range_grid = min(latitude_range_grid, max_latitude_cell_range)
        window_latitude = 2 * latitude_range_grid + 1
        
        longitude_range_grid = numpy.int( numpy.ceil( longitude_delta(range_km, centre_latitude) / longitude_grid_resolution ) )
        longitude_range_grid = min(longitude_range_grid, max_longitude_cell_range)
        window_longitude = 2 * longitude_range_grid + 1

        print "centre_latitude: ", centre_latitude
        print window_latitude, window_longitude

        # extract extended block with latitude margin for window width
        extended_lat_ind_0 = max(lat_ind_0 - latitude_range_grid, 0 )
        extended_lat_ind_1 = min((block_ind+1) * latitude_block_size+latitude_range_grid, len(latitudes) )
        print extended_lat_ind_0, extended_lat_ind_1
        block_data_array = data_array[extended_lat_ind_0:extended_lat_ind_1,:]
        
        # compute and store statistics for window
        print block_data_array.shape
        block_lower_quartile, block_median, block_upper_quartile = compute_window_median_quartiles(block_data_array, window_latitude, window_longitude, padmode_latitude='reflect', padmode_longitude='wrap')
        print "block_lower_quartile.shape:", block_lower_quartile.shape
        
        extract_lat_ind_0 = lat_ind_0-extended_lat_ind_0
        extract_lat_ind_1 = extract_lat_ind_0 + (lat_ind_1 - lat_ind_0)
        
        print lat_ind_0, lat_ind_1
        print extended_lat_ind_0, extended_lat_ind_1
        print extract_lat_ind_0, extract_lat_ind_1
        
        lower_quartile[lat_ind_0:lat_ind_1,:] = block_lower_quartile[extract_lat_ind_0:extract_lat_ind_1,:]
        median[lat_ind_0:lat_ind_1,:] = block_median[extract_lat_ind_0:extract_lat_ind_1,:]
        upper_quartile[lat_ind_0:lat_ind_1,:] = block_upper_quartile[extract_lat_ind_0:extract_lat_ind_1,:]
    
    return lower_quartile, median, upper_quartile

"""

Core methods for flag calculation

"""

def altitude_adjustment(analysis_field):
    """Remove the effects of altitude before running the areal QC"""
    
    variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == ALTITUDE_FILE_VAR))
    altitude_map = iris.load_cube( ALTITUDE_FILE, variable_constraint ).data[:,:] * METRES_TO_KM
    adjusted_analysis_field = analysis_field - LAPSE_RATE * altitude_map 
    
    return adjusted_analysis_field

def extremes_checks(analysis_directory, output_directory, iteration, processdate, extreme_lower_limit, extreme_upper_limit):
    """check upper and lower extreme temperature limits"""
    
    # date time filename formated datestring for day on interest
    datestring = "{:04d}{:02d}{:02d}".format(processdate.year, processdate.month, processdate.day )
    print processdate
    
    # get the component field for this date
    try:
        analysis, uncertainty, influence, ensemble = load_analysis(analysis_directory, iteration, datestring)
    except:
        print "Error loading large scale file" 

    analysis_field = numpy.squeeze( analysis.data )

    # detect exceedances
    extreme_low_value_exceedance = analysis_field < extreme_lower_limit
    extreme_high_value_exceedance = analysis_field > extreme_upper_limit
    
    # extend flagged regions spatially
    extreme_low_value_exceedance = numpy.int64( 1.0 - compute_window_mean(extreme_low_value_exceedance, OUTLIER_MASK_WINDOW, OUTLIER_MASK_WINDOW, padmode_latitude='reflect', padmode_longitude='wrap') ) == 0
    extreme_high_value_exceedance = numpy.int64( 1.0 - compute_window_mean(extreme_high_value_exceedance, OUTLIER_MASK_WINDOW, OUTLIER_MASK_WINDOW, padmode_latitude='reflect', padmode_longitude='wrap') ) == 0
    
    # set flag values for output
    flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )    
    
    flag_values[extreme_low_value_exceedance] = flag_values[extreme_low_value_exceedance] | EXTREME_LOW_FLAG
    flag_values[extreme_high_value_exceedance] = flag_values[extreme_high_value_exceedance] | EXTREME_HIGH_FLAG
    
    return flag_values

def areal_checks(analysis_directory, output_directory, iteration, processdate, quartile_limit):
    """Compute a quantile based areal qc check to detect regional outliers"""
    
    # date time filename formated datestring for day on interest
    datestring = "{:04d}{:02d}{:02d}".format(processdate.year, processdate.month, processdate.day )
    print processdate
    
    # get the component field for this date
    try:
        analysis, uncertainty, influence, ensemble = load_analysis(analysis_directory, iteration, datestring)
    except:
        print "Error loading large scale file" 

    analysis_field = numpy.squeeze( analysis.data )
    
    latitudes = analysis.coord('latitude').points
    longitudes = analysis.coord('longitude').points
    
    # remove the altitude effect before masking to reduce detection of altitude related effects
    adjusted_analysis_field = altitude_adjustment(analysis_field)
    
    # compute median and quantiles across the field
    lower_quartile, median, upper_quartile = block_iqr_by_range(latitudes,
                                                                longitudes,
                                                                adjusted_analysis_field,
                                                                SEARCH_RANGE_KM,
                                                                max_latitude_cell_range=MAX_LATITUDE_SEARCH_RANGE_GRID_CELLS,
                                                                max_longitude_cell_range = MAX_LONGITUDE_SEARCH_RANGE_GRID_CELLS,
                                                                latitude_block_size = LATITUDE_BLOCK_SIZE)
    
    # compute detection limits at each location in the field
    areal_lower_limit = lower_quartile - quartile_limit * (median - lower_quartile)
    areal_upper_limit = upper_quartile + quartile_limit * (upper_quartile - median)
    
    # flag exceedances
    areal_low_value_exceedance = adjusted_analysis_field < areal_lower_limit
    areal_high_value_exceedance = adjusted_analysis_field > areal_upper_limit
    
    # extend the flagged regions spatially
    areal_low_value_exceedance = numpy.int64( 1.0 - compute_window_mean(areal_low_value_exceedance, OUTLIER_MASK_WINDOW, OUTLIER_MASK_WINDOW, padmode_latitude='reflect', padmode_longitude='wrap') ) == 0
    areal_high_value_exceedance = numpy.int64( 1.0 - compute_window_mean(areal_high_value_exceedance, OUTLIER_MASK_WINDOW, OUTLIER_MASK_WINDOW, padmode_latitude='reflect', padmode_longitude='wrap') ) == 0
    
    # set flag values for output
    flag_values = numpy.zeros( definitions.GLOBAL_FIELD_SHAPE[1:], FLAG_TYPE )
    
    flag_values[areal_low_value_exceedance] = flag_values[areal_low_value_exceedance] | AREAL_LOW_FLAG
    flag_values[areal_high_value_exceedance] = flag_values[areal_high_value_exceedance] | AREAL_HIGH_FLAG
    
    return flag_values

def flagging_operation(  reference_time_string,
                         operation_index,
                         analysis_directory,
                         iteration,
                         output_directory, ):
    
    # get dates in the month to be processed in this operation
    reference_time = parser.parse(reference_time_string)
    processing_dates = operation_dates(reference_time, operation_index)
    
    # derive flags for each day
    for processdate in processing_dates:
        
        # check that directory exists and if not then make it
        if not os.path.exists(os.path.join(output_directory, str(processdate.year))):
            os.makedirs(os.path.join(output_directory, str(processdate.year)))
        
        # compute flags
        areal_flag_values = areal_checks(analysis_directory, output_directory, iteration, processdate, QUARTILE_LIMIT)
        extreme_flag_values = extremes_checks(analysis_directory, output_directory, iteration, processdate, EXTREME_LOWER_LIMIT, EXTREME_UPPER_LIMIT)
        
        # join flags
        flag_values = areal_flag_values | extreme_flag_values
        
        #save
        outputfile = os.path.join(output_directory, '{:04d}'.format(processdate.year), 'eustace_analysis_{:d}_qc_flags_{:04d}{:02d}{:02d}.nc'.format(iteration, processdate.year, processdate.month, processdate.day))
        save_flag_file(flag_values, processdate, outputfile)

def main():

    print 'Submission of advanced standard analysis jobs'
    
    argparser = argparse.ArgumentParser(description='Monthly batches of value QC flagging operations')
    
    argparser.add_argument('--reference_time_string', type=str, default="1850-01-01", help='reference first day to run formatted YYY-MM-DD. Should be first day of month.')
    argparser.add_argument('--operation_index',  type=int, default=0, help='the number of months since the reference_time_string that this masking operation corresponds to')
    argparser.add_argument('--analysis_directory', type=str, default="/work/scratch/cmorice/advanced_standard/", help='root directory of the analysis')
    argparser.add_argument('--iteration',  type=int, default=9, help='operation index at which the analysis grid produced')
    argparser.add_argument('--output_directory',  type=str, default="/gws/nopw/j04/eustace_vol2/masking/", help='operation index at which the analysis grid produced')
    
    args = argparser.parse_args()

    flagging_operation(  args.reference_time_string,
                         args.operation_index,
                         args.analysis_directory,
                         args.iteration,
                         args.output_directory)

if __name__ == '__main__':
    
    main()
    
    
