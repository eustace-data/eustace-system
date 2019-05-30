"""Definitions for creating output files according to EUSTACE specification."""

from outputvariable import OutputVariable
from outputvariable import OutputVariableTemplate
from eustace.timeutils.epoch import EPOCH
from netCDF4 import default_fillvals
from netCDF4 import stringtochar
import numpy

__version__ = "$Revision: 1369 $"
__author__ = "Joel R. Mitchelson"

FILE_FORMAT = 'NETCDF4'
FILE_CONVENTIONS = 'CF-1.6'

TITLE_PREFIX = 'EUSTACE Surface Air Temperature Estimates '
SYSTEM_PREFIX = 'EUSTACE System '

FRAMEWORK_NAME = 'eustace'

GLOBAL_FIELD_SUBDIRECTORY_PATTERN = '{YYYY}'
GLOBAL_FIELD_FILENAME_PATTERN = '{variable}_{collection}_{framework}_{realisation}_{YYYY}{mm}{dd}.nc'
GLOBAL_FIELD_SHAPE = (1, 720, 1440)
GLOBAL_SAMPLE_SHAPE = (1, 720, 1440, 10)
GLOBAL_FIELD_RESOLUTION = 0.25
GLOBAL_FIELD_OUTPUT_FLAGS = ['MAP', 'post_STD', 'prior_STD']

STATIONS_TEMPERATURE_FILENAME_PATTERN = 'eustace_stations_{region}_{version}_temperature.nc'
STATIONS_STATUS_FILENAME_PATTERN = 'eustace_stations_{region}_{version}_status.nc'
STATIONS_TEMPERATURE_TITLE_SUFFIX = 'at Stations'
STATIONS_STATUS_TITLE_SUFFIX = 'Station Status'
STATIONS_NAME_STRLEN = 32
STATIONS_TEMPERATURE_FEATURE_TYPE_NAME = 'featureType'
STATIONS_TEMPERATURE_FEATURE_TYPE_VALUE = 'timeSeries'

DIMENSION_NAME_TIME = 'time'
DIMENSION_NAME_LATITUDE = 'latitude'
DIMENSION_NAME_LONGITUDE = 'longitude'
DIMENSION_NAME_STATION = 'station'
DIMENSION_NAME_NAME_STRLEN = 'name_strlen'
DIMENSION_NAME_DETECTION_TIME = 'detection_time'
DIMENSION_NAME_TASMIN_BREAK = 'tasmin_break'
DIMENSION_NAME_TAS_BREAK = 'tas_break'
DIMENSION_NAME_TASMAX_BREAK = 'tasmax_break'
DIMENSION_NAME_MERGED_BREAK = 'merged_break'
DIMENSION_NAME_BOUNDS = 'bounds'
DIMENSION_NAME_ELEVATION = 'elevation'

TEMPERATURE_DATA_TYPE = numpy.int16
TEMPERATURE_FILL_VALUE = -32768
TEMPERATURE_SCALE_FACTOR = 0.005
TEMPERATURE_ADD_OFFSET = 273.15

TEMPERATURE_UNCERTAINTY_ADD_OFFSET = 0.0
TEMPERATURE_UNCERTAINTY_SCALE_FACTOR = 0.001
TEMPERATURE_UNCERTAINTY_FILL_VALUE = -32768

CLIMATOLOGY_FRACTION_DATA_TYPE = numpy.int16
CLIMATOLOGY_FRACTION_FILL_VALUE = -32768
CLIMATOLOGY_FRACTION_SCALE_FACTOR = 1.5625e-05
CLIMATOLOGY_FRACTION_ADD_OFFSET = 0.5

OBSERVATION_INFLUENCE_DATA_TYPE = numpy.int16
OBSERVATION_INFLUENCE_FILL_VALUE = -32768
OBSERVATION_INFLUENCE_SCALE_FACTOR = 1.5625e-05
OBSERVATION_INFLUENCE_ADD_OFFSET = 0.5

ENSEMBLE_DATA_TYPE = numpy.int16
ENSEMBLE_FILL_VALUE = -32768
ENSEMBLE_SCALE_FACTOR = 0.005
ENSEMBLE_ADD_OFFSET = 273.15

PERTURBATION_DATA_TYPE = numpy.int16
PERTURBATION_FILL_VALUE = -32768
PERTURBATION_SCALE_FACTOR = 0.005
PERTURBATION_ADD_OFFSET = 0.0

TEMPLATE_TEMPERATURE = OutputVariableTemplate(
    dtype=TEMPERATURE_DATA_TYPE,
    fill_value=TEMPERATURE_FILL_VALUE,
    template_long_name='{quantity} daily surface air temperature',
    scale_factor=TEMPERATURE_SCALE_FACTOR,
    add_offset=TEMPERATURE_ADD_OFFSET,
    standard_name='air_temperature',
    units='K')

TEMPLATE_TEMPERATURE_UNCERTAINTY = OutputVariableTemplate(
    dtype=TEMPERATURE_DATA_TYPE,
    fill_value=TEMPERATURE_UNCERTAINTY_FILL_VALUE,
    template_long_name='Total uncertainty in {quantity} daily surface air temperature',
    scale_factor=TEMPERATURE_UNCERTAINTY_SCALE_FACTOR,
    add_offset=TEMPERATURE_UNCERTAINTY_ADD_OFFSET,
    units='K')
    
TEMPLATE_PERTURBATION = OutputVariableTemplate(
    dtype=PERTURBATION_DATA_TYPE,
    fill_value=PERTURBATION_FILL_VALUE,
    template_long_name='{quantity} daily surface air temperature',
    scale_factor=PERTURBATION_SCALE_FACTOR,
    add_offset=PERTURBATION_ADD_OFFSET,
    standard_name='air_temperature',
    units='K')

TEMPLATE_TEMPERATURE_UNCERTAINTY_TOTAL = OutputVariableTemplate.extend(
    TEMPLATE_TEMPERATURE_UNCERTAINTY,
    template_long_name='Total uncertainty in {quantity} daily surface air temperature')

TEMPLATE_CLIMATOLOGY_FRACTION = OutputVariableTemplate(
    dtype=CLIMATOLOGY_FRACTION_DATA_TYPE,
    fill_value=CLIMATOLOGY_FRACTION_FILL_VALUE,
    template_long_name='Climatology fraction for {quantity} daily surface air temperature',
    scale_factor=CLIMATOLOGY_FRACTION_SCALE_FACTOR,
    add_offset=CLIMATOLOGY_FRACTION_ADD_OFFSET)

TEMPLATE_OBSERVATION_INFLUENCE = OutputVariableTemplate(
    dtype=OBSERVATION_INFLUENCE_DATA_TYPE,
    fill_value=OBSERVATION_INFLUENCE_FILL_VALUE,
    template_long_name='Observation influence for {quantity} daily surface air temperature',
    scale_factor=OBSERVATION_INFLUENCE_SCALE_FACTOR,
    add_offset=OBSERVATION_INFLUENCE_ADD_OFFSET)

TEMPLATE_ENSEMBLE = OutputVariableTemplate(
    dtype=ENSEMBLE_DATA_TYPE,
    fill_value=ENSEMBLE_FILL_VALUE,
    template_long_name='{quantity} daily surface air temperature ensemble',
    scale_factor=ENSEMBLE_SCALE_FACTOR,
    add_offset=ENSEMBLE_ADD_OFFSET,
    units='K')

TIME_BOUNDS = OutputVariable(
    name='timebounds',
    dtype=numpy.float32,
    fill_value=None)

TIME_UNITS_DAYS_SINCE_EPOCH = 'days since {0:04d}-{1:02d}-{2:02d}T00:00:00Z'.format(EPOCH.year, EPOCH.month, EPOCH.day)

TIME = OutputVariable(
    name='time',
    dtype=numpy.float32,
    fill_value=None,
    standard_name='time',
    long_name='Time at zero longitude',
    units=TIME_UNITS_DAYS_SINCE_EPOCH,
    calendar='gregorian',
    bounds=TIME_BOUNDS.name,
    ancillary_variables='timeoffset',
    axis='T')

TIME_OFFSET = OutputVariable(
    name='timeoffset',
    dtype=numpy.float32,
    fill_value=None,
    long_name='Local time offset from UTC (days)',
    units='days')

LATITUDE = OutputVariable(
    name='latitude',
    dtype=numpy.float32,
    fill_value=None,
    standard_name='latitude',
    long_name='Latitude (deg)',
    units='degrees_north',
    axis='Y')

LONGITUDE = OutputVariable(
    name='longitude',
    dtype=numpy.float32,
    fill_value=None,
    standard_name='longitude',
    long_name='Longitude (deg)',
    units='degrees_east',
    axis='X')

ELEVATION = OutputVariable(
    name='elevation',
    dtype=numpy.float32,
    fill_value=-999.9,
    standard_name='surface_altitude',
    long_name='Height above the geoid (m)',
    units='m')

TAS = OutputVariable.from_template(TEMPLATE_TEMPERATURE, 'tas', quantity='average', cell_methods='time: mean')
TASMIN = OutputVariable.from_template(TEMPLATE_TEMPERATURE, 'tasmin', quantity='minimum',  cell_methods='time: minimum')
TASMAX = OutputVariable.from_template(TEMPLATE_TEMPERATURE, 'tasmax', quantity='maximum', cell_methods='time: maximum')

TASUNCERTAINTY = OutputVariable.from_template(TEMPLATE_TEMPERATURE_UNCERTAINTY_TOTAL, 'tasuncertainty', quantity='average')
TASMINUNCERTAINTY = OutputVariable.from_template(TEMPLATE_TEMPERATURE_UNCERTAINTY_TOTAL, 'tasminuncertainty', quantity='mimimum')
TASMAXUNCERTAINTY = OutputVariable.from_template(TEMPLATE_TEMPERATURE_UNCERTAINTY_TOTAL, 'tasmaxuncertainty', quantity='maximum')

TAS_CLIMATOLOGY_FRACTION = OutputVariable.from_template(TEMPLATE_CLIMATOLOGY_FRACTION, 'tasclimatologyfraction', quantity='average')
TASMIN_CLIMATOLOGY_FRACTION = OutputVariable.from_template(TEMPLATE_CLIMATOLOGY_FRACTION, 'tasminclimatologyfraction', quantity='minimum')
TASMAX_CLIMATOLOGY_FRACTION = OutputVariable.from_template(TEMPLATE_CLIMATOLOGY_FRACTION, 'tasmaxclimatologyfraction', quantity='maximum')

TAS_OBSERVATION_INFLUENCE = OutputVariable.from_template(TEMPLATE_OBSERVATION_INFLUENCE, 'tasobservationinfluence', quantity='average')
TASMIN_OBSERVATION_INFLUENCE = OutputVariable.from_template(TEMPLATE_OBSERVATION_INFLUENCE, 'tasminobservationinfluence', quantity='minimum')
TASMAX_OBSERVATION_INFLUENCE = OutputVariable.from_template(TEMPLATE_OBSERVATION_INFLUENCE, 'tasmaxobservationinfluence', quantity='maximum')

TASENSEMBLE = OutputVariable.from_template(TEMPLATE_ENSEMBLE, 'tasensemble', quantity='average', cell_methods='time: mean')
TASMINENSEMBLE = OutputVariable.from_template(TEMPLATE_ENSEMBLE, 'tasminensemble', quantity='minimum',  cell_methods='time: minimum')
TASMAXENSEMBLE = OutputVariable.from_template(TEMPLATE_ENSEMBLE, 'tasmaxensemble', quantity='maximum', cell_methods='time: maximum')

TASPERTURBATION = OutputVariable.from_template(TEMPLATE_PERTURBATION, 'tas_perturbation', quantity='average', cell_methods='time: mean')
TASMINPERTURBATION = OutputVariable.from_template(TEMPLATE_PERTURBATION, 'tasmin_perturbation', quantity='average', cell_methods='time: mean')
TASMAXPERTURBATION = OutputVariable.from_template(TEMPLATE_PERTURBATION, 'tasmax_perturbation', quantity='average', cell_methods='time: mean')

SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY = OutputVariableTemplate.extend(
    TEMPLATE_TEMPERATURE_UNCERTAINTY,
    template_long_name='{component} uncertainty on {quantity} daily surface air temperature')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMIN = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_rand',
    component='random', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMIN = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_corr_atm',
    component='Locally correlated atmospheric', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMIN = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_corr_sfc',
    component='Locally correlated surface', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCSYSTASMIN = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_sys',
    component='systematic', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMAX = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_rand',
    component='random', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMAX = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_corr_atm',
    component='Locally correlated atmospheric', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMAX = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_corr_sfc',
    component='Locally correlated surface', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_LAND_UNCSYSTASMAX = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_LAND_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_sys',
    component='systematic', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_LAND_TMINMODELNUMBER = OutputVariable(
    name='tasmin_model_number',
    long_name='Model number used for estimating Tmin from satellite data',
    dtype=numpy.int16,
    fill_value=numpy.int16(0),
    scale_factor=numpy.int16(1),
    add_offset=numpy.int16(0),
    flag_values=numpy.array([1, 2, 3], numpy.uint8),
    flag_meanings="""Model 1: include both LST-day and LST-night to predict tasmin\n
		   Model 2: include LST-night to predict tasmin\n
		   Model 3: include LST-day to predict tasmin\n""")

SURFACEAIRMODEL_ANCILLARY_LAND_TMAXMODELNUMBER = OutputVariable(
    name='tasmax_model_number',
    long_name='Model number used for estimating Tmax from satellite data',
    dtype=numpy.int16,
    fill_value=numpy.int16(0),
    scale_factor=numpy.int16(1),
    add_offset=numpy.int16(0),
    flag_values=numpy.array([1, 2, 3], numpy.uint8),
    flag_meanings="""Model 1: include both LST-day and LST-night to predict tasmax\n
		   Model 2: include LST-day to predict tasmax\n
		   Model 3: include LST-night to predict tasmax\n""")

OLD_LAND_ANCILLARY_NAMES = ['unc_rand_tasmin', 'unc_corr_atm_tasmin', 'unc_corr_sfc_tasmin', 'unc_sys_tasmin', 
                            'unc_rand_tasmax', 'unc_corr_atm_tasmax', 'unc_corr_sfc_tasmax', 'unc_sys_tasmax']

SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY = OutputVariableTemplate.extend(
    TEMPLATE_TEMPERATURE_UNCERTAINTY,
    template_long_name='{component} on {quantity} daily surface air temperature')

SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_EXCLUDINGCLOUD = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tas_unc_no_cloud',
    component='Total uncertainty excluding cloud', quantity='average')

SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_UNCORRELATED = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tas_unc_rand',
    component='Random uncertainty', quantity='average')

SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_SYNOPTICALLYCORRELATED = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tas_unc_corr_local',
    component='Locally correlated uncertainty', quantity='average')

SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_SYSTEMATIC = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tas_unc_sys',
    component='Systematic uncertainty', quantity='average')

SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_CLOUD =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tas_unc_cloud',
    component='Cloud component of uncertainty', quantity='average')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_EXCLUDINGCLOUD = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_no_cloud',
    component='Total uncertainty excluding cloud', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_UNCORRELATED = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_rand',
    component='Random uncertainty', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_SYNOPTICALLYCORRELATED = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_corr_local',
    component='Locally correlated uncertainty', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_SYSTEMATIC = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_sys',
    component='Systematic uncertainty', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_CLOUD =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmin_unc_cloud',
    component='Cloud component of uncertainty', quantity='minimum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_EXCLUDINGCLOUD = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_no_cloud',
    component='Total uncertainty excluding cloud', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_UNCORRELATED = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_rand',
    component='Random uncertainty', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_SYNOPTICALLYCORRELATED = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_corr_local',
    component='Locally correlated uncertainty', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_SYSTEMATIC = OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_sys',
    component='Systematic uncertainty', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_CLOUD =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY,
    'tasmax_unc_cloud',
    component='Cloud component of uncertainty', quantity='maximum')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY = OutputVariableTemplate.extend(
    TEMPLATE_TEMPERATURE_UNCERTAINTY,
    template_long_name='{component} on average daily surface air temperature')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_RANDOM =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_rand',
    component='Random uncertainty')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_corr_sat',
    component='Locally correlated uncertainty (from satellite retrieval)')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_sys',
    component='Systematic uncertainty')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED2 =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_corr_mod',
    component='Locally correlated uncertainty (from surface-air model)')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC2 =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_sys_mod',
    component='Systematic uncertainty (from surface-air model)')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER0 =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_parameter_0',
    component='Systematic uncertainty mean offset')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER1 =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_parameter_1',
    component='Systematic uncertainty first fourier component')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER2 =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_parameter_2',
    component='Systematic uncertainty second fourier component')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER3 =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_parameter_3',
    component='Systematic uncertainty third fourier component')

SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER4 =  OutputVariable.from_template(
    SURFACEAIRMODEL_ANCILLARY_OCEAN_TEMPLATE_UNCERTAINTY,
    'tas_unc_parameter_4',
    component='Systematic uncertainty fourth fourier component')

OLD_OCEAN_ANCILLARY_NAMES = ['unc_rand_tas', 'unc_corr_tas', 'unc_syst_tas', 'unc_corr2_tas', 'unc_syst2_tas', 
                             'unc_parameter_0_tas', 'unc_parameter_1_tas', 'unc_parameter_2_tas', 'unc_parameter_3_tas', 'unc_parameter_4_tas']

STATION_NAME = OutputVariable(
    name='station_name',
    dtype=numpy.dtype('S1'),
    fill_value=numpy.int8(0),
    long_name='Station name',
    cf_role='timeseries_id')

TASMIN_QC = OutputVariable(
    name='tasmin_qc',
    dtype=numpy.uint8,
    fill_value=default_fillvals['u1'],
    standard_name='air_temperature status_flag',
    long_name='Quality control flags',
    valid_range=numpy.zeros((2,), numpy.uint8),
    flag_masks=numpy.empty((0,), numpy.uint8),
    flag_meanings='')

TASMAX_QC = OutputVariable(
    name='tasmax_qc',
    dtype=numpy.uint8,
    fill_value=default_fillvals['u1'],
    standard_name='air_temperature status_flag',
    long_name='Quality control flags',
    valid_range=numpy.zeros((2,), numpy.uint8),
    flag_masks=numpy.empty((0,), numpy.uint8),
    flag_meanings='')

DETECTION_TIME = OutputVariable(
    name='detection_time',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    standard_name='time',
    long_name='Start time of period for break detection status report (days)',
    units=TIME_UNITS_DAYS_SINCE_EPOCH)

TASMIN_DETECTION_QC = OutputVariable(
    name='tasmin_detection_qc',
    dtype=numpy.uint8,
    fill_value=default_fillvals['u1'],
    long_name='Quality control flags for break detection',
    standard_name='air_temperature status_flag',
    valid_range=numpy.array([ 0, 2 ], numpy.uint8),
    flag_values=numpy.array([ 0, 1, 2 ], numpy.uint8),
    flag_meanings='not_possible possible_but_unreliable reliable')

TASMAX_DETECTION_QC = OutputVariable(
    name='tasmax_detection_qc',
    dtype=numpy.uint8,
    fill_value=default_fillvals['u1'],
    long_name='Quality control flags for break detection',
    standard_name='air_temperature status_flag',
    valid_range=numpy.array([ 0, 2 ], numpy.uint8),
    flag_values=numpy.array([ 0, 1, 2 ], numpy.uint8),
    flag_meanings='not_possible possible_but_unreliable reliable')

TASMIN_BREAK_STATION = OutputVariable(
    name='tasmin_break_station',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'],
    long_name='Index of station at which break was detected')

TASMAX_BREAK_STATION = OutputVariable(
    name='tasmax_break_station',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'],
    long_name='Index of station at which break was detected')

TASMIN_BREAK_AMPLITUDE = OutputVariable(
    name='tasmin_break_amplitude',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    long_name='Inhomogeneity in minimum surface air temperature',
    units='K')

TASMAX_BREAK_AMPLITUDE = OutputVariable(
    name='tasmax_break_amplitude',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    long_name='Inhomogeneity in maximum surface air temperature',
    units='K')

TASMIN_BREAK_TIME_BOUNDS = OutputVariable(
    name='tasmin_break_time_bounds',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    long_name='Relative time bounds of break evaluation period (days)',
    units='days')

TASMAX_BREAK_TIME_BOUNDS = OutputVariable(
    name='tasmax_break_time_bounds',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    long_name='Relative time bounds of break evaluation period (days)',
    units='days')

TASMIN_BREAK_TIME_AFFECTED_BOUNDS = OutputVariable(
    name='tasmin_break_time_affected_bounds',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    long_name='Relative time bounds of period affected by break (days)',
    units=TIME_UNITS_DAYS_SINCE_EPOCH)

TASMAX_BREAK_TIME_AFFECTED_BOUNDS = OutputVariable(
    name='tasmax_break_time_affected_bounds',
    dtype=numpy.float32,
    fill_value=default_fillvals['f4'],
    long_name='Relative time bounds of period affected by break (days)',
    units=TIME_UNITS_DAYS_SINCE_EPOCH)

MERGED_BREAK_TIME = OutputVariable(
    name='merged_break_time',
    long_name='Breakpoint time (days)',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'],
    units=TIME_UNITS_DAYS_SINCE_EPOCH)

MERGED_BREAK_STATION = OutputVariable(
    name='merged_break_station',
    long_name='Index of station at which break was detected',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'])

MERGED_BREAK_LIKELIHOOD = OutputVariable(
    name='merged_break_likelihood',
    long_name='Likelihood of breakpoint',
    dtype=numpy.int32,
    fill_value=-128,
    valid_range=numpy.array([ 1, 127 ], numpy.uint8))

DETECTION_FEASIBILITY = OutputVariable(
    name='detection_feasibility',
    long_name='Feasibility of detection flag',
    dtype=numpy.int8,
    fill_value=-128,
    valid_range=numpy.array([ 0, 2 ], numpy.uint8),
    flag_values=numpy.array([ 0, 1, 2 ], numpy.uint8),
    flag_meanings = 'not_possible only_absolute_possible all_tests_possible')

TASMIN_BREAK_TIME = OutputVariable(
    name='tasmin_break_time',
    long_name='Breakpoint time (days)',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'],
    units=TIME_UNITS_DAYS_SINCE_EPOCH)

TASMIN_BREAK_STATION = OutputVariable(
    name='tasmin_break_station',
    long_name='Index of station at which break was detected',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'])

TASMIN_BREAK_LIKELIHOOD = OutputVariable(
    name='tasmin_break_likelihood',
    long_name='Likelihood of breakpoint',
    dtype=numpy.int32,
    fill_value=-128,
    valid_range=numpy.array([ 1, 127 ], numpy.uint8))

TASMIN_BREAK_TYPE = OutputVariable(
    name='tasmin_break_type',
    long_name='Type of breakpoint',
    dtype=numpy.int8,
    fill_value=-128,
    valid_range=numpy.array([ 0, 2 ], numpy.uint8),
    flag_values=numpy.array([ 0, 1, 2 ], numpy.uint8),
    flag_meanings='relative_test absolute_test merging_metadata')

TASMIN_DETECTION_SCORE = OutputVariable(
    name='tasmin_detection_score',
    long_name='Score for break detection',
    dtype=numpy.int8,
    fill_value=-128,
    valid_range=numpy.array([ 0, 8 ], numpy.uint8))

TAS_BREAK_TIME = OutputVariable(
    name='tas_break_time',
    long_name='Breakpoint time (days)',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'],
    units=TIME_UNITS_DAYS_SINCE_EPOCH)

TAS_BREAK_STATION = OutputVariable(
    name='tas_break_station',
    long_name='Index of station at which break was detected',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'])

TAS_BREAK_LIKELIHOOD = OutputVariable(
    name='tas_break_likelihood',
    long_name='Likelihood of breakpoint',
    dtype=numpy.int32,
    fill_value=-128,
    valid_range=numpy.array([ 1, 127 ], numpy.uint8))

TAS_BREAK_TYPE = OutputVariable(
    name='tas_break_type',
    long_name='Type of breakpoint',
    dtype=numpy.int8,
    fill_value=-128,
    valid_range=numpy.array([ 0, 2 ], numpy.uint8),
    flag_values=numpy.array([ 0, 1, 2 ], numpy.uint8),
    flag_meanings='relative_test absolute_test merging_metadata')

TAS_DETECTION_SCORE = OutputVariable(
    name='tas_detection_score',
    long_name='Score for break detection',
    dtype=numpy.int8,
    fill_value=-128,
    valid_range=numpy.array([ 0, 8 ], numpy.uint8))

TASMAX_BREAK_TIME = OutputVariable(
    name='tasmax_break_time',
    long_name='Breakpoint time (days)',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'],
    units=TIME_UNITS_DAYS_SINCE_EPOCH)

TASMAX_BREAK_STATION = OutputVariable(
    name='tasmax_break_station',
    long_name='Index of station at which break was detected',
    dtype=numpy.int32,
    fill_value=default_fillvals['i4'])

TASMAX_BREAK_LIKELIHOOD = OutputVariable(
    name='tasmax_break_likelihood',
    long_name='Likelihood of breakpoint',
    dtype=numpy.int32,
    fill_value=-128,
    valid_range=numpy.array([ 1, 127 ], numpy.uint8))

TASMAX_BREAK_TYPE = OutputVariable(
    name='tasmax_break_type',
    long_name='Type of breakpoint',
    dtype=numpy.int8,
    fill_value=-128,
    valid_range=numpy.array([ 0, 2 ], numpy.uint8),
    flag_values=numpy.array([ 0, 1, 2 ], numpy.uint8),
    flag_meanings='relative_test absolute_test merging_metadata')

TASMAX_DETECTION_SCORE = OutputVariable(
    name='tasmax_detection_score',
    long_name='Score for break detection',
    dtype=numpy.int8,
    fill_value=-128,
    valid_range=numpy.array([ 0, 8 ], numpy.uint8))

BIAS_SOURCES_PATTERN_DICTIONARY = {'insitu_land':'insitu_land_%A_%Y%m%d.bin',
			   'insitu_ocean':'insitu_ocean_%A_%Y%m%d.bin',
			   'surfaceairmodel_land_global':'surfaceairmodel_land_%A_%Y%m%d.bin',
			   'surfaceairmodel_ocean_global':'surfaceairmodel_ocean_%A_%Y%m%d.bin',
			   'surfaceairmodel_ice_global':'surfaceairmodel_ice_%A_%Y%m%d.bin'}
