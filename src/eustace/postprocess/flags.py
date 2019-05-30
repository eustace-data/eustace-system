""""

Flags are computed in the following modules:

threshold_check:   CLIMATOLOGY_UNC_FLAG, LARGE_SCALE_UNC_FLAG, DAY_FLAG
window_checks:     PRIOR_WINDOW_FLAG, POST_WINDOW_FLAG
calendar_checks:   CALENDAR_DAY_FLAG, PRIOR_CALENDAR_FLAG, POST_CALENDAR_FLAG
outlier_check:     AREAL_LOW_FLAG, AREAL_HIGH_FLAG, EXTREME_LOW_FLAG, EXTREME_HIGH_FLAG

"""

import numpy

# bit mask flags used in output file
FLAG_TYPE = numpy.uint16
TYPE_NAME = 'u2'
FLAG_MAX_N = 16
FLAG_MAX_USED = 13
NULL_FLAG = FLAG_TYPE(0)


# time window checks
DAY_FLAG             = FLAG_TYPE(1) << 0 # not observed today (observation inluence derived)
PRIOR_WINDOW_FLAG    = FLAG_TYPE(1) << 1 # not observed is last few days (observation inluence derived)
POST_WINDOW_FLAG     = FLAG_TYPE(1) << 2 # not observed in next few days (observation inluence derived)
# calendar day checks
PRIOR_CALENDAR_FLAG  = FLAG_TYPE(1) << 3 # not observed in previous years (all of DAY_FLAG, PRIOR_WINDOW_FLAG and POST_WINDOW_FLAG set for all previous years)
POST_CALENDAR_FLAG   = FLAG_TYPE(1) << 4 # not observed in later years (all of DAY_FLAG, PRIOR_WINDOW_FLAG and POST_WINDOW_FLAG set for all previous years)
# slow component uncertatinty thresholds
CLIMATOLOGY_UNC_FLAG = FLAG_TYPE(1) << 5 # climatology component not constrained (component uncertainty based)
LARGE_SCALE_UNC_FLAG = FLAG_TYPE(1) << 6 # large scale component not constrained (component uncertainty based)
# extreme value checks
AREAL_LOW_FLAG       = FLAG_TYPE(1) << 7 # inconsitent cold in region
AREAL_HIGH_FLAG      = FLAG_TYPE(1) << 8 # inconsitent warmth in region
EXTREME_LOW_FLAG     = FLAG_TYPE(1) << 9 # unrealistically hot exceeding fixed maximum temperature threshold
EXTREME_HIGH_FLAG    = FLAG_TYPE(1) << 10 # unrealistically cold exceeding fixed minimum temperature threshold
PERSISTANT_OUTLIER_FLAG = FLAG_TYPE(1) << 11 # persistantly failed any combination of AREAL_LOW_FLAG, AREAL_HIGH_FLAG, EXTREME_LOW_FLAG or EXTREME_HIGH_FLAG in time window

# omitted data source flags
MISSING_MARINE_FLAG  = FLAG_TYPE(1) << 12 # marine data source unavailable

#CLIMATOLOGY_LATENT_FLAG = FLAG_TYPE(1) << 13 # poorly constrained variables contribute to this locations value
#LARGE_SCALE_LATENT_FLAG = FLAG_TYPE(1) << 14 # poorly constrained variables contribute to this locations value

# location based percentile exceedences
LOCATION_LOW_FLAG     = FLAG_TYPE(1) << 13 # unrealistically hot for location
LOCATION_HIGH_FLAG    = FLAG_TYPE(1) << 14 # unrealistically cold for location

# missing data indicator
MISSING_FLAG_FLAG    = FLAG_TYPE(1) << 15 # calendar day not constrained on a sufficient number of years

# useful derived combinations of elemental flags
CALENDAR_DAY_FLAG = DAY_FLAG | PRIOR_WINDOW_FLAG | POST_WINDOW_FLAG # helper combination that says that there is no observation nearby in time

FLAG_MEANINGS = ['day_constraint',
                 'prior_window_constraint',
                 'post_window_constraint',
                 'prior_calendar_constraint',
                 'post_calendar_constraint',
                 'climatology_uncertainty_qc',
                 'large_scale_uncertainty_qc',
                 'areal_low_value_qc',
                 'areal_high_value_qc',
                 'extreme_low_value_qc',
                 'extreme_high_value_qc',
                 'marine_location_flag',
                 'location_low_value_qc',
                 'location_high_value_qc',
                 'unconstrained_calendar_day_qc',
                 ]

