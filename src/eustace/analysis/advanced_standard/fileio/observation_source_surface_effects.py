"""
Collection of methods describing covariate effect for a given source
"""

import numpy

BREAK_STATION_INDEX_SHIFT = 1       # Adjustment for indexing from 1 in input breakpoint file

def insitu_land_covariate_effect(time, station_indices, observed_breakpoints):
    """Define covariate effect related to temperature series breakpoints. It assigns a specific bias term for each breaking point period detected, sparing only the most recent
    interval of measurements"""
    
    
    # Consider the following arrays, labeling breakpoints and stations that have detected such breakpoints
    #
    #   B = [t1(s1), t2(s1), ..., t1(si), t2(si), t3(si), ..., tk(sn)]
    #   S = [   s1,     s1,  ...,    si,     si,     si, ...,     sn]
    #
    # Apply the following transformation
    #
    #   B -> B' = [0, t1(s1), t2(s1), ..., 0, t1(si), t2(si), t3(si), ..., 0, tk(sn)]
    #
    # Now compute B'-t*, where t* is the date time (expressed as days elapsed from 1850/01/01)
    #
    #   B' - t* = [-t*, t1(s1)-t*, t2(s1)-t*, ..., -t*, t1(si)-t*, t2(si)-t*, t3(si)-t*, ..., -t*, tk(sn)-t*]
    #
    #   (B' - t*) elements can be either positive or negative, e.g.:
    #
    #   sign(B' - t*) = [-, -, +, ..., -, +, -, -, ..., -, +]
    #
    # To detected the right break point interval, look for -/+ ("minus to plus") transitions. In the example above, they would be
    #
    #   B_valid = t2(s1), t1(si),tk(sn)
    #   S_valid =    s1,     si,    sn
    #
    # Then, for t=t*, you would need to associate a bias term to station s1, si and sn. The associated bias term for each station is given by the index I representing the position of
    # tj(si) inside the array B, e.g. I(t2(s1) = 1)
    #
    #   effect = [[s1, I(t2(s1)], [si, I(t1(si)], [sn, I(tk(sn)]]
    #
    # Now we just need to remap the station index into the index of available observations, which is determined by using the array of valid indices
    #
    #   observations =  [O[s0], O[s1], O[s2], ..., O[sn]] ,          valid_indices = [i1, i2, ..., ik]
    #   observations[valid_indices] = [O[i1], O[i2], ..., O[ik]]
    #
    #   valid_effect = [[s1(valid_indices), I(t2(s1)], [si(valid_indices), I(t1(si)], [sn(valid_indices), I(tk(sn)]]
    #
    #   where si(valid_indices) is the mapping from the value of the station index to the index position inside valid_indices, e.g.
    #
    #   observations = [0[0], O[1], ...,O[n]],       valid_indices = [1, 3, 4, 7, n-1]
    #   observations[valid_indices] = [O[1], O[3], O[4], O[7], O[n-1]]
    #   effect = [[1, 4], [2, 6], [7, 1]]
    #   valid_effect = [[0,4], [3, 1]]
    #
    
    if isinstance(time, numpy.ndarray):
        time = time[0]
    
    # station indices shifted to the analysis indexing that starts at zero
    break_stations = observed_breakpoints.break_station[:] - BREAK_STATION_INDEX_SHIFT
    new_adjustment_times    = observed_breakpoints.break_time[:]
    
    # subset to only those adjutments with adjustment times later than the current time
    time_differences = new_adjustment_times - time
    
    adjustment_indices = numpy.arange(len(break_stations))[time_differences > 0]
    break_stations     = break_stations[time_differences > 0]
    time_differences   = time_differences[time_differences > 0]

    # The most recent adjustment period is that corresponding to a change of station number in the subsetted break_stations
    most_recent_adjustment = break_stations!=numpy.roll(break_stations,1)
    if len( most_recent_adjustment ) > 0 and numpy.sum( most_recent_adjustment ) == 0:
        # only one station selected so find the most recent adjustment as the smallest time difference
        # rather than from changes in station indexes
        most_recent_adjustment = [numpy.argmin(time_differences),]
    
    # extract the indices to the most recent active adjustments/state variables and the corresponding stations indices
    active_adjustment_indices = adjustment_indices[most_recent_adjustment]
    corresponding_stations_indices     = break_stations[most_recent_adjustment]

    # map affected station inds in corresponding_stations_indices to the input station_indices
    station_restriction = numpy.nonzero( numpy.isin(station_indices, corresponding_stations_indices) )[0]
    break_restriction = numpy.nonzero( numpy.isin(corresponding_stations_indices, station_indices) )[0]
    
    # construct bias effect
    effect = numpy.column_stack( ( station_restriction, active_adjustment_indices[break_restriction] ) )
    
    if(len(effect)):
      message = 'INSITU LAND: {0} bias terms for "insitu_land" group added for date {1}'.format(len(effect),time)
      print(message)

      return effect
    else:
      return None
      
def global_satellite_effect(valid_indices):
    """
    Define covariate effects related to global biases for a given satellite source.
    """
    return numpy.column_stack((numpy.arange(len(valid_indices)), numpy.repeat(0,len(valid_indices))))
    
def hemispheric_satellite_effect(valid_indices, valid_latitudes):
    """
    Define covariate effects related to hemispheric biases for a given satellite source.
    """
    return numpy.column_stack((numpy.arange(len(valid_indices)), valid_latitudes < 0.0 ))