"""Describe input observation sources to analyse."""

import numpy
from datetime import datetime
from eumopps.timeutils.timebase import TimeBaseDays

class Observations(object):
    """References to NumPy arrays describing an observation source."""

    def __init__(self, mask, time, location, measurement, uncorrelatederror, locallycorrelatederror):
        """Build observation set from NumPy arrays."""
        self.mask = mask
        self.time = time
        self.location = location
        self.measurement = measurement
        self.uncorrelatederror = uncorrelatederror
        self.locallycorrelatederror = locallycorrelatederror

    def number_of_observations(self):
        """Total number of observations."""

        return self.measurement.shape[0]

    @staticmethod
    def mean(a, b):
        """Compute mean of two observation sets."""

        # Use locations in common
        a_restriction = numpy.nonzero( numpy.isin(a.location, b.location) )[0]
        b_restriction = numpy.nonzero( numpy.isin(b.location, a.location) )[0]

        # locations should be identical
        # (need this check because code below assumes original locations were in same order)
        if not numpy.all(a.location[a_restriction] == b.location[b_restriction]):
            raise ValueError('Mismatch in location identifiers.')
        common_locations = a.location[a_restriction]

        # Mask as invalid if a OR b are invalid
        mask = numpy.logical_or(a.mask[a_restriction], b.mask[b_restriction])

        # times should be identical
        if not numpy.all(a.time == b.time):
            raise ValueError('Mismatch in times.')
        common_time = a.time
        
        # Mean measurements
        measurement = ((numpy.ma.masked_array(a.measurement[a_restriction], mask) + numpy.ma.masked_array(b.measurement[b_restriction], mask)) * numpy.float32(0.5)).data

        # RMS uncorrelated error
        uncmax = numpy.ma.masked_array(a.uncorrelatederror[a_restriction], mask)
        uncmin = numpy.ma.masked_array(b.uncorrelatederror[b_restriction], mask)
        uncorrelatederror = (numpy.ma.sqrt((uncmax*uncmax + uncmin*uncmin) * numpy.float32(0.5))).data

        # RMS correlated error components (must have same length scale)
        if len(a.locallycorrelatederror) != len(b.locallycorrelatederror):
            raise ValueError('Mismatch in number of locally correlated error components.')
        locallycorrelatederror = [ ]
        for index in range(len(a.locallycorrelatederror)):
            corrmax = numpy.ma.masked_array((a.locallycorrelatederror[index])[a_restriction], mask)
            corrmin = numpy.ma.masked_array((b.locallycorrelatederror[index])[b_restriction], mask)
            locallycorrelatederror.append( (numpy.ma.sqrt((corrmax*corrmax + corrmin*corrmin) * numpy.float32(0.5))).data )
        
        # done
        return Observations(mask, common_time, common_locations, measurement, uncorrelatederror, locallycorrelatederror)

    def append(self, other, locationoffset):

        self.mask = numpy.hstack((self.mask.ravel(), other.mask.ravel()))
        if self.time != other.time:
            raise ValueError('Mismatch in time')
        self.location = numpy.hstack((self.location, (other.location+locationoffset)))
        self.measurement = numpy.hstack((self.measurement.ravel(), other.measurement.ravel()))
        self.uncorrelatederror = numpy.hstack((self.uncorrelatederror.ravel(), other.uncorrelatederror.ravel()))
        if len(self.locallycorrelatederror) != len(other.locallycorrelatederror):
            raise ValueError('Mismatch in number of locally correlated error components.')
        for index in range(len(self.locallycorrelatederror)):
            self.locallycorrelatederror[index] = numpy.hstack((self.locallycorrelatederror[index].ravel(), other.locallycorrelatederror[index].ravel()))

class ObservationsBreakPoints(object):
    """References to NumPy arrays describing an observation break point source."""

    MAX_LIKELIHOOD = 127

    def __init__(self, break_time, break_station, break_likelihood, break_type = None, detection_feasibility = None, detection_score = None):
        """Build observation set from NumPy arrays."""
        self.break_time = break_time
	self.break_station = break_station
	self.break_likelihood = break_likelihood
	self.POLICY_CODES = ['HARD_CUTOFF', 'LAPLACE_KERNEL']

	# optional, depending on the observable
	self.break_type = break_type
	self.detection_feasibility = detection_feasibility
	self.detection_score =  detection_score

    def number_of_observations(self):
        """Total number of observations."""

        return self.break_time.shape[0]

    def total_number_of_stations(self):
	"""Total number of stations"""

	return self.stations_collection().size

    def stations_collection(self):
	"""Stations that detected breaking points"""

	return self.station_statistics()[0]

    def stations_count(self):
	"""Number of occurrencies for each station"""

	return self.station_statistics()[1]

    def station_statistics(self):
	"""Station indices"""

	return numpy.unique(self.break_station, return_counts = True)

    def filtered_break_points(self, station_index):
      """Breaking points values filtered down to a specific station index"""

      return self.break_time[self.stations_break_point_indices(station_index)]

    def stations_break_point_indices(self, station_index):
      """Indices of the first occurrences of the station values in the original array"""
  
      return numpy.nonzero(numpy.isin(self.break_station, station_index))[0]

    def break_points_timestamp(self, epoch):
	"""Transform breaking point value into a date of the from YYYYMMDD """
	converter = TimeBaseDays(epoch)
	return numpy.array(['{0.year:4d}-{0.month:02d}-{0.day:02d}'.format(converter.number_to_datetime(int(number))) for number in self.break_time])#, dtype='datetime64[D]')
	
    def apply_policy(self, policy_code, threshold = 0., decay_constant = 0.):
	"""Reject temperature series breaking points according to a specified policy.
	  Possible poplicy values are:
	  1) BASE: only breaking points with detection feasibility flags equal to 0 gets removed.
	  2) COMBINED: extension of the BASE policy, where additional removals can be applied according to the \'likelihood\' value of the breaking point"""

	if policy_code not in self.POLICY_CODES:
	    raise ValueError('Wrong value for \"policy_code\" parameter: possible values are '+str(self.POLICY_CODES))
	elif policy_code == self.POLICY_CODES[0]:
	    self.apply_base_policy(threshold)
	elif policy_code == self.POLICY_CODES[1]:
	    self.apply_kernel_policy(threshold, decay_constant)

    def apply_base_policy(self, threshold):
	"""Filter information on breaking points by removing observations with likelihood smaller than a given threshold"""

	indices_to_keep = numpy.argwhere(self.break_likelihood >= threshold)
	indices_to_keep = indices_to_keep.reshape(indices_to_keep.shape[0],)
        self.break_time = self.break_time[indices_to_keep]
	self.break_station = self.break_station[indices_to_keep]
	self.break_likelihood = self.break_likelihood[indices_to_keep]

    def apply_kernel_policy(self, threshold, decay_constant):
	"""Filter information on breaking points by using a Gaussian smoothing kernel: the kernel depends on the distance between the maximal likelihood value (127), and the actual likelihood value.
           The rejection is then applied by filtering out those points with kernel smaller than a given threshold"""

	gaussian_kernels = numpy.exp(-decay_constant*numpy.abs(self.break_likelihood-127))
	indices_to_keep = numpy.argwhere(gaussian_kernels >= threshold)
	indices_to_keep = indices_to_keep.reshape(indices_to_keep.shape[0],)
        self.break_time = self.break_time[indices_to_keep]
	self.break_station = self.break_station[indices_to_keep]
	self.break_likelihood = self.break_likelihood[indices_to_keep]

class ObservationSource(object):
    """Define an observation source."""

    TMEAN = 'Tmean'
    TMIN = 'Tmin'
    TMAX = 'Tmax'

    def __init__(self):
        """Empty constructor in bass class"""
        pass

    def observables(self):
        """The names of variables estimated from this source."""
        raise NotImplementedError

    def observation_location_lookup(self):
        """NumPy array in which column number corresponds to location id and rows are latitude and longitude."""
        raise NotImplementedError

    def observations(self, observable):
        """Retrieve observations for specified observable quantity."""
        raise NotImplementedError

    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty."""
        raise NotImplementedError
