"""
In-situ land source
-------------------
Observation source interface to output of in-situ land preprocessing.
"""

import numpy
import netCDF4
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from datetime import datetime
from eustace.timeutils.epoch import days_since_epoch
from eustace.outputformats import definitions

ALTITUDE_SCALE_FACTOR = 0.001 # conversion scale factor from m to km

class ObservationSourceInsituLand(ObservationSource):
    """Provide observation source interface to output of pre-processed in-situ land data.
       Optionally restrict to a single day."""

    # Map observation source keys to file variable definitions
    OBSERVATIONMAPS = { ObservationSource.TMIN: definitions.TASMIN, ObservationSource.TMAX: definitions.TASMAX }
    QCOBSERVATIONMAPS = { ObservationSource.TMIN: definitions.TASMIN_QC, ObservationSource.TMAX: definitions.TASMAX_QC }

   
    def __init__(self, filename, daynumber=None, hold_out_filename=None, altitude_adjustment = None):
        """Load and parse specified filename according to specified format instance."""

        super(ObservationSource, self).__init__()
        self.dataset = netCDF4.Dataset(filename, 'r')
        self.daynumber = daynumber
        
        self.altitude_adjustment = altitude_adjustment
        
        if hold_out_filename != None:
            self.hold_out_dataset = netCDF4.Dataset(hold_out_filename, 'r')
        else:
            self.hold_out_dataset = False
            
    def identify_validation_stations(self):
	"""Determine those indices referring to stations that will be removed from the training set"""
	
	variable = 'station_code'
	
	size_dataset = self.dataset.variables[variable].shape[1]
	view_dataset = '|S'+str(size_dataset)
	available_stations = self.dataset.variables[variable][:].data.view(view_dataset).ravel()

	size_holdout = self.hold_out_dataset.variables[variable].shape[1]
	view_holdout = '|S'+str(size_holdout)
	filtered_stations = self.hold_out_dataset.variables[variable][:].data.view(view_holdout).ravel()

	filtering_mask = numpy.in1d(available_stations, filtered_stations).reshape(-1,1).T
	
	return filtering_mask
	
    def observables(self):
        """The names of variables estimated from this source."""
        return ObservationSourceInsituLand.OBSERVATIONMAPS.keys()

    def number_of_days(self):
        """Total number of days."""

        num_times = len( self.dataset.dimensions[definitions.DIMENSION_NAME_TIME] )
        if self.daynumber is None:
            return num_times
        else:
            return 1

    def number_of_stations(self):
        """Total number of stations."""

        num_stations = len( self.dataset.dimensions[definitions.DIMENSION_NAME_STATION] )
        return num_stations

    def number_of_observations(self):
        """Total number of observations."""

        num_times = len( self.dataset.dimensions[definitions.DIMENSION_NAME_TIME] )
        num_stations = len( self.dataset.dimensions[definitions.DIMENSION_NAME_STATION] )
        if self.daynumber is None:
            return num_times * num_stations
        else:
            return num_stations

    def observation_location_lookup(self):
        """NumPy array of 3D coordinates of latitude, longitude, reference time."""

        # station coordinates
        longitude = self.dataset.variables[definitions.LONGITUDE.name][:]
        latitude = self.dataset.variables[definitions.LATITUDE.name][:]

        # allow for longitudes expressed as [0,360) instead of [-180,180)
        wrapped_longitudes = numpy.nonzero(longitude >= numpy.float(180.0))
        longitude[wrapped_longitudes] -= 360.0

        # combined coordinate
        location = numpy.vstack( [ latitude, longitude ] )

        return location

    def observations(self, observable, qc_mask):
        """Retrieve observations for specified observable (tmin, tmax) quantity.
           Optionally filter out observation based on elements that fail qc checks
        """

        # Retrieve measurement and variable info
        variable = self.dataset.variables[ ObservationSourceInsituLand.OBSERVATIONMAPS[observable].name ]
        qc_var = self.dataset.variables[ ObservationSourceInsituLand.QCOBSERVATIONMAPS[observable].name ]
       
        variable.set_auto_maskandscale(True)
        time = self.dataset.variables[definitions.TIME.name][:]
        
        # Restrict by time if requested
        if self.daynumber is None:
            measurement = variable[:]
        
        else:
            # argwhere - Finds the indices of array elements that are non-zero, grouped by element
            # returns indices of the time array where current day number is specified
            dayindex = numpy.argwhere(time == self.daynumber)
            if len(dayindex) == 1:
                # numpy.ravel returns a contiguous flattened array.
                # extract measurement by advanced index and slice  
                measurement = variable[dayindex.ravel(), :]
                time = time[dayindex]

                # the mask is updated so True areas also represent areas where 
                # qc flags have failed checks  
                if qc_mask:
                    # replicate method used for measurement to return data relating to specified day indices 
                    qc_var_day = qc_var[dayindex.ravel(), :]
                    # create a mask for areas where the qc flag values are not 0
                    # this mask - True values - relates to all data that failed qc checks
                    qc_data_mask = numpy.ma.masked_not_equal(qc_var_day, 0)
                    # update the measurement mask to additionally mask out areas from the qc data mask, 
                    # in addition to those areas already masked out
                    measurement.mask = numpy.logical_or(
                                       measurement.mask, qc_data_mask)

            else:
                return None                        

        if self.hold_out_dataset:
            # mask out those stations not used for training the learning algorithm
            measurement.mask = numpy.logical_or(measurement.mask, self.identify_validation_stations())

            # the mask is updated so True areas also represent areas where 
            # qc flags have failed checks  
            if qc_mask:
                qc_data_mask = numpy.ma.masked_not_equal(qc_var, 0)
                measurement.mask = numpy.logical_or(
                                   measurement.mask, qc_data_mask)
                                   
        # Formulate output
        all_time = numpy.tile(time, (measurement.shape[1], 1)).transpose()
        all_stations = numpy.tile(numpy.array(range(measurement.shape[1]), numpy.uint64), (measurement.shape[0], 1))
        uncorrelatederror = numpy.empty((measurement.shape), numpy.float32)
        uncorrelatederror.fill(numpy.float32(0.3))
        locallycorrelatederror = [ ]

        # Apply an adjustment for altitude if an AltitudeAdjustment object has been provided
        if self.altitude_adjustment is not None:
            
            locations  = self.observation_location_lookup().T
            elevations = self.dataset.variables[definitions.ELEVATION.name][:] * ALTITUDE_SCALE_FACTOR
            
            self.altitude_adjustment.load_reference_map()
            adjustment = self.altitude_adjustment.adjustment(elevations, locations)
            adjustment_uncertainty = self.altitude_adjustment.adjustment_uncertainty(elevations, locations)

            measurement[~measurement.mask] = (adjustment + measurement)[~measurement.mask]
            uncorrelatederror[~measurement.mask] = numpy.sqrt(adjustment_uncertainty**2 + uncorrelatederror**2).astype(numpy.float32)[~measurement.mask]
            
            if hasattr(elevations, 'mask'):
                measurement.mask = numpy.logical_or(measurement.mask, elevations.mask)
            

        return Observations(
            measurement.mask.ravel(),
            all_time.ravel(),
            all_stations.ravel(),
            measurement.data.ravel().astype(numpy.float32),
            uncorrelatederror.ravel(),
            locallycorrelatederror)

    def local_correlation_length_scale(self, observable):
        """Length scale for locally correlated component of uncertainty (empty array in this case)."""
        
        return numpy.array([ ], numpy.float32)
