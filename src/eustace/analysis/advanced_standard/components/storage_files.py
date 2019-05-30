"""Storage to file."""

import numpy
import cPickle
from spacetime import SpaceTimeComponentSolutionStorage
from spatial import SpatialComponentSolutionStorage
from eustace.outputformats.ensuredirectory import ensuredirectory

import os

def writepickle(quantity, filename):
    ensuredirectory(filename)
    with open(filename, 'wb') as fh:
        cPickle.dump(quantity, fh, protocol=2)

def readpickle(filename):
    with open(filename, 'rb') as fh:
        return cPickle.load(fh)

class SpaceTimeComponentSolutionStorage_Files(SpaceTimeComponentSolutionStorage):

    def __init__(self, 
                 statefilename_read=None,
                 marginal_std_filename_read = None,
                 sample_filename_read = None,
                 prior_sample_filename_read = None,
                 statefilename_write=None,
                 marginal_std_filename_write = None,
                 sample_filename_write = None,
                 prior_sample_filename_write = None,
                 measurement_time_index_write=None,
                 measurementfilename_write=None,
                 measurementfilelist_read=None):

        super(SpaceTimeComponentSolutionStorage_Files, self).__init__()
        
        self.statefilename_read = statefilename_read
        self.marginal_std_filename_read = marginal_std_filename_read
        self.sample_filename_read = sample_filename_read
        self.prior_sample_filename_read = prior_sample_filename_read
        
        self.statefilename_write = statefilename_write
        self.marginal_std_filename_write = marginal_std_filename_write
        self.sample_filename_write = sample_filename_write
        self.prior_sample_filename_write = prior_sample_filename_write
        self.measurement_time_index_write = measurement_time_index_write
        
        self.measurementfilename_write = measurementfilename_write
        self.measurementfilelist_read = measurementfilelist_read

    def state_write(self, state):

        if self.statefilename_write is None:
            raise ValueError('no state write file specified')

        writepickle(state, self.statefilename_write)

    def state_marginal_std_write(self, marginal_std):

        if self.marginal_std_filename_write is None:
            raise ValueError('no marginal variances write file specified')

        writepickle(marginal_std, self.marginal_std_filename_write)

    def state_sample_write(self, sample):

        if self.sample_filename_write is None:
            raise ValueError('no sample write file specified')

        writepickle(sample, self.sample_filename_write)
        
    def state_prior_sample_write(self, sample):

        if self.prior_sample_filename_write is None:
            raise ValueError('no prior sample write file specified')

        writepickle(sample, self.prior_sample_filename_write)

    def state_read(self):

        if self.statefilename_read is None:
            return None

        return readpickle(self.statefilename_read)

    def state_marginal_std_read(self):

        if self.marginal_std_filename_read is None:
            return None

        return readpickle(self.marginal_std_filename_read)

    def state_sample_read(self):

        if self.sample_filename_read is None:
            return None

        return readpickle(self.sample_filename_read)

    def state_prior_sample_read(self):

        if self.prior_sample_filename_read is None:
            return None

        return readpickle(self.prior_sample_filename_read)

    def measurement_write(self, measurement, time_index):

        if self.measurementfilename_write is None:
            raise ValueError('no measurement write file specified')

        if self.measurement_time_index_write != time_index:
            raise ValueError('time index mismatch')

        writepickle(measurement, self.measurementfilename_write)

    def measurement_read(self, time_index):

        if self.measurementfilelist_read is None:
            raise ValueError('no measurement read list')

        if (time_index < 0) or (time_index >= len(self.measurementfilelist_read)):
            raise ValueError('time index of range of measurement read list')

        return readpickle(self.measurementfilelist_read[time_index])

    def timeindices_read(self):

        return range(len(self.measurementfilelist_read)) if self.measurementfilelist_read else None

class SpatialComponentSolutionStorage_Files(SpatialComponentSolutionStorage):

    def __init__(self, time_index=None, 
                       statefilename_read=None, 
                       marginal_std_filename_read = None, 
                       sample_filename_read = None,
                       prior_sample_filename_read = None,
                       statefilename_write=None, 
                       marginal_std_filename_write = None,
                       sample_filename_write = None,
                       prior_sample_filename_write = None):

        super(SpatialComponentSolutionStorage_Files, self).__init__()
        
        self.time_index = time_index
        
        self.statefilename_read = statefilename_read
        self.marginal_std_filename_read = marginal_std_filename_read
        self.sample_filename_read = sample_filename_read
        self.prior_sample_filename_read = prior_sample_filename_read
        
        self.statefilename_write = statefilename_write
        self.marginal_std_filename_write = marginal_std_filename_write
        self.sample_filename_write = sample_filename_write
        self.prior_sample_filename_write = prior_sample_filename_write
        
    def partial_state_write(self, state, time_index):

        if self.statefilename_write is None:
            raise ValueError('no state write file specified')

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        writepickle(state, self.statefilename_write)
                
    def partial_state_marginal_std_write(self, marginal_std, time_index):

        if self.marginal_std_filename_write is None:
            raise ValueError('no marginal variances write file specified')

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        writepickle(marginal_std, self.marginal_std_filename_write)

    def partial_state_sample_write(self, sample, time_index):

        if self.sample_filename_write is None:
            raise ValueError('no marginal variances write file specified')

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        writepickle(sample, self.sample_filename_write)

    def partial_state_prior_sample_write(self, sample, time_index):

        if self.prior_sample_filename_write is None:
            raise ValueError('no marginal variances write file specified')

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        writepickle(sample, self.prior_sample_filename_write)

    def partial_state_read(self, time_index):
        if self.statefilename_read is None:
            return None

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        return readpickle(self.statefilename_read)

    def partial_state_marginal_std_read(self, time_index):

        if self.marginal_std_filename_read is None:
            return None

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        return readpickle(self.marginal_std_filename_read)
    
    def partial_state_sample_read(self, time_index):

        if self.sample_filename_read is None:
            return None

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        return readpickle(self.sample_filename_read)
        
    def partial_state_prior_sample_read(self, time_index):

        if self.prior_sample_filename_read is None:
            return None

        if self.time_index != time_index:
            raise ValueError('time index mismatch')

        return readpickle(self.prior_sample_filename_read)

class DelayedSpatialComponentSolutionStorage_Files(SpaceTimeComponentSolutionStorage):
    """
    
    Class for reading and writing states/measurements from disk using dictionaries of filenames.
    
    """

    def __init__(self, 
                 statefiledictionary_read={},
                 marginal_std_filedictionary_read = {},
                 sample_filedictionary_read = {},
                 prior_sample_filedictionary_read = {},
                 measurementfiledictionary_read={},
                 statefiledictionary_write={},
                 marginal_std_filedictionary_write = {},
                 sample_filedictionary_write = {},
                 prior_sample_filedictionary_write = {},
                 measurementfiledictionary_write={}):

        super(DelayedSpatialComponentSolutionStorage_Files, self).__init__()

        self.statefiledictionary_read = self.parse_lookup(statefiledictionary_read)
        self.marginal_std_filedictionary_read = self.parse_lookup(marginal_std_filedictionary_read)
        self.sample_filedictionary_read = self.parse_lookup(sample_filedictionary_read)
        self.prior_sample_filedictionary_read = self.parse_lookup(prior_sample_filedictionary_read)
        self.measurementfiledictionary_read = self.parse_lookup(measurementfiledictionary_read)
        
        self.statefiledictionary_write = self.parse_lookup(statefiledictionary_write)
        self.marginal_std_filedictionary_write = self.parse_lookup(marginal_std_filedictionary_write)
        self.sample_filedictionary_write = self.parse_lookup(sample_filedictionary_write)
        self.prior_sample_filedictionary_write = self.parse_lookup(prior_sample_filedictionary_write)
        self.measurementfiledictionary_write = self.parse_lookup(measurementfiledictionary_write)

    def parse_lookup(self, lookup):
        """Allows input as list that is enumerated and used to create the filename lookup dictionary"""
        if lookup is None:
            return None
        elif type(lookup) is dict:
            return lookup
        elif type(lookup) is list:
            return {index: filename for index, filename in enumerate(lookup)}
        else:
            raise ValueError

    def partial_state_read(self, state_key):
        return self.state_read(state_key)
        
    def partial_state_write(self, state, state_key):
        return self.state_write(state, state_key)

    def partial_state_marginal_std_write(self, state_marginal_std, time_index):
        return self.state_marginal_std_write(state_marginal_std, time_index)

    def partial_state_marginal_std_read(self, time_index):
        return self.state_marginal_std_read(time_index)

    def state_write(self, state, state_key):

        if self.statefiledictionary_write is None:
            raise ValueError('no state write file specified')

        if state_key not in self.statefiledictionary_write:
            raise ValueError('time key mismatch')

        writepickle(state, self.statefiledictionary_write[state_key])
                
    def state_marginal_std_write(self, marginal_std, state_key):

        if self.marginal_std_filedictionary_write is None:
            raise ValueError('no marginal variances write file specified')

        if state_key not in self.marginal_std_filedictionary_write:
            raise ValueError('time key mismatch')

        writepickle(marginal_std, self.marginal_std_filedictionary_write[state_key])

    def state_sample_write(self, sample, state_key):

        if self.samplefiledictionary_write is None:
            raise ValueError('no state write file specified')

        if state_key not in self.samplefiledictionary_write:
            raise ValueError('time key mismatch')

        writepickle(state, self.samplefiledictionary_write[state_key])

    def state_read(self, state_key):
        
        if self.statefiledictionary_read is None:
            return None

        if state_key not in self.statefiledictionary_read:
            return None #raise ValueError('time key mismatch')

        return readpickle(self.statefiledictionary_read[state_key])

    def state_marginal_std_read(self, std_key):

        if self.marginal_std_filedictionary_read is None:
            return None

        if std_key not in self.marginal_std_filedictionary_read:
            return None #raise ValueError('time key mismatch')

        return readpickle(self.marginal_std_filedictionary_read[std_key])

    def state_sample_read(self, state_key):
        
        if self.samplefiledictionary_read is None:
            return None

        if state_key not in self.samplefiledictionary_read:
            return None #raise ValueError('time key mismatch')

    def measurement_write(self, measurement, measurement_key):

        if self.measurementfiledictionary_write is None:
            raise ValueError('no measurement write file specified')

        if measurement_key not in self.measurementfiledictionary_write:
            raise ValueError('time key mismatch')

        writepickle(measurement, self.measurementfiledictionary_write[measurement_key])

    def measurement_read(self, measurement_key):

        if self.measurementfiledictionary_read is None:
            raise ValueError('no measurement read file specified')

        if measurement_key not in self.measurementfiledictionary_read:
            raise ValueError('time key mismatch')

        return readpickle(self.measurementfiledictionary_read[measurement_key])

    def timeindices_read(self):

        return sorted(self.measurementfiledictionary_read.keys()) if self.measurementfiledictionary_read else None
        

from eustace.timeutils import epoch
from eumopps.timeutils import datetime_numeric

class DelayedSpatialComponentSolutionStorageFlexible_Files(SpaceTimeComponentSolutionStorage):
    """
    
    Class for reading and writing states/measurements from disk through expansion of filename patterns.
    
    """

    def __init__(self, 
                 state_read_pattern=None,
                 marginal_std_read_pattern = None,
                 sample_read_pattern=None,
                 prior_sample_read_pattern=None,
                 measurement_read_pattern=None,
                 state_write_pattern=None,
                 marginal_std_write_pattern = None,
                 sample_write_pattern=None,
                 prior_sample_write_pattern=None,
                 measurement_write_pattern=None):

        super(DelayedSpatialComponentSolutionStorageFlexible_Files, self).__init__()

        self.state_read_pattern = state_read_pattern
        self.marginal_std_read_pattern = marginal_std_read_pattern
        self.sample_read_pattern = sample_read_pattern
        self.prior_sample_read_pattern = prior_sample_read_pattern
        self.measurement_read_pattern = measurement_read_pattern

        self.state_write_pattern = state_write_pattern
        self.marginal_std_write_pattern = marginal_std_write_pattern
        self.sample_write_pattern = sample_write_pattern
        self.prior_sample_write_pattern = prior_sample_write_pattern
        self.measurement_write_pattern = measurement_write_pattern

    
    def file_write(self, values, pattern, key):
        """Write values to pickle file with filename pattern % key."""
        if pattern is None:
            raise ValueError('no filename pattern specified')
        
        filename = self.filename_from_patterns(pattern, key)
        
        writepickle(values, filename)
    
    
    def file_read(self, pattern, key):
        """Read values from pickle file with filename pattern % key.
        
        Raises a warning if file does not exist and returns None.
        
        """
        if pattern is None:
            raise ValueError('no filename pattern specified')
        
        filename = self.filename_from_patterns(pattern, key)
        
        values = None        
        if os.path.exists(filename):
            values = readpickle(filename)
        else:
            Warning('File does not exist: '+filename+' - Outputting "None"')
            
        return values

    def filename_from_patterns(self, patterns, key):

        if type(key) is int:
            t = epoch.epoch_plus_days(key)

        paths = [ datetime_numeric.build_from_pattern(pattern, t) for pattern in patterns ]
        name = os.path.join(*paths)
        return name

    def state_write(self, state, state_key):
        self.file_write(state, self.state_write_pattern, state_key)
                
    def state_marginal_std_write(self, marginal_std, std_key):
        self.file_write(marginal_std, self.marginal_std_write_pattern, std_key)

    def state_sample_write(self, sample, sample_key):
        self.file_write(sample, self.sample_write_pattern, sample_key)
        
    def prior_state_sample_write(self, sample, sample_key):
        self.file_write(sample, self.prior_sample_write_pattern, sample_key)

    def measurement_write(self, measurement, measurement_key):
        self.file_write(measurement, self.measurement_write_pattern, measurement_key)

    def state_read(self, state_key):
        return self.file_read( self.state_read_pattern, state_key )

    def state_marginal_std_read(self, std_key):
        return self.file_read( self.marginal_std_read_pattern, std_key )

    def state_sample_read(self, sample_key):
        return self.file_read( self.sample_read_pattern, sample_key )

    def prior_state_sample_read(self, sample, sample_key):
        self.file_write(sample, self.prior_sample_read_pattern, sample_key)

    def measurement_read(self, measurement_key):
        return self.file_read(self.measurement_read_pattern, measurement_key)

    def timeindices_read(self):
        raise NotImplemented
        
    def partial_state_read(self, state_key):
        return self.state_read(state_key)
        
    def partial_state_write(self, state, state_key):
        return self.state_write(state, state_key)
    
    def partial_state_marginal_std_write(self, marginal_std, std_key):
        return self.state_marginal_std_write(marginal_std, std_key)

    def partial_state_marginal_std_read(self, std_key):
        return self.state_marginal_std_read(std_key)
        
    def partial_state_sample_write(self, sample_key):
        return self.state_sample_write(sample_key)
    
    def partial_state_sample_read(self, sample_key):
        return self.state_sample_read(sample_key)
        
    def partial_state_prior_sample_write(self, sample_key):
        return self.prior_state_sample_write(sample_key)
    
    def partial_state_prior_sample_read(self, sample_key):
        return self.prior_state_sample_write(sample_key)
