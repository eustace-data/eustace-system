"""Storage to file."""

import cPickle
from spacetime import SpaceTimeComponentSolutionStorage
from storage_files import SpaceTimeComponentSolutionStorage_Files
from eustace.outputformats.ensuredirectory import ensuredirectory

import os

def writepickle(quantity, filename):
    ensuredirectory(filename)
    with open(filename, 'wb') as fh:
        cPickle.dump(quantity, fh, protocol=2)

def readpickle(filename):
    with open(filename, 'rb') as fh:
        my_data = cPickle.load(fh)
    return my_data

class SpatialComponentSolutionStorageIndexed_Files(SpaceTimeComponentSolutionStorage):
    """
    
    Class for reading and writing states/measurements for spatial components 
    using time sorted lists of filenames.
    
    """

    def __init__(self, 
                 statefilelist_read=None,
                 marginal_std_filelist_read = None,
                 sample_filelist_read = None,
                 prior_sample_filelist_read = None,
                 measurementfilelist_read=None,
                 statefilelist_write=None,
                 marginal_std_filelist_write = None,
                 sample_filelist_write = None,
                 prior_sample_filelist_write = None,
                 measurementfilelist_write=None,
                 time_list = None):

        super(SpatialComponentSolutionStorageIndexed_Files, self).__init__()

        self.statefilelist_read = statefilelist_read
        self.marginal_std_filelist_read = marginal_std_filelist_read
        self.sample_filelist_read = sample_filelist_read
        self.prior_sample_filelist_read = prior_sample_filelist_read
        self.measurementfilelist_read = measurementfilelist_read
        
        self.statefilelist_write = statefilelist_write
        self.marginal_std_filelist_write = marginal_std_filelist_write
        self.sample_filelist_write = sample_filelist_write
        self.prior_sample_filelist_write = prior_sample_filelist_write
        self.measurementfilelist_write = measurementfilelist_write
        
        self.time_list = time_list
        
        self.file_lookup = None
        
        #print self.__dict__
        

    def number_of_substeps(self):
        """Return number of time steps represented"""
        
        if self.time_list is not None:
            return len(self.time_list)
        elif self.statefilelist_read is not None:
            return len(self.statefilelist_read)
        elif self.statefilelist_write is not None:
            return len(self.statefilelist_write)
        elif self.measurementfilelist_read is not None:
            return len(self.measurementfilelist_read)
        elif self.measurementfilelist_write is not None:
            return len(self.measurementfilelist_write)
        else:
            return 0

    def make_file_lookup(self):
        """Create lookup table of file list indices indexed by the keys provided in time_list"""

        if self.time_list is not None:
            # lookup for file index by key in time_list
            self.file_lookup = dict( zip( self.time_list, range(len(self.time_list)) ) )
        else:
            # trivial lookup where the key is the file index, i.e. {0: 0, 1: 1, ... , N: N}
            n_substeps = self.number_of_substeps()
            self.file_lookup = dict( zip(range(n_substeps), range(n_substeps) ) )
            
        #print self.__dict__
        

    def file_write(self, values, filelist, time_key):
        """Write values to pickle file with filename pattern % key."""
        
        if self.file_lookup is None:
            # Construct the file lookup table if it has not yet been produced
            self.make_file_lookup()
        
        if time_key not in self.file_lookup:
            raise ValueError('time_key not in file_lookup')
        
        if filelist is not None:
            index = self.file_lookup[time_key]
            
            if index is not None:
                filename = filelist[index]
                writepickle(values, filename) 
            else:
                raise ValueError('matching file not found')
        else:
            raise ValueError('file list not initialised')
    
    def file_read(self, filelist, time_key):
        
        if self.file_lookup is None:
            # Construct the file lookup table if it has not yet been produced
            self.make_file_lookup()
        
        if time_key not in self.file_lookup:
            raise ValueError('time_key not in file_lookup')
        
        if filelist is not None:
            index = self.file_lookup[time_key]
            
            if index is not None:
                filename = filelist[index]
                return readpickle(filename) 
                
        return None
    
    def partial_state_read(self, time_key):
        return self.state_read(time_key)
        
    def partial_state_write(self, values, time_key):
        self.state_write(values, time_key)

    def partial_state_marginal_std_read(self, time_key):
        return self.state_marginal_std_read(time_key)

    def partial_state_marginal_std_write(self, values, time_key):
        self.state_marginal_std_write(values, time_key)

    def partial_state_sample_read(self, time_key):
        return self.state_sample_read(time_key)

    def partial_state_sample_write(self, values, time_key):
        self.state_sample_write(values, time_key)
        
    def partial_state_prior_sample_read(self, time_key):
        return self.state_prior_sample_read(time_key)

    def partial_state_prior_sample_write(self, values, time_key):
        self.state_prior_sample_write(values, time_key)
        
    def state_write(self, values, time_key):
        
        if self.statefilelist_write is None:
            raise ValueError('no state write file specified')
        
        self.file_write(values, self.statefilelist_write, time_key)
                
    def state_marginal_std_write(self, values, time_key):

        if self.marginal_std_filelist_write is None:
            raise ValueError('no marginal variances write file specified')
        
        self.file_write(values, self.marginal_std_filelist_write, time_key)

    def state_sample_write(self, values, time_key):

        if self.sample_filelist_write is None:
            raise ValueError('no state write file specified')

        self.file_write(values, self.sample_filelist_write, time_key)

    def state_prior_sample_write(self, values, time_key):

        if self.prior_sample_filelist_write is None:
            raise ValueError('no state write file specified')

        self.file_write(values, self.prior_sample_filelist_write, time_key)

    def state_read(self, time_key):
        
        if self.statefilelist_read is None:
            return None

        return self.file_read(self.statefilelist_read, time_key)

    def state_marginal_std_read(self, time_key):

        if self.marginal_std_filelist_read is None:
            return None

        return self.file_read(self.marginal_std_filelist_read, time_key)

    def state_sample_read(self, time_key):
        
        if self.sample_filelist_read is None:
            return None

        return self.file_read(self.sample_filelist_read, time_key)

    def state_prior_sample_read(self, time_key):
        
        if self.prior_sample_filelist_read is None:
            return None

        return self.file_read(self.prior_sample_filelist_read, time_key)

    def measurement_write(self, values, time_key):

        if self.measurementfilelist_write is None:
            raise ValueError('no measurement write file specified')

        self.file_write(values, self.measurementfilelist_write, time_key)

    def measurement_read(self, time_key):

        if self.measurementfilelist_read is None:
            raise ValueError('no measurement read file specified')

        return self.file_read(self.measurementfilelist_read, time_key)

    def timeindices_read(self):
        return self.time_list
        
        
class SpaceTimeComponentSolutionStorageBatched_Files(SpaceTimeComponentSolutionStorage):
    """Storage of a space time component solution with measurements merged into a single measurement update"""

    def __init__(self, 
                 statefilename_read=None,
                 marginal_std_filename_read = None,
                 sample_filename_read = None,
                 prior_sample_filename_read = None,
                 statefilename_write=None,
                 marginal_std_filename_write = None,
                 sample_filename_write = None,
                 prior_sample_filename_write = None,
                 measurementfilename_write=None,
                 measurementfilelist_read=None,
                 keep_in_memory = False):

        super(SpaceTimeComponentSolutionStorageBatched_Files, self).__init__()
        
        self.statefilename_read = statefilename_read
        self.marginal_std_filename_read = marginal_std_filename_read
        self.sample_filename_read = sample_filename_read
        self.prior_sample_filename_read = prior_sample_filename_read
        
        self.statefilename_write = statefilename_write
        self.marginal_std_filename_write = marginal_std_filename_write
        self.sample_filename_write = sample_filename_write
        self.prior_sample_filename_write = prior_sample_filename_write
        
        self.measurementfilename_write = measurementfilename_write
        self.measurementfilelist_read = measurementfilelist_read

        # Storage in memory to allow reuse on read without disk access
        self.keep_in_memory = keep_in_memory
        
        self.state = None
        self.marginal_std = None
        self.sample = None
        self.prior_sample = None

    def state_write(self, state):

        if self.statefilename_write is None:
            raise ValueError('no state write file specified')

        writepickle(state, self.statefilename_write)
        
        if self.keep_in_memory:
            self.state = state

    def state_marginal_std_write(self, marginal_std):

        if self.marginal_std_filename_write is None:
            raise ValueError('no marginal variances write file specified')

        writepickle(marginal_std, self.marginal_std_filename_write)
        
        if self.keep_in_memory:
            self.marginal_std = marginal_std

    def state_sample_write(self, sample):

        if self.sample_filename_write is None:
            raise ValueError('no sample write file specified')

        writepickle(sample, self.sample_filename_write)
        
        if self.keep_in_memory:
            self.sample = sample
        
    def state_prior_sample_write(self, sample):

        if self.prior_sample_filename_write is None:
            raise ValueError('no prior sample write file specified')

        writepickle(sample, self.prior_sample_filename_write)

        if self.keep_in_memory:
            self.prior_sample = sample

    def state_read(self):

        if self.statefilename_read is None:
            return None
        #print self.statefilename_read
        if self.keep_in_memory is True:
            if self.state is None:
                self.state = readpickle(self.statefilename_read)    
            return self.state
        else:
            return readpickle(self.statefilename_read)

    def state_marginal_std_read(self):

        if self.marginal_std_filename_read is None:
            return None

        if self.keep_in_memory is True:
            if self.marginal_std is None:
                self.marginal_std = readpickle(self.marginal_std_filename_read)    
            return self.marginal_std
        else:
            return readpickle(self.marginal_std_filename_read)

    def state_sample_read(self):

        if self.sample_filename_read is None:
            return None

        if self.keep_in_memory is True:
            if self.sample is None:
                self.sample = readpickle(self.sample_filename_read)    
            return self.sample
        else:
            return readpickle(self.sample_filename_read)

    def state_prior_sample_read(self):

        if self.prior_sample_filename_read is None:
            return None

        if self.keep_in_memory is True:
            if self.prior_sample is None:
                self.prior_sample = readpickle(self.prior_sample_filename_read)    
            return self.prior_sample
        else:
            return readpickle(self.prior_sample_filename_read)

    def measurement_write(self, measurement, time_index):
        # time index is redundant - only a single measurement file is supported
        
        if self.measurementfilename_write is None:
            raise ValueError('no measurement write file specified')

        #if self.measurement_time_index_write != time_index:
            #raise ValueError('time index mismatch')

        writepickle(measurement, self.measurementfilename_write)

    def measurement_read(self, time_index):
        # time index is redundant - only a single measurement file is supported

        if self.measurementfilelist_read is None:
            raise ValueError('no measurement read list')

        #if (time_index < 0) or (time_index >= len(self.measurementfilelist_read)):
            #raise ValueError('time index of range of measurement read list')

        return readpickle(self.measurementfilename_read)

    def timeindices_read(self):
        raise NotImplementedError('Time indexing is not supported')
        #return range(len(self.measurementfilelist_read)) if self.measurementfilelist_read else None