"""Compute a consistent diagnostic from all sources by aggregation onto a grid."""

import numpy
from eustace.satgrid_iris.gridbin.gridfieldlist import GridFieldList
from eustace.satgrid_iris.gridbin.gridfieldlist import GridFieldDescriptor
from diagnosticgrid_aggregate import diagnosticgrid_aggregate

class DiagnosticGridInput(object):
    """Input for diagnostic calculation: day number, grid indices, observation, (uncorrelated) uncertainty."""

    def __init__(self, observable, daynumber, indices, observation, uncertainty):
        self.observable = observable
        self.daynumber = daynumber
        self.indices = indices
        self.observation = observation
        self.uncertainty = uncertainty

class DiagnosticGridObservationConnector(object):
    """Extract SimpleSparse objects per day from observation source."""

    def __init__(self, axes, source):
        """Construct for given observation source and list of day numbers requested."""

        if (len(axes) != 2) and (len(axes) != 3):
            raise ValueError('Must have either two or three axes.')

        if (len(axes) == 3) and (axes[0].count != 1):
            raise ValueError('When using three axes the first should have one data point only')
        
        # Store params
        self.axes = axes
        self.source = source

        # Retrieve observation coordinates (common to all observables)
        location_lookup = self.source.observation_location_lookup()
        
        # Observation grid indices based on axes (common to all observables)
        indices_lat = self.axes[-2].compute_indices(location_lookup[0,:]).astype(numpy.int32)
        indices_lon = self.axes[-1].compute_indices(location_lookup[1,:]).astype(numpy.int32)
        self.grid_indices_for_location_ids = (indices_lat * numpy.int32(self.axes[-1].count)) + indices_lon


    def get_day(self, observable, daynumber):
        """Extract the specified day number from observation source."""

        # get observable
        obs = self.source.observations(observable)

        # a mask which is only zero for observations of this day
        sourcedays = obs.time.astype(numpy.int32)
        daymask = (sourcedays != daynumber).astype(numpy.bool)

        # combined mask (exclude items not on this day or unobserved)
        mask = numpy.logical_or(daymask, obs.mask)

        # array elements without mask
        valid = numpy.nonzero(numpy.logical_not(mask))

        # masked measurements and uncertainty
        observation = obs.measurement[valid].astype(numpy.float64)
        uncertainty = obs.uncorrelatederror[valid].astype(numpy.float64)

        # indices to grid (via location id)
        location_ids = obs.location[valid]
        indices = self.grid_indices_for_location_ids[location_ids]

        # create simplesparse object
        return DiagnosticGridInput(observable, daynumber, indices, observation, uncertainty)


class DiagnosticGridBins(GridFieldList):
    """Represent gridded diagnostics computed from SimpleSparse objects."""

    AGGREGATE_INITIAL_VALUE = { 
        'count': numpy.int64(0),
        'observation_max': numpy.float64(-1.0E30), 
        'observation_min': numpy.float64(1.0E30), 
        'uncertainty_max': numpy.float64(-1.0E30),
        'uncertainty_min': numpy.float64(1.0E30),
        'observation_weightedsum': numpy.float64(0.0),
        'precision_sum': numpy.float64(0.0) 
        }
    FIELDNAMEFORMAT = '{sourcename}_{observable}_{diagnostic}'

    def __init__(self, axes):
        """Create as a CubeList with specified axes to be used for all fields."""

        # The base class does what we need
        super(DiagnosticGridBins, self).__init__(axes)

    
    def create_fields_from_sparse_observations(self, sourcename, obs):
        """The obs parameter should have the SimpleSparse interface."""

        aggregate_parameters = { }
        aggregate_parameters['input_indices'] = obs.indices.astype(numpy.int32)
        aggregate_parameters['input_observation'] = obs.observation.astype(numpy.float64)
        aggregate_parameters['input_uncertainty'] = obs.uncertainty.astype(numpy.float64)
        
        # Build output fields that we will grid onto
        for diagnostic, initialvalue in DiagnosticGridBins.AGGREGATE_INITIAL_VALUE.iteritems():

            # Make a blank field (with mask)
            fieldname = DiagnosticGridBins.FIELDNAMEFORMAT.format(sourcename=sourcename, observable=obs.observable, diagnostic=diagnostic)
            descriptor = GridFieldDescriptor(var_name=fieldname, usemask=True, dtype=numpy.dtype(initialvalue))
            emptyfield = self.get_or_create_field(descriptor, initialvalue)

            # Flag as valid for those obs which are valid
            emptyfield.data.mask.ravel()[obs.indices] = False

            # Create dictionary of the data part
            aggregate_parameters[diagnostic] = emptyfield.data.data.ravel()

        # run aggregation (Cython native code)
        diagnosticgrid_aggregate(**aggregate_parameters)

    def compute_weighted_mean(self, sourcename, observable):
        """Use weighted sum of observations and sum of precision to compute weighted mean field."""

        weightedsum = self.get_field(DiagnosticGridBins.FIELDNAMEFORMAT.format(sourcename=sourcename, observable=observable, diagnostic='observation_weightedsum'))
        precision_sum = self.get_field(DiagnosticGridBins.FIELDNAMEFORMAT.format(sourcename=sourcename, observable=observable, diagnostic='precision_sum'))
        newfieldname = DiagnosticGridBins.FIELDNAMEFORMAT.format(sourcename=sourcename, observable=observable, diagnostic='observation_weightedmean')
        descriptor = GridFieldDescriptor(var_name=newfieldname, usemask=True, dtype=numpy.float64)
        contents = weightedsum.data / precision_sum.data
        self.create_field(descriptor, contents)
