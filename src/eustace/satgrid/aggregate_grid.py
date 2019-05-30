"""Aggregate onto grids."""

__version__ = "$Revision: 472 $"
__author__ = "Joel R. Mitchelson"

import numpy
import aggregate

class AggregateGrid(object):
    """Base class for aggregation."""

    def __init__(self):
        """Base class constructor."""
        pass

class AggregateGridCount(AggregateGrid):
    """Count number of repetitions of each grid index."""

    def __init__(self, dimensions):
        super(AggregateGridCount, self).__init__()
        self.dimensions = dimensions
        self.count = numpy.zeros(dimensions, numpy.int32)

    def aggregate_at_indices(self, indices):
        """Aggregate this set of indices."""

        aggregate.aggregate_count(self.count, indices)


class AggregateGridSum(AggregateGrid):
    """Sum observations onto grid boxes and keep count to allow calculation of mean."""

    def __init__(self, dimensions):
        super(AggregateGridSum, self).__init__()
        self.dimensions = dimensions
        self.count = numpy.zeros(dimensions, numpy.int32)
        self.sum = numpy.zeros(dimensions, numpy.float32)

    def aggregate_at_indices(self, indices, observations):
        """Sum to aggregate the given observations onto grid boxes determined by the given grid box indices."""

        count_array = self.count.flat
        sum_array = self.sum.flat
        obs_bin_iter = numpy.nditer([observations, indices])
        numberone = numpy.int32(1)

        # iterate over observations and corresponding bins
        for (obs, measurementbin) in obs_bin_iter:

            # accumulate obs
            sum_array[measurementbin] += obs

            # accumulate count at each bin
            # slightly quicker using numpy.bincount but seems to have bugs
            count_array[measurementbin] += numberone

        # bincount = numpy.bincount(bin_indices_flat)
        # count_flat[0:bincount.shape[0]] += bincount

    def get_mean(self, flag_invalid):
        """Use sum and count to compute mean,
           except where count is zero in which case result is set to invalid_flag."""

        result = numpy.empty(shape=self.dimensions, dtype=numpy.float32)
        result.fill(flag_invalid)
        valid_indices = numpy.nonzero(self.count)
        mean = self.sum[valid_indices] / self.count[valid_indices]
        result[valid_indices] = mean
        return result


class AggregateGridSumDevSq(AggregateGrid):
    """Implementation of AggregateGrid to sum onto grid boxes to facilitate calculation of mean."""

    def __init__(self, dimensions):
        super(AggregateGridSumDevSq, self).__init__()
        self.dimensions = dimensions
        self.sum_devsq = numpy.zeros(self.dimensions, numpy.float32)

    def aggregate_at_indices(self, grid_mean_observations, indices, observations):
        """Sum to aggregate the given observations onto grid boxes determined by the given grid box indices."""

        aggregate.aggregate_sum_dev_sq(self.sum_devsq, grid_mean_observations, observations, indices)


class AggregateGridArray(AggregateGrid):
    """Aggregate onto grid boxes, holding all data in-memory."""

    def __init__(self, dimensions, max_bin_obs):
        super(AggregateGridArray, self).__init__()
        self.dimensions = dimensions
        self.max_bin_obs = max_bin_obs
        arrayshape = (dimensions[0], dimensions[1], max_bin_obs)
        self.count = numpy.zeros(dimensions, numpy.int32)
        self.values = numpy.zeros(arrayshape, numpy.float32)
        self.uncertainty_random_sum_sq = numpy.zeros(dimensions, numpy.float32)
        self.uncertainty_local_atm_sum_sq = numpy.zeros(dimensions, numpy.float32)
        self.uncertainty_local_sfc_sum_sq = numpy.zeros(dimensions, numpy.float32)
        self.uncertainty_systematic_sum_sq = numpy.zeros(dimensions, numpy.float32)

    def aggregate_at_indices(self, indices, observations, uncertainty_random, uncertainty_local_atm, uncertainty_local_sfc, uncertainty_systematic):
        """Put given observations into arrays per grid-box."""

        # Allow use of 5 keyword arguments to this method
        # pylint: disable=too-many-arguments

        # And use of lots of local variables is helpful for efficiency
        # pylint: disable=too-many-locals

        values_flat = self.values.flat
        count_flat = self.count.flat
        uncertainty_random_sum_sq_flat = self.uncertainty_random_sum_sq.flat
        uncertainty_local_atm_sum_sq_flat = self.uncertainty_local_atm_sum_sq.flat
        uncertainty_local_sfc_sum_sq_flat = self.uncertainty_local_sfc_sum_sq.flat
        uncertainty_systematic_sum_sq_flat = self.uncertainty_systematic_sum_sq.flat

        bin_iter = numpy.nditer([indices, observations, uncertainty_random, uncertainty_local_atm, uncertainty_local_sfc])
        numberone = numpy.int32(1)
        uncertainty_systematic_sq = uncertainty_systematic * uncertainty_systematic

        # iterate over observations and corresponding bins
        for (measurementbin, obs, unc_ran, unc_loc_atm, unc_loc_sfc) in bin_iter:

            # get latest count
            latest_count = count_flat[measurementbin]

            # convert to index to array where all obs are kept
            arrayindex = measurementbin * self.max_bin_obs + latest_count

            # store each obs
            values_flat[arrayindex] = obs

            # accumulate count at each bin
            count_flat[measurementbin] = latest_count + numberone

            # accumulate unc_ran^2, unc_loc^2
            uncertainty_random_sum_sq_flat[measurementbin] += unc_ran * unc_ran
            uncertainty_local_atm_sum_sq_flat[measurementbin] += unc_loc_atm * unc_loc_atm
            uncertainty_local_sfc_sum_sq_flat[measurementbin] += unc_loc_sfc * unc_loc_sfc

            # sum the systematic uncertainty (effectively weighted according to
            # number of observations)
            uncertainty_systematic_sum_sq_flat[measurementbin] += uncertainty_systematic_sq

    def compute_fields(self, flag_invalid):
        """Output dictionary of computed fields."""

        # fields which will be output
        field_names = ['tsmean', 'tsmin', 'tsmax', 'tsmean_unc_ran', 'tsmean_unc_loc_atm', 'tsmean_unc_loc_sfc', 'tsmean_unc_sys']
        field_names_count = ['ts_number_of_observations', ]

        # create empty fields filled with invalid flag
        fields = {}
        for name in field_names:
            fields[name] = numpy.empty(shape=self.dimensions, dtype=numpy.float32)
            fields[name].fill(flag_invalid)

        # create zero fields for counts
        for name in field_names_count:
            fields[name] = numpy.zeros(shape=self.dimensions, dtype=numpy.int32)

        # indices with count > 0
        valid_indices = numpy.nonzero(self.count)

        # compute stats
        values_mask = (self.values == numpy.float32(0.0))
        values = numpy.ma.masked_array(self.values, mask=values_mask)

        fields['tsmean'][valid_indices] = numpy.divide(numpy.sum(self.values, axis=2)[valid_indices], self.count[valid_indices])

        # numpy.max(self.values, axis=2)[valid_indices]
        fields['tsmax'][valid_indices] = values.max(axis=2)[valid_indices]

        # numpy.min(self.values, axis=2)[valid_indices]
        fields['tsmin'][valid_indices] = values.min(axis=2)[valid_indices]

        fields['ts_number_of_observations'][
            valid_indices] = self.count[valid_indices]

        fields['tsmean_unc_ran'][valid_indices] = \
            numpy.sqrt(self.uncertainty_random_sum_sq[valid_indices] /
                       (self.count[valid_indices] * self.count[valid_indices]))

        fields['tsmean_unc_loc_atm'][valid_indices] = \
            numpy.sqrt(self.uncertainty_local_atm_sum_sq[valid_indices] /
                       self.count[valid_indices])

        fields['tsmean_unc_loc_sfc'][valid_indices] = \
            numpy.sqrt(self.uncertainty_local_sfc_sum_sq[valid_indices] /
                       self.count[valid_indices])

        fields['tsmean_unc_sys'][valid_indices] = \
            numpy.sqrt(self.uncertainty_systematic_sum_sq[valid_indices] /
                       self.count[valid_indices])

        return fields


class AggregateGridCountSumMinMax(AggregateGrid):
    """Implementation of AggregateGrid to sum onto grid boxes, holding whole array in-memory."""

    # This class has many instance attributes by nature of the computation
    # pylint: disable=too-many-instance-attributes

    def __init__(self, dimensions):
        super(AggregateGridCountSumMinMax, self).__init__()
        self.dimensions = dimensions
        self.count = numpy.zeros(dimensions, numpy.int32)
        self.values_sum = numpy.zeros(dimensions, numpy.float32)
        self.values_sum_sq = numpy.zeros(dimensions, numpy.float32)
        self.values_min = numpy.empty(dimensions, numpy.float32)
        # initialise to impossibly large value
        self.values_min.fill(numpy.finfo('f').max)
        self.values_max = numpy.empty(dimensions, numpy.float32)
        # initialise to impossibly small value
        self.values_max.fill(-numpy.finfo('f').max)
        self.uncertainty_random_sum_sq = numpy.zeros(dimensions, numpy.float32)
        self.uncertainty_local_atm_sum_sq = numpy.zeros(dimensions, numpy.float32)
        self.uncertainty_local_sfc_sum_sq = numpy.zeros(dimensions, numpy.float32)
        self.uncertainty_systematic_sum_sq = numpy.zeros(dimensions, numpy.float32)

    def aggregate_at_indices(self, indices, observations, uncertainty_random, uncertainty_local_atm, uncertainty_local_sfc, uncertainty_systematic):
        """Put given observations into arrays per grid-box."""

        # Allow use of 5 keyword arguments to this method
        # pylint: disable=too-many-arguments

        aggregate.aggregate_count_sum_min_max(self.count, self.values_sum, self.values_min, self.values_max, observations, indices)
        aggregate.aggregate_sum(self.values_sum_sq, observations * observations, indices)
        aggregate.aggregate_sum(self.uncertainty_random_sum_sq, uncertainty_random * uncertainty_random, indices)
        aggregate.aggregate_sum(self.uncertainty_local_atm_sum_sq, uncertainty_local_atm * uncertainty_local_atm, indices)
        aggregate.aggregate_sum(self.uncertainty_local_sfc_sum_sq, uncertainty_local_sfc * uncertainty_local_sfc, indices)
        aggregate.aggregate_constant(self.uncertainty_systematic_sum_sq, uncertainty_systematic * uncertainty_systematic, indices)

    def compute_fields(self, flag_invalid):
        """Output dictionary of computed fields."""

        fields = {}

        # fields which will be output
        field_names_float32 = ['tsmean', 'tsmin', 'tsmax', 'tsmean_unc_ran', 'tsmean_unc_loc_atm', 'tsmean_unc_loc_sfc', 'tsmean_unc_sys']
        field_names_count = ['ts_number_of_observations', ]

        # create empty fields filled with invalid flag
        for name in field_names_float32:
            fields[name] = numpy.empty(shape=self.dimensions, dtype=numpy.float32)
            fields[name].fill(flag_invalid)

        # create zero fields for counts
        for name in field_names_count:
            fields[name] = numpy.zeros(shape=self.dimensions, dtype=numpy.int32)

        # indices with non-zero count
        valid_indices = numpy.nonzero(self.count)

        fields['tsmean'][valid_indices] = numpy.divide(
            self.values_sum[valid_indices],
            self.count[valid_indices])

        fields['tsmax'][valid_indices] = self.values_max[valid_indices]

        fields['tsmin'][valid_indices] = self.values_min[valid_indices]

        fields['ts_number_of_observations'][valid_indices] = self.count[valid_indices]

        fields['tsmean_unc_ran'][valid_indices] = \
            numpy.sqrt(self.uncertainty_random_sum_sq[valid_indices] /
                       (self.count[valid_indices] * self.count[valid_indices]))

        fields['tsmean_unc_loc_atm'][valid_indices] = \
            numpy.sqrt(self.uncertainty_local_atm_sum_sq[valid_indices] /
                       self.count[valid_indices])

        fields['tsmean_unc_loc_sfc'][valid_indices] = \
            numpy.sqrt(self.uncertainty_local_sfc_sum_sq[valid_indices] /
                       self.count[valid_indices])

        fields['tsmean_unc_sys'][valid_indices] = \
            numpy.sqrt(self.uncertainty_systematic_sum_sq[valid_indices] /
                       self.count[valid_indices])

        return fields
