"""Utilities for processing catalogues in batches, to be called from Jinja2 code in cylc suites."""

from eumopps.catalogue.fileio.formatnetcdf import CatalogueReaderNetCDF
import os

class BatchSpecification(object):
    """Specify batch-processing of a specified module of operations."""

    def __init__(self, catalogue, modulename, batchsize, options=''):
        """Initialise."""

        self.catalogue = catalogue
        self.modulename = modulename
        self.batchsize = batchsize
        self.options = options

    def operation_count(self):
        """Query catalogue to find operations."""

        return CatalogueReaderNetCDF().quickcount(self.catalogue, self.modulename, 'operations')

    def cycle_points_required(self):
        """Get cycle points required to cover required operations, assuming number of operations per cycle point is self.batchsize."""

        return CatalogueReaderNetCDF().quickcount(self.catalogue, self.modulename, 'operations', batchsize=self.batchsize)

    def scheduling(self):
        """String for scheduling item."""

        return 'R/P1/{n}'.format(n=(self.cycle_points_required()-1))

    def subtask_index_name(self):
        """Name of index used for subtasks per cycle point. Like subtask_mymodulename"""

        return 'subtask_{modulename}'.format(modulename=self.modulename)

    def subtask_name(self):
        """Name of an indexed subtask, like mymodulename<subtask_mymodulename>."""

        return '{task}<{index}>'.format(task=self.modulename, index=self.subtask_index_name())

    def subtask_parameters(self):
        """Parameters to place in parameters section of suite definition."""

        return '{parametername} = 0..{indexmax}'.format(parametername=self.subtask_index_name(), indexmax=(self.batchsize - 1))

    def subtask(self):
        """Get cylc task section.  

        Will be like:
        [[mymodulename<subtask_mymodulename>]]
        script = python -m eumopps.catalogue.commandrun mycatalogue.nc mymodule $CYLC_TASK_PARAM_subtask_mymodulename --use_batch_size 100 --batch $CYLC_TASK_CYCLE_POINT --alow_unversioned_code
        """

        pattern = \
        '[[{subtask_name}]]\n' + \
        '        script = python -m eumopps.catalogue.commandrun {catalogue} {modulename} $CYLC_TASK_PARAM_{subtask_index_name} --use_batch_size {batchsize} --batch $CYLC_TASK_CYCLE_POINT {options}'

        return pattern.format(subtask_name=self.subtask_name(), catalogue=self.catalogue, modulename=self.modulename, subtask_index_name=self.subtask_index_name(), batchsize=self.batchsize, options=self.options)

        
class BatchSpecificationCollection(object):
    """Specify batch procesing across a collection of modules."""

    def __init__(self, specs):
        """Initialise.  Expect a dictionary with keys which are module (dataset) names and entries contain catalogue name and number of batches."""

        self.modules = { modulename: BatchSpecification(modulename=modulename, **spec) for modulename, spec in specs.iteritems() }

    def cycle_points_required(self):
        """Total batches required to cover all operations for all modules."""
        
        return max([ module.cycle_points_required() for module in self.modules.values() ])

    def initial_cycle_point(self):
        """Initial cycle point (always zero)."""

        return 0

    def final_cycle_point(self):
        """Get final cycle point (should be total - 1 )."""

        return self.cycle_points_required() - 1

    def subtask_parameters(self):
        """Parameters to use for all subtasks."""

        text = ''
        for module in self.modules.values():
            text += module.subtask_parameters() + '\n        '
        return text

    def subtasks(self):
        """Specification of all subtasks."""
        
        text = ''
        for module in self.modules.values():
            text += module.subtask() + '\n        '
        return text
