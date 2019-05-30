"""Utilities for processing catalogues in batches, to be called from Jinja2 code in cylc suites."""

from eumopps.catalogue.jinja.batches import BatchSpecification
from eumopps.catalogue.jinja.batches import BatchSpecificationCollection

def eumopps_cylc(config, methodname, modulename=None):

    globalmethodname = 'eumopps_{methodname}'.format(methodname=methodname)
    runmethod = globals()[globalmethodname]
    return runmethod(config, modulename) if (modulename is not None) else runmethod(config)

def eumopps_cycle_points_required(specs):
    
    return BatchSpecificationCollection(specs).cyle_points_required()

def eumopps_initial_cycle_point(specs):

    return BatchSpecificationCollection(specs).initial_cycle_point()

def eumopps_final_cycle_point(specs):

    return BatchSpecificationCollection(specs).final_cycle_point()

def eumopps_parameters(specs):

    return BatchSpecificationCollection(specs).subtask_parameters()

def eumopps_subtasks(specs):

    return BatchSpecificationCollection(specs).subtasks()

def eumopps_module_scheduling(specs, modulename):

    return BatchSpecification(modulename=modulename, **(specs[modulename])).scheduling()

def eumopps_module_subtask(specs, modulename):

    return BatchSpecification(modulename=modulename, **(specs[modulename])).subtask()

def eumopps_module_subtask_name(specs, modulename):

    return BatchSpecification(modulename=modulename, **(specs[modulename])).subtask_name()
