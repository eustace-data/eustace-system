"""Compare global fields between NetCDF files."""

import netCDF4
import numpy
import argparse
import os

DATA_TOLERANCE = 0.0025

class Issue(object):

    ERROR = 'ERROR'
    WARNING = 'WARNING'

    def __init__(self, status, fieldname, message):
        self.status = status
        self.fieldname = fieldname
        self.message = message

    def __repr__(self):
        return '{status}: [{fieldname}] {message}'.format(**self.__dict__)

class IssueWarning(Issue):

    def __init__(self, fieldname, message):
        
        super(IssueWarning, self).__init__(Issue.WARNING, fieldname, message)

class IssueError(Issue):

    def __init__(self, fieldname, message):
        
        super(IssueError, self).__init__(Issue.ERROR, fieldname, message)

def comparefields(baseline_filename, comparison_filename, fieldmap, rowstart=None, rowend=None, offsets=None):
    """Compare file with baseline. Dictionary fieldmap has keys which are the comparison fields to check and map to the corresponding baseline fields."""

    issues = [ ]

    if not os.path.isfile(baseline_filename):
        issues.append(IssueError(baseline_filename, 'file does not exist'))

    if not os.path.isfile(comparison_filename):
        issues.append(IssueError(comparison_filename, 'file does not exist'))

    if issues:
        return issues

    baseline_dataset = netCDF4.Dataset(baseline_filename, 'r')
    comparison_dataset = netCDF4.Dataset(comparison_filename, 'r')

    for (comparison_fieldname, baseline_fieldname) in fieldmap.iteritems():

        baseline_field = numpy.squeeze(baseline_dataset.variables[baseline_fieldname][:])
        comparison_field = numpy.squeeze(comparison_dataset.variables[comparison_fieldname][:])

        if rowend:
            comparison_field = comparison_field[0:rowend,:]

        if rowstart:
            comparison_field = comparison_field[rowstart:,:]

        if baseline_field.shape == comparison_field.shape:

            mask_diff_extra = numpy.nonzero(comparison_field.mask.ravel() & ~baseline_field.mask.ravel())[0]
            if len(mask_diff_extra) > 0:
                issues.append( IssueWarning(comparison_fieldname, 'mask mismatch {0} additional masks versus baseline'.format(len(mask_diff_extra))) )

            mask_diff_removed = numpy.nonzero(~comparison_field.mask.ravel() & baseline_field.mask.ravel())[0]
            if len(mask_diff_removed) > 0:
                issues.append( IssueError(comparison_fieldname, 'mask mismatch {0} not masked where baseline was masked'.format(len(mask_diff_removed))) )

            if offsets and comparison_fieldname in offsets:
                baseline_field += offsets[comparison_fieldname]

            valid_indices = numpy.nonzero(~comparison_field.mask.ravel())[0]
            diff_indices = numpy.nonzero(numpy.abs(baseline_field.data.ravel()[valid_indices] - comparison_field.data.ravel()[valid_indices]) > DATA_TOLERANCE)[0]
            if len(diff_indices) > 0:
                issues.append( IssueError(comparison_fieldname, 'found {0} data differences out of {1} (baseline: {2} comparison: {3})'.format(
                        len(diff_indices),
                        len(valid_indices),
                        (baseline_field.data.ravel()[valid_indices])[diff_indices],
                        (comparison_field.data.ravel()[valid_indices])[diff_indices])) )
        else:

            issues.append( IssueError(comparison_fieldname, 'shape mismatch {0} != {1}'.format(comparison_field.shape, baseline_field.shape)) )

    return issues

def main():
    """Run from command line."""

    # Parse commandline
    parser = argparse.ArgumentParser('comparefields')
    parser.add_argument('baseline_filename', help='baseline NetCDF file')
    parser.add_argument('baseline_variable', help='baseline variable name')
    parser.add_argument('comparison_filename', help='comparison NetCDF file')
    parser.add_argument('comparison_variable', help='comparison variable name')
    parser.add_argument('--rowstart', type=int, required=False, help='comparison variable name')
    parser.add_argument('--rowend', type=int, required=False, help='comparison variable name')
    args = parser.parse_args()

    # Get all arguments
    parameters = args.__dict__

    # Define fieldmap
    parameters['fieldmap'] = { args.comparison_variable: args.baseline_variable }
        
    # Remove fieldmap inputs
    del parameters['baseline_variable']
    del parameters['comparison_variable']

    # Run comparison
    issues = comparefields(**args.__dict__)

    # Results
    if issues:
        print issues
    else:
        print 'OK'

if __name__=='__main__':
    main()
