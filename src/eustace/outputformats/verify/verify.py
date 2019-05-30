"""Verify output format against a defined specification."""

from netCDF4 import Dataset


class OutcomeDescription(object):
    """Description string and true/false outcome for a given verification check."""

    def __init__(self, description, outcome):
        self.description = description
        self.outcome = outcome


class Verifier(object):
    """Perform verification against specification."""

    # Use short name 'v' as concise way to refer to variable
    # and 'nc' as standard way to refer to NetCDF
    # pylint: disable=invalid-name

    def __init__(self):
        self.outcomes = []

    def log_outcome(self, description, outcome):
        """Add item to the outcomes table for later printing."""
        self.outcomes.append(OutcomeDescription(description, outcome))

    def get_outcomes_table(self):
        """Get the outcomes table as a string."""
        horizontal_line = '--------------------------------------------------------------------------------'
        log = ''
        log += horizontal_line + '\n'
        for result in self.outcomes:
            log += '| ' + result.description.ljust(69, ' ') + ' | ' + (
                'Pass' if result.outcome else 'Fail').ljust(4, ' ') + ' |\n'
        log += horizontal_line
        return log

    def has_passed(self):
        """Returns True if all outcomes were True, False otherwise."""
        return all([result.outcome for result in self.outcomes])

    def check_open(self, message, pathname):
        """Check file could be opened and append result to outcomes table."""

        # Attempt open
        # Catching all exceptions is sensible here as we would like
        # to append any errors to the log, no matter what they are
        # pylint: disable=broad-except
        try:
            nc = Dataset(pathname, 'r')
        except Exception:
            nc = None

        # Log any errors
        self.log_outcome(message, (nc is not None))

        # Return NetCDF instance or None if failed to open
        return nc

    def check_attributes(self, nc, expected):
        """Check attributes comply with specification and append result to outcomes table."""

        # check for attributes
        # where expected dictionary value is True then just check existence
        # where expected dictionary value is a string then check the string
        # corresponds
        for (key, value) in expected.iteritems():
            keyexists = key in nc.__dict__
            if isinstance(value, bool):
                if value is True:
                    self.log_outcome(description=(
                        'Attribute exists: ' + key), outcome=keyexists)
            elif isinstance(value, basestring):
                self.log_outcome(description=(
                    'Attribute exists: ' + key), outcome=keyexists)
                if keyexists:
                    self.log_outcome(description=(
                        'Attribute value: ' + key + '=' + value), outcome=(nc.__dict__[key] == value))

    def check_variables(self, nc, expected):
        """Check existence of at least one variable and attributes on all variables recgonised."""

        # one or more must exist
        self.log_outcome('One or more recognised variables', any([v.name in nc.variables for v in expected]))

        # check known variables
        for v in expected:
            if v.name in nc.variables:
                ncvar = nc.variables[v.name]
                if v.dtype:
                    self.log_outcome(v.name + ' ' + str(v.dtype), ncvar.dtype == v.dtype)
                for (meta_key, meta_value) in v.metadata.iteritems():
                    self.log_outcome(v.name + '.' + meta_key + ' exists', meta_key in ncvar.__dict__)
                    if meta_value is not None:
                        outcome = False
                        if meta_key in ncvar.__dict__:
                            test_value = ncvar.__getattribute__(meta_key)
                            if isinstance(meta_value, list):
                                outcome = all(meta_value == test_value)
                            else:
                                outcome = (meta_value == test_value)
                        self.log_outcome(v.name + '.' + meta_key + '=' + str(meta_value), outcome)

    def check_dimensions(self, nc, expected):
        """Check global field dimensions comply with specification."""

        # perform check
        for (key, value) in expected.iteritems():
            if value is None:
                self.log_outcome('Dimension ' + key + '=UNLIMITED',
                                 (key in nc.dimensions) and (nc.dimensions[key].isunlimited()))
            else:
                self.log_outcome('Dimension ' + key + '=' + str(value),
                                 (key in nc.dimensions) and (not nc.dimensions[key].isunlimited())
                                 and (len(nc.dimensions[key]) == value))

    def check_file(self, pathname, specification, message):
        """Check file against spec."""

        nc = self.check_open(message, pathname=pathname)
        if nc:
            self.check_attributes(nc, specification.attributes)
            self.check_variables(nc, specification.variables)
            self.check_dimensions(nc, specification.dimensions)

    def check_all(self, args):
        """Method to run all checks required - must override in derived class."""
        raise NotImplementedError

    def run(self, args, outputstream):
        """Run the check_all method using specified args and write result to output stream.  Return True if all ok, False otherwise."""
        self.check_all(args)
        outputstream.write(self.get_outcomes_table())
        return self.has_passed()
