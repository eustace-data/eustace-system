"""Data structures for specification of output formats
   for verification purposes."""

__version__ = "$Revision: 550 $"
__author__ = "Joel R. Mitchelson"


class VariableSpecification(object):
    """Describes a variable in the specification."""

    def __init__(self, name, dtype=None):
        self.name = name
        self.dtype = dtype
        self.metadata = {}

    def add_metadata(self, name, value):
        """Append metadata to dictionary."""
        self.metadata[name] = value
        return self

    def __str__(self):
        """Represent in a way suitable for readable log of specification."""
        result = '--' + self.name + ' ' + str(self.dtype)
        for (meta_key, meta_value) in self.metadata.iteritems():
            result += '\n    ' + meta_key + ': ' + meta_value
        return result


class FileSpecification(object):
    """Describes expected attributes, variables, and dimensions to be found in an output file."""

    def __init__(self):
        """Initialise with attributes only."""

        self.attributes = {
            'title': True,
            'institution': True,
            'source': True,
            'history': True,
            'Conventions': 'CF-1.6'
        }

        self.dimensions = {}

        self.variables = []
