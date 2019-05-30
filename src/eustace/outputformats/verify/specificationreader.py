"""Read specification information from documentation.  Used in tests to show that verification is done according to the latest spec."""

__version__ = "$Revision: 550 $"
__author__ = "Joel R. Mitchelson"

import numpy
import ast
from eustace.parsedoc import parsedoc
from specification import VariableSpecification

VARIABLE_TYPES = {'(float32)': numpy.float32, '(int32)': numpy.int32, '(int16)': numpy.int16, '(char)': numpy.dtype('S1'), '(uint8)': numpy.uint8}
DIMENSION_UNLIMITED = 'UNLIMITED'


class SpecificationReader(object):
    """Read specification details from document."""

    def __init__(self):
        self.paragraphs_all = []

    def load_document(self, pathname):
        """Open specified document (MS Word format) and read all paragraphs."""

        spec = open(pathname, 'r')
        spec_xml = parsedoc.retrieve_document_xml(spec)
        self.paragraphs_all = parsedoc.parse_docx_paragraphs(spec_xml)

    def read_attributes(self, startmarker, endmarker):
        """Read attributes table from document and return as dictionary."""

        # get relevant paragraphs
        paragraphs_attributes = parsedoc.get_paragraphs_between(
            self.paragraphs_all, startmarker, endmarker)

        # parse attributes table
        attributes = {}
        for paragraph in paragraphs_attributes:
            column = paragraph.split(': ')
            key = column[0]
            value = column[1]
            attributes[key] = value

        # return result
        return attributes

    def read_variables(self, startmarker, endmarker):
        """Read variables table from document and return as list of VariableSpecification objects."""

        # get relevant paragraphs
        paragraphs_variables = parsedoc.get_paragraphs_between(
            self.paragraphs_all, startmarker, endmarker)

        # build list
        variables = []
        current_variable = None
        for paragraph in paragraphs_variables:
            if paragraph.startswith('Variable: '):
                fields = paragraph.split(' ')
                name = fields[1]
                dtype = VARIABLE_TYPES[fields[2]]
                current_variable = VariableSpecification(name, dtype)
                variables.append(current_variable)
            elif paragraph:
                key = paragraph.split(': ')[0]
                value = paragraph[(len(key) + 2):]
                # allow arrays of integers as needed for flag meanings
                if value and (value[0] == '['):
                    value = ast.literal_eval(value)
                current_variable.add_metadata(key, value)

        # return result
        return variables

    def read_dimensions(self, startmarker, endmarker):
        """Read required dimensions from document and return as dictionary of required values (or None to indicate UNLIMITED)."""

        # get relevant paragraphs
        paragraphs_dimensions = parsedoc.get_paragraphs_between(
            self.paragraphs_all, startmarker, endmarker)

        # parse dimensions required
        dimensions_required = {}
        for paragraph in paragraphs_dimensions[1:]:
            column = paragraph.split(': ')
            if column[1] == DIMENSION_UNLIMITED:
                dimensions_required[column[0]] = None
            else:
                dimensions_required[column[0]] = int(column[1])

        # return result
        return dimensions_required

    def read_name_list(self, startmarker, endmarker):
        """Read a comma-separated list."""

        # Read relevant line (comma-separated list)
        namelines = parsedoc.get_paragraphs_between(
            self.paragraphs_all, startmarker, endmarker)

        # Parse
        return namelines[0].strip().split(', ')
