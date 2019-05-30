"""Check that global field specification is up-to-date with latest documentation."""

# pylint: disable=missing-docstring,invalid-name

import unittest
import eustaceconfig
import os
from ..specificationreader import SpecificationReader
from ..specification_globalfield import SpecificationGlobalField

FORMAT_SPECIFICATION_DOCUMENT_RELATIVE_PATH = 'data/docs/WP2_MS12_recommendations_data_formats.docx'


class TestSpecificationGlobalField(unittest.TestCase):

    def test_equal(self):

        # Load spec
        pathname = os.path.normpath(os.path.join(
            eustaceconfig.SYSTEM_PATH, FORMAT_SPECIFICATION_DOCUMENT_RELATIVE_PATH))
        reader = SpecificationReader()
        reader.load_document(pathname)

        # Load info from written spec
        attributes = reader.read_attributes(
            'The following global attributes specified by CF conventions should exist',
            'See section')
        variables = reader.read_variables(
            'and so 16-bit is recommended here.',
            'NetCDF dimensions and dependent variables')
        dimensions = reader.read_dimensions(
            'The following layout of dimensions is recommended per daily file at 0.25 degree resolution',
            'Dependencies')

        # Load hard-coded spec
        spec = SpecificationGlobalField()

        # Check that written spec matches hard-coded spec used for verification

        # Attributes
        self.assertEqual(attributes.keys(), spec.attributes.keys())
        self.assertEqual(
            [attributes[key] for key in attributes if spec.attributes[key] is not True],
            [spec.attributes[key] for key in attributes if spec.attributes[key] is not True])

        # Dimensions
        self.assertEqual(dimensions, spec.dimensions)

        # Variables
        self.assertEqual([v.name for v in variables], [v.name for v in spec.variables])
        self.assertEqual([v.dtype for v in variables], [v.dtype for v in spec.variables])

        # Use loop just as a quick check of all metadata
        self.assertTrue(len(variables) > 0)
        for v_written in variables:
            v_spec = next(v for v in spec.variables if v.name == v_written.name)
            self.assertEqual(v_written.metadata.keys(), v_spec.metadata.keys())
            keys_to_check = [key for key in v_spec.metadata.keys() if v_spec.metadata[key] is not None]
            self.assertTrue((v_spec.name == 'timebounds') or (len(keys_to_check) > 0))
            meta_spec = [v_spec.metadata[key] for key in keys_to_check]
            meta_written = [v_written.metadata[key] for key in keys_to_check]
            self.assertTrue((v_spec.name == 'timebounds') or (len(meta_written) > 0))
            self.assertEqual(meta_written, meta_spec)
