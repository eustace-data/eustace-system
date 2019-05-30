"""Check that station specification is up-to-date with latest documentation."""

# pylint: disable=missing-docstring,invalid-name

import unittest
import eustaceconfig
import os
from ..specificationreader import SpecificationReader
from ..specification_stations import SpecificationStationTemperature
from ..specification_stations import SpecificationStationStatus

FORMAT_SPECIFICATION_DOCUMENT_RELATIVE_PATH = 'data/docs/WP2_MS12_recommendations_data_formats.docx'


class TestSpecificationStations(unittest.TestCase):

    def test_equal_temperaturespec(self):

        # Load spec
        pathname = os.path.normpath(os.path.join(
            eustaceconfig.SYSTEM_PATH, FORMAT_SPECIFICATION_DOCUMENT_RELATIVE_PATH))
        reader = SpecificationReader()
        reader.load_document(pathname)

        # Load attributes as per global field
        attributes = reader.read_attributes(
            'The following global attributes specified by CF conventions should exist',
            'See section')

        # These are all variables as described in global fields
        globalfield_variables = reader.read_variables(
            'and so 16-bit is recommended here.',
            'NetCDF dimensions and dependent variables')

        # Filter subset
        filternames = reader.read_name_list(
            'The temperature file should contain the following variables with attributes as described',
            'The status file should contain the time variable')
        variables = [v for v in globalfield_variables if v.name in filternames]

        # Extend with station-specific ones
        variables.extend(reader.read_variables(
            'In addition to variables already mentioned, the following variables are recommended for expression of station temperatures',
            'Additional variables for station data status file'))

        # Check that they're the same as hard-coded spec used for verification
        spec = SpecificationStationTemperature()

        # Attributes
        self.assertEqual(attributes.keys(), spec.attributes.keys())
        self.assertEqual(
            [attributes[key] for key in attributes if spec.attributes[key] is not True],
            [spec.attributes[key] for key in attributes if spec.attributes[key] is not True])

        # Variables
        self.assertEqual([v.name for v in variables], [v.name for v in spec.variables])
        self.assertEqual([v.dtype for v in variables], [v.dtype for v in spec.variables])

        # Use loop just as a quick check of all metadata
        self.assertTrue(len(variables) > 0)
        for v_written in variables:
            v_spec = ([v for v in spec.variables if v.name == v_written.name])[0]
            # Note only subset of keys are used for stations
            self.assertTrue(all(key in v_written.metadata.keys() for key in v_spec.metadata.keys()))
            keys_to_check = [key for key in v_spec.metadata.keys() if v_spec.metadata[key] is not None]
            self.assertTrue(len(keys_to_check) > 0)
            meta_spec = [v_spec.metadata[key] for key in keys_to_check]
            meta_written = [v_written.metadata[key] for key in keys_to_check]
            self.assertTrue(len(meta_written) > 0)
            self.assertEqual(meta_written, meta_spec)

    def test_equal_statusspec(self):

        # Load spec
        pathname = os.path.normpath(os.path.join(eustaceconfig.SYSTEM_PATH, FORMAT_SPECIFICATION_DOCUMENT_RELATIVE_PATH))
        reader = SpecificationReader()
        reader.load_document(pathname)

        # Load attributes as per global field
        attributes = reader.read_attributes(
            'The following global attributes specified by CF conventions should exist',
            'See section')

        # Extend with station-specific ones
        variables = reader.read_variables(
            'In addition to variables already mentioned, the following variables are suggested ways to express station status',
            'Dimensions and dependencies for station temperature file')

        # Check that they're the same as hard-coded spec used for verification
        spec = SpecificationStationStatus()

        # Variables
        self.assertEqual([v.name for v in variables], [v.name for v in spec.variables])
        self.assertEqual([v.dtype for v in variables], [v.dtype for v in spec.variables])

        # Use loop just as a quick check of all metadata
        self.assertTrue(len(variables) > 0)
        for v_written in variables:
            v_spec = ([v for v in spec.variables if v.name == v_written.name])[0]
            self.assertEqual(v_written.metadata.keys(), v_spec.metadata.keys())
            keys_to_check = [key for key in v_spec.metadata.keys() if v_spec.metadata[key] is not None]
            self.assertTrue(len(keys_to_check) > 0)
            meta_spec = [v_spec.metadata[key] for key in keys_to_check]
            meta_written = [v_written.metadata[key] for key in keys_to_check]
            self.assertTrue(len(meta_written) > 0)
            self.assertEqual(meta_written, meta_spec)
