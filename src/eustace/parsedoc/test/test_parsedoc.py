""" Unit tests for docx parsing (require test_parsedoc.docx in current working directory)."""

# pylint: disable=missing-docstring

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

import unittest
import os
from eustace.parsedoc import parsedoc
import eustaceconfig


class TestParseDoc(unittest.TestCase):

    def setUp(self):
        pathname = os.path.join(eustaceconfig.CODE_PATH, 'eustace/parsedoc/test/test_parsedoc.docx')
        inputfile = open(pathname, 'rb')
        self.rawxml = parsedoc.retrieve_document_xml(inputfile)

    def test_retrieve_document_xml(self):
        self.assertTrue(self.rawxml)

    def test_parse_docx_paragraphs(self):
        paragraphs = parsedoc.parse_docx_paragraphs(self.rawxml)
        self.assertEqual('Test Document for Parsing with Python', paragraphs[0])
        self.assertEqual('Very large tomatoes', paragraphs[5])

    def test_get_paragraphs_between(self):
        paragraphs = parsedoc.parse_docx_paragraphs(self.rawxml)
        items = parsedoc.get_paragraphs_between(
            paragraphs, 'Main items', 'Other stuff')
        self.assertEqual(['Very large tomatoes', 'Huge apricots', 'Mediocre marrows'], items)

    def test_get_paragraphs_containing(self):
        paragraphs = parsedoc.parse_docx_paragraphs(self.rawxml)
        items = parsedoc.get_paragraphs_between(paragraphs, 'Main items', 'Other stuff')
        measurements = parsedoc.get_paragraphs_containing(
            paragraphs, items, 'Measurement:')
        self.assertEqual({
            'Very large tomatoes': 'Measurement: 20cm',
            'Huge apricots': 'Measurement: 15cm',
            'Mediocre marrows': 'Measurement: 12cm'
        }, measurements)

if __name__ == '__main__':
    unittest.main()
