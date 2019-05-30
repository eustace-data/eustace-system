"""Parsing of EUSTACE deliverable document 4.3."""

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

DELIVERABLE_D4_3_RELATIVE_PATH = 'data/docs/WP4_D4_3_data_management_plan.docx'
DATASET_CONTENTS_START = '8. Input Data Sets'
DATASET_CONTENTS_END = 'Output Data'

from eumopps.catalogue.dataset import CatalogueDataSet
from eustace.parsedoc.parsedoc import retrieve_document_xml
from eustace.parsedoc.parsedoc import parse_docx_paragraphs
from eustace.parsedoc.parsedoc import get_paragraphs_between
from eustace.parsedoc.parsedoc import get_paragraphs_containing
import eustaceconfig
import os


def get_dataset_list():
    """Parse and return expected array of descriptors per dataset."""
    pathname = os.path.normpath(os.path.join(
        eustaceconfig.SYSTEM_PATH, DELIVERABLE_D4_3_RELATIVE_PATH))
    inputfile = open(pathname, 'rb')
    rawxml = retrieve_document_xml(inputfile)
    paragraphs = parse_docx_paragraphs(rawxml, ignore_trailing_digits=True)
    names = get_paragraphs_between(
        paragraphs, DATASET_CONTENTS_START, DATASET_CONTENTS_END)
    paths = get_paragraphs_containing(
        paragraphs, names, eustaceconfig.WORKSPACE_PATH)
    result = [CatalogueDataSet(name, paths.get(name)) for name in names]
    return result
