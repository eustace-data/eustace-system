"""Download and parsing of DOCX format documents."""

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"

import zipfile
import io
import xml.etree.cElementTree

WORD_NAMESPACE = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
PARA = WORD_NAMESPACE + 'p'
TEXT = WORD_NAMESPACE + 't'


def retrieve_document_xml(inputfile):
    """Retrieve XML contents of MS Word document from the given input stream."""
    inputbin = io.BytesIO(inputfile.read())
    docx = zipfile.ZipFile(inputbin)
    rawxml = docx.read('word/document.xml')
    return rawxml


def parse_docx_paragraphs(rawxml, ignore_trailing_digits=False):
    """Extract paragraphs from XML of MS Word document and return as list of strings."""
    tree = xml.etree.cElementTree.XML(rawxml)
    paragraphs = []
    for paragraph in tree.getiterator(PARA):
        this_para = ''
        for node in paragraph.getiterator(TEXT):
            if node.text and (not ignore_trailing_digits or not node.text.isdigit()):
                try:
                    this_para += node.text.encode('ascii')
                except UnicodeError:
                    pass
        paragraphs.append(this_para)
    return paragraphs


def get_paragraphs_between(paragraphs, startmarker, endmarker):
    """Filter the paragraphs list and return only those between the
       first occurence of startmarker and the subsequent occurence of endmarker."""
    result = []
    mode_contents = False
    for paragraph in paragraphs:
        if mode_contents:
            if endmarker in paragraph:
                mode_contents = False
                break
            else:
                result.append(paragraph)
        elif startmarker in paragraph:
            mode_contents = True
    return result


def get_paragraphs_containing(paragraphs, datacontents, workspaceprefix):
    """Read structured data.  Assumes a set of known headings as listed in datacontents
       and a corresponding prefix identifying the paragraph following each heading
       which contains data of interest.  Returns a dictionary of heading and corresponding value."""
    result = {}
    current_data = None
    for paragraph in paragraphs:
        if paragraph in datacontents:
            current_data = paragraph
        elif workspaceprefix in paragraph:
            if current_data not in result:
                result[current_data] = paragraph
    return result
