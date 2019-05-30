"""Input/Output of catalogue in JSON format."""

from eumopps.catalogue.dataset import CatalogueDataSet
from eumopps.catalogue.catalogue import Catalogue
from eumopps.jsonobjects import jsonobjects
from reader import CatalogueReader
from writer import CatalogueWriter

class CatalogueReaderJSON(CatalogueReader):
    """Read catalogue from JSON."""

    def load(self, pathname):
        """Load datasets where pathname corresponds to a JSON file."""

        inputstream = open(pathname, 'r')
        return self.loadstream(inputstream)

    def loadstream(self, inputstream):
        """Load datasets from the given stream."""

        return jsonobjects.load(inputstream)


class CatalogueWriterJSON(CatalogueWriter):
    """Write catalogue to JSON."""

    def __init__(self):

        super(CatalogueWriterJSON, self).__init__()

    def save(self, pathname, catalogue):
        """Open file and write pre-amble."""

        outputstream = open(pathname, 'w')
        jsonobjects.save(outputstream, catalogue)
        outputstream.flush()
        outputstream.close()
