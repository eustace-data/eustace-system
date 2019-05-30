"""Read information from input catalogue."""

__version__ = "$Revision: 396 $"
__author__ = "Joel R. Mitchelson"


class CatalogueReader(object):
    """Abstract base class to read data set information from input catalogue."""

    def load(self, pathname):
        """Read catalogue from specified pathname."""
        raise NotImplementedError
