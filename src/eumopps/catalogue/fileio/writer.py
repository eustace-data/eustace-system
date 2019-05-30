"""Write catalogue information to file."""

__version__ = "$Revision: 966 $"
__author__ = "Joel R. Mitchelson"


class CatalogueWriter(object):
    """Abstract base class to write catalogue information to file."""

    def save(self, pathname, catalogue):
        """Save all."""

        raise NotImplementedError

