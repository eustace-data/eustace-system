"""Top-level catalogue representation."""

import uuid


class Catalogue(object):
    """
    Represent top-level catalogue of data sets and operations.

    NOTE: Although the datasets are named, they are in a list not a dictionary
    to preserve ordering, so that they can be efficiently referenced using integers.
    """

    def __init__(self, datasets=None, identifier=None, operations=None):
        """Initialise"""

        self.datasets = datasets if datasets else []
        self.operations = operations if operations else []

        if identifier:
            self.identifier = identifier
        else:
            self.newidentity()

    def newidentity(self):
        """Update catalogue identifier."""

        self.identifier = str(uuid.uuid1())

    def datasetindex(self, datasetname):
        """Index of dataset with specified name or None if no such index exists."""

        try:
            return next(index for index, dataset in enumerate(self.datasets) if dataset.name == datasetname)
        except StopIteration:
            pass

    def dataset(self, datasetname):
        """Dataset with specified name or None if no such datset exists."""

        try:
            return next(dataset for dataset in self.datasets if dataset.name == datasetname)
        except StopIteration:
            pass

    def datasubset(self, datasetname, subsetindex):
        """Specified subset of datset or None if no such index exists."""

        dataset = self.dataset(datasetname)
        if dataset:
            return dataset.subset(subsetindex)

    def number_of_subsets(self, datasetname):
        """Find number of subsets in the named data set."""

        # Search for matching name and return number of subsets
        for dataset in self.datasets:
            if dataset.name == datasetname:
                return len(dataset.subsets)

        # Not found - return zero
        return 0
