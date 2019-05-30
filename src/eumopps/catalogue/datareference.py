"""Internal reference to datasets within a catalogue."""

class DataReference(object):
    """Trio of integer indices used to refer to an entry in a catalogue."""

    def __init__(self, dataset, subset, entry):
        self.dataset = dataset
        self.subset = subset
        self.entry = entry

    def aslist(self):
        """Represent as list of 3 values."""
        return [self.dataset, self.subset, self.entry]

def reference_list_to_integers(reflist):
    """Convert references to trios of integers."""
    result = []
    for ref in reflist:
        result.extend(ref.aslist())
    return result

def reference_integers_to_list(refnumbers):
    """Convert flat list of trios of integers to DataReference objects."""
    return [DataReference(refnumbers[index], refnumbers[index + 1], refnumbers[index + 2]) for index in range(0, len(refnumbers), 3)]

