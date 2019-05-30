"""Define API for binning operations."""

import numpy

class GridBinOperation(object):
    """Abstract base class representing a binning operation."""

    def __init__(self, outputids, inputids):
        """Initialise and store identifiers for input and output ids"""
        self.inputdescriptors = { identifier: None for identifier in inputids }
        self.outputdescriptors = { identifier: None for identifier in outputids }
        self.outputinitialvalue = { identifier: 0 for identifier in outputids }

    def set_output_descriptor(self, identifier, descriptor, initialvalue=0):
        """Set details of output fields to create for the specified output variable id."""
        if not identifier in self.outputdescriptors:
            raise ValueError('identifier is not a recognised operation output')
        if not hasattr(descriptor, 'var_name'):
            raise TypeError('descriptor does not have the required attributes of class GridBinFieldDescriptor')
        self.outputdescriptors[identifier] = descriptor
        self.outputinitialvalue[identifier] = initialvalue
        return self

    def set_input_descriptor(self, identifier, descriptor):
        """Set details of input fields to use for the specified input variable id."""
        if not identifier in self.inputdescriptors:
            raise ValueError('identifier is not a recognised operation input')        
        if not hasattr(descriptor, 'var_name'):
            raise TypeError('descriptor does not have the required attributes of class GridBinFieldDescriptor')
        self.inputdescriptors[identifier] = descriptor
        return self

    def run(self, binmap):
        """Run the operate method on the given input-to-output map."""

        # check inputs and outputs all specified
        empty_identifiers = [ ]
        empty_identifiers += [ identifier for (identifier, descriptor) in self.inputdescriptors.iteritems() if descriptor is None ]
        empty_identifiers += [ identifier for (identifier, descriptor) in self.outputdescriptors.iteritems() if descriptor is None ]
        if empty_identifiers:
            raise ValueError('Missing descriptors for identifiers: ' + str(empty_identifiers))

        # Dictionary of parameters to build
        parameters = { }

        # List of output masks (will each be set to zero)
        outputmasks = [ ]

        # Build inputs from input cube
        for identifier, descriptor in self.inputdescriptors.iteritems():
            parameters[identifier] = binmap.get_input_array(descriptor).ravel()

        # Build (or make new blank) outputs
        for identifier, descriptor in self.outputdescriptors.iteritems():
            initialvalue = self.outputinitialvalue[identifier]
            field = binmap.outputgrid.get_or_create_field(descriptor, initialvalue).data
            if hasattr(field, 'mask'):
                parameters[identifier] = field.data.ravel()
                outputmasks.append(field.mask.ravel())
            else:
                parameters[identifier] = field.ravel()

        # perform the operation
        self.operate(binmap.mapfrom, binmap.mapto, **parameters)

        # flag as valid on those outputs which were produced
        for mask in outputmasks:
            mask[binmap.mapto] = False
    
    def operate(self, mapfrom, mapto, **kwargs):
        """Run the operation.  The keyword arguments should correspond to input and output ids."""
        raise NotImplementedError
