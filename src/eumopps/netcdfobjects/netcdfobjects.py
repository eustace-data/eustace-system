"""Load and save configuration objects using NetCDF"""

from netCDF4 import Dataset
from netCDF4 import stringtoarr
from netCDF4 import chartostring, stringtochar
from eumopps.timeutils.timebase import TimeBaseSeconds
from eumopps.jsonobjects import jsonobjects
from datetime import datetime
import numpy
import importlib

# Attributes used for saving
NETCDFVERSION = 'NETCDF4'
ATTRIBUTE_CONVENTIONS = 'Conventions'
CFSTANDARD = 'CF-1.6'

# Names used for storing lists of complex objects
LISTITEM_OBJECT_SUFFIX = '{0}_{1:08d}'
LIST_COUNT_NAME='list_count'
STRING_LENGTH_DIMENSION_FORMAT = 'strlen_{0}'

# Timebase used for saving
TIMEBASE = TimeBaseSeconds(datetime(1850, 1, 1))

# Recognised types for saving
STRING_TYPES = [ basestring ]
INTEGER_TYPES = [ int, long, numpy.int32, numpy.int64 ]
FLOAT_TYPES = [ float, numpy.float32, numpy.float64 ]
DATETIME_TYPES = [ datetime ]
SIMPLE_TYPES = STRING_TYPES + INTEGER_TYPES + FLOAT_TYPES + DATETIME_TYPES

def nc_open_for_load(pathname):

    nc = Dataset(pathname, 'r', format=NETCDFVERSION) # pylint: disable=invalid-name
    return nc

def nc_open_for_save(pathname):

    nc = Dataset(pathname, 'w', format=NETCDFVERSION)  # pylint: disable=invalid-name
    nc.setncattr(ATTRIBUTE_CONVENTIONS, CFSTANDARD)
    return nc

def load_object_list(parentgroup, listname):
    """Load list of objects from groups in NetCDF file."""

    numlistitems = len(parentgroup.dimensions[listname])
    listitems = []
    for listindex in range(numlistitems):
        group = parentgroup.groups[LISTITEM_OBJECT_SUFFIX.format(listname, listindex)]
        listitems.append(load_object(group))
    return listitems

def read_variable(variable):
    """Read contents of variable (will be simple type)."""

    result = None

    if variable.dtype == numpy.dtype('S1'):

        # Produces table of unicode
        result = chartostring(variable[:])

        # Convert to local string types
        if variable.ndim == 1:

            result = str(result)

        elif variable.ndim == 2:

            result = [ str(result[index][:]) for index in range(variable.shape[0]) ]

        elif variable.ndim == 3:

            result = [ [ str(result[i][j][:]) for j in range(variable.shape[1]) ] for i in range(variable.shape[0]) ]

        else:

            raise ValueError('unsupported string dimension for variable {0}'.format(str(variable)))

    elif (variable.dtype == numpy.int64) or (variable.dtype == numpy.float64):

        # Retrieve values
        result = variable[:]

        if (variable.ndim == 0) and ('units' in variable.ncattrs()) and (variable.units == TIMEBASE.units()):

            return TIMEBASE.number_to_datetime(variable[:])

        elif variable.ndim == 1:

            # Remove any masked ones
            if hasattr(result, 'mask'):
                result = result.data[~result.mask]

            # Convert to list
            result = result.tolist()

            # Convert to datetime if requested
            if ('units' in variable.ncattrs()) and (variable.units == TIMEBASE.units()):
                result = [ TIMEBASE.number_to_datetime(element) for element in result ]

        elif variable.ndim == 2:

            # Remove any masked per-row and convert to lists
            if hasattr(result, 'mask'):
                result = [ (result.data[rowindex,~result.mask[rowindex,:]]).tolist() for rowindex in range(variable.shape[0]) ]
            else:
                result = [ (result[rowindex,:]).tolist() for rowindex in range(variable.shape[0]) ]

    else:

        raise ValueError('unsupported variable {0}'.format(str(variable)))

    return result

def load_object(group):
    """Parse a group object as read from file."""

    # Empty result to fill
    result = None

    # Append all attributes
    parameters = { }
    for key in group.ncattrs():

        # Get value
        value = group.getncattr(key)

        # Use local strings (not unicode)
        if isinstance(value, unicode):
            value = str(value)

        # Set parameter
        parameters[key] = value

    # The reserved dimension name indicates it's an object list
    if LIST_COUNT_NAME in group.dimensions.keys():

        # Build list of empty dictionaries
        numrecords = group.dimensions[LIST_COUNT_NAME].size
        records = [{ k: v for (k, v) in parameters.iteritems() } for index in range(numrecords)]

        # Populate each dictionary with parameters of class to build
        for variable in group.variables.values():

            # Read all entries for this variable
            values = read_variable(variable)

            # Populate relevant class parameters
            for index, value in enumerate(values):
                records[index][variable.name] = value

        # Store as objects
        return [ makeobject(record_parameters) for record_parameters in records ]

    else:

        # Keep track of subgroups processed
        subgroups_not_processed = group.groups.keys()

        # Append variables
        for variable in group.variables.values():
            parameters[variable.name] = read_variable(variable)

        # Append subgroups corresponding to sub-lists of objects
        for dimension in group.dimensions.values():

            # dimension name is for a group list - recurse into it, and keep track of groups processed
            if LISTITEM_OBJECT_SUFFIX.format(dimension.name, 0) in group.groups:

                parameters[dimension.name] = []
                for index in range(dimension.size):

                    subgroupname = LISTITEM_OBJECT_SUFFIX.format(dimension.name, index)
                    subgroup = group.groups[subgroupname]
                    newobject = load_object(subgroup)
                    parameters[dimension.name].append(newobject)
                    subgroups_not_processed.remove(subgroupname)

        # assume any unprocessed groups are objects to be added individually
        for subgroupname in subgroups_not_processed:
            subgroup = group.groups[subgroupname]
            parameters[subgroupname] = load_object(subgroup)
	    
        # Build object
        return makeobject(parameters)

def makeobject(parameters):

    if jsonobjects.CLASSID in parameters:

        # Recognised as a class - build the class
        [ modulename, classname ] = parameters[jsonobjects.CLASSID].rsplit('.', 1)
        class_to_make = getattr(importlib.import_module(modulename), classname)
        del parameters[jsonobjects.CLASSID]
        if ATTRIBUTE_CONVENTIONS in parameters:
            del parameters[ATTRIBUTE_CONVENTIONS]
        try:
            return class_to_make(**parameters)
        except Exception as this_exception:
            print "Failed to construct: "+repr(class_to_make)
            print "with inputs: ", parameters
            raise this_exception
    
    else:

        # Not recognised as a class - assume it's just a dictionary
        return parameters

def is_type(x, typeclasses):

    return any([ isinstance(x, typeclass) for typeclass in typeclasses ])

def is_list_of(x, typeclasses):

    if isinstance(x, list):

        for typeclass in typeclasses:

            if all([ (element is None) or isinstance(element, typeclass) for element in x ]):

                return True

    return False

def is_list2d_of(x, typeclasses):

    return \
        isinstance(x, list) and \
        all([ (element is None) or is_list_of(element, typeclasses) for element in x ])

def is_simple_type_or_list_of_simple_types(x):

    return \
        is_type(x, SIMPLE_TYPES) or \
        is_list_of(x, SIMPLE_TYPES) or \
        is_list2d_of(x, STRING_TYPES)

def object_dictionary(x):

    if hasattr(x, '__dict__'):

        return x.__dict__

    elif isinstance(x, dict):

        return x

def is_simple_python_object(x):

    y = object_dictionary(x)
    return (y is not None) and all([ (value is None) or is_simple_type_or_list_of_simple_types(value) for value in y.values() ])

def is_list_of_simple_python_objects(x):

    return isinstance(x, list) and all([ (element is None) or is_simple_python_object(element) for element in x ])

def is_list_of_python_objects(x):

    return isinstance(x, list) and all([ object_dictionary(element) is not None for element in x ])

def save_simple_python_object(parentgroup, itemname, item):
    """Store a Python object by introspection.
       Supports: objects with members which are ints, floats, strings, datetimes, or lists of those.
       Definitely does not support: objects which have other objects as members.
    """

    # Make subgroup
    subgroup = parentgroup.createGroup(itemname)

    # Store object class info so that we know how to reconstruct
    if hasattr(item, '__dict__'):
        subgroup.setncattr(jsonobjects.CLASSID, jsonobjects.get_class_id(item))

    # Store class members
    for (subkey, subvalue) in object_dictionary(item).iteritems():

        if subvalue is None:

            # Nothing to write for empty object members
            pass

        elif is_type(subvalue, DATETIME_TYPES):

            # Special case for saving single date value (because must save units)
            newvariable = subgroup.createVariable(varname=subkey, datatype='i8')
            newvariable[:] = TIMEBASE.datetime_to_number(subvalue)
            newvariable.units = TIMEBASE.units()

        elif is_type(subvalue, SIMPLE_TYPES):

            # Set an individual value attribute
            subgroup.setncattr(subkey, subvalue)

        elif is_list_of(subvalue, SIMPLE_TYPES) or \
                is_list2d_of(subvalue, STRING_TYPES):

            # Store lists as member variable
            save_list_of_simple_values(subgroup, subkey, subvalue)

        else:

            raise ValueError('unsupported type for member {0}:{1}'.format(itemname, subkey))


def save_list_of_simple_python_objects(parentgroup, membername, objectlist):
    """Store a list of python objects in their own group, using introspection to extract member variables.
       All objects are expected to have the same member variables as the first one in the list.
       There will be one array variable in the group for each field of the object."""

    if (objectlist is not None) and (len(objectlist) > 0):

        # Make subgroup
        subgroup = parentgroup.createGroup(membername)

        # Dimension
        subgroup.createDimension(LIST_COUNT_NAME, len(objectlist))

        # Store object class info so that we know how to reconstruct
        if hasattr(objectlist[0], '__dict__'):
            subgroup.setncattr(jsonobjects.CLASSID, jsonobjects.get_class_id(objectlist[0]))

        # Store class members
        for subkey in object_dictionary(objectlist[0]).keys():

            # Save the list
            save_list_of_simple_values(subgroup, subkey, [ object_dictionary(item)[subkey] for item in objectlist ], dimensionname=LIST_COUNT_NAME)

def save_complex_python_object(parentgroup, membername, item):

    if membername:

        # Make subgroup
        subgroup = parentgroup.createGroup(membername)

    else:

        # Input values directly into parent group
        subgroup = parentgroup

    # Store object class info so that we know how to reconstruct
    if hasattr(item, '__dict__'):
        subgroup.setncattr(jsonobjects.CLASSID, jsonobjects.get_class_id(item))

    # Save the member variables
    for (key, value) in object_dictionary(item).iteritems():
        save_object(subgroup, key, value)

def save_list_of_complex_python_objects(parentgroup, itemname, item):

    # Need a dimensionname for the list
    parentgroup.createDimension(itemname, len(item))

    # Store each member as separate subgroup
    for index, element in enumerate(item):

        # Make element name by appending a suffix
        elementname = LISTITEM_OBJECT_SUFFIX.format(itemname, index)

        # Save this group
        save_complex_python_object(parentgroup, elementname, element)

def save_object(group, key, value):
    """Store a value or list of values."""

    if (value is None) or (isinstance(value, list) and len(value) == 0):

        # Just skip None values or empty lists
        pass

    elif is_type(value, SIMPLE_TYPES):

        # Set an individual value attribute
        group.setncattr(key, value)

    elif is_list_of(value, SIMPLE_TYPES) or \
            is_list2d_of(value, STRING_TYPES) or \
            is_list2d_of(value, INTEGER_TYPES) or \
            is_list2d_of(value, FLOAT_TYPES):

        # Store lists as member variable
        save_list_of_simple_values(group, key, value)

    elif is_list_of_simple_python_objects(value):

        # Store as group with 2D member variables (leading dimension is index in list)
        save_list_of_simple_python_objects(group, key, value)

    elif is_simple_python_object(value):

        # Store as group with member variables
        save_simple_python_object(group, key, value)

    elif is_list_of_python_objects(value):

        # Store as several subgroups (one per list element)
        save_list_of_complex_python_objects(group, key, value)

    elif hasattr(value, '__dict__') or isinstance(value, dict):

        # Store as complex python object
        save_complex_python_object(group, key, value)
    else:

        print '', key, ' = ', value

        raise ValueError('object type not supported for member {0}'.format(key))

def save_list_of_simple_values(group, membername, outputlist, dimensionname=None):
    """Store a list of values as own member variable."""

    # Leading dimension is member name if none specified
    if dimensionname is None:
        group.createDimension(membername, len(outputlist) if outputlist else 0)
        dimensionname = membername

    if (outputlist is None) or (len(outputlist) == 0) or all([ element is None for element in outputlist ]):

        pass

    elif is_list_of(outputlist, STRING_TYPES):

        save_string_list(group, membername, outputlist, dimensionname)

    elif is_list_of(outputlist, INTEGER_TYPES):

        save_numeric_list(group, membername, outputlist, dimensionname, numpy.int64)

    elif is_list_of(outputlist, FLOAT_TYPES):

        save_numeric_list(group, membername, outputlist, dimensionname, numpy.float64)

    elif is_list_of(outputlist, DATETIME_TYPES):

        save_datetime_list(group, membername, outputlist, dimensionname)

    elif is_list2d_of(outputlist, STRING_TYPES):

        save_string_list2d(group, membername, outputlist, dimensionname)

    elif is_list2d_of(outputlist, INTEGER_TYPES):

        save_numeric_list2d(group, membername, outputlist, dimensionname, numpy.int64)

    elif is_list2d_of(outputlist, FLOAT_TYPES):

        save_numeric_list2d(group, membername, outputlist, dimensionname, numpy.float64)

    else:

        print '{0}: {1}\n'.format(membername, outputlist)

        raise ValueError('type not supported for {0}'.format(membername))

def save_numeric_list(group, membername, outputlist, dimensionname, dtype):

    newvariable = group.createVariable(membername, dtype, [ dimensionname ])
    newvariable[:] = numpy.ma.masked_array(
        data=[ element if (element is not None) else 0 for element in outputlist ],
        mask=[ element is None for element in outputlist ],
        dtype=dtype)

def save_datetime_list(group, membername, outputlist, dimensionname):

    newvariable = group.createVariable(membername, numpy.int64, [ dimensionname ])
    newvariable.units = TIMEBASE.units()
    newvariable[:] = numpy.ma.masked_array(
        data=[ TIMEBASE.datetime_to_number(element) if (element is not None) else 0 for element in outputlist ],
        mask=[ element is None for element in outputlist ],
        dtype=numpy.int64)

def save_string_list(group, membername, stringarray, dimensionname):

    # Compute length of strings
    max_length = max([ len(element) if (element is not None) else 0 for element in stringarray ]) + 1

    # Name used to store length of strings
    lengthname = STRING_LENGTH_DIMENSION_FORMAT.format(membername)

    # Build new variable
    group.createDimension(lengthname, max_length)
    group.createVariable(membername, 'S1', [dimensionname, lengthname])

    # Get rid of the behaviour in which None elements map to the string 'None'
    # (we want them just to be empty strings)
    stringarray_filtered = [ s if (s is not None ) else '' for s in stringarray ]
    
    # Convert to NumPy at max_length (which includes a null terminator)
    stringrray_numpy = numpy.array(stringarray_filtered, dtype=('S' + str(max_length)))

    # Set in output
    group.variables[membername][:] = stringtochar(stringrray_numpy)

def save_string_list2d(group, membername, array2d, dimensionname):

    # lengths of all elements
    element_length = [ ]
    for element in array2d:
        if element:
            element_length.extend([ len(subelement) for subelement in element if subelement ])

    # if nothing then don't bother saving
    if len(element_length) == 0:
        return

    # Compute max length over all strings in array        
    max_length = max(element_length) + 1

    # Name used to store length of strings
    lengthname = STRING_LENGTH_DIMENSION_FORMAT.format(membername)

    # Build new variable
    group.createDimension(membername, max([ len(element) if element else 0 for element in array2d]))
    group.createDimension(lengthname, max_length)
    group.createVariable(membername, 'S1', [dimensionname, membername, lengthname])

    # populate contents
    for index, element in enumerate(array2d):
        if element is not None:
            listcontents = numpy.zeros((group.variables[membername].shape[1], max_length), 'S1')
            for subindex, subelement in enumerate(element):
                 if subelement is not None:
                     listcontents[subindex] = stringtoarr(subelement, max_length)
            group.variables[membername][index] = listcontents


def save_numeric_list2d(group, membername, array2d, dimensionname, dtype):

    numrows = len(array2d)
    maxcols = max([ len(element) if element else 0 for element in array2d])
    dimensionname2 = dimensionname + "2"
    group.createDimension(dimensionname2, maxcols)
    newvariable = group.createVariable(membername, dtype, [ dimensionname, dimensionname2 ])
    contents = numpy.ma.masked_all((numrows, maxcols), dtype=dtype)
    for index, element in enumerate(array2d):
        if element is not None:
            for subindex, subelement in enumerate(element):
                if subelement is not None:
                    contents.mask[index][subindex] = False
                    contents.data[index][subindex] = subelement
    newvariable[:] = contents
