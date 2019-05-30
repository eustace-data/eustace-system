"""Load and save configuration objects using JSON."""

import importlib
import json
from datetime import datetime
from eumopps.timeutils import datetime_numeric

CLASSID = 'python_class'
"""Metadata to append to JSON output to describe python class."""

DATETIMEVALUEKEY = 'timestring'
"""Key for special case of datetime objects."""

def load(inputstream):
    """Load objects from JSON inputstream."""
    return json.load(fp=inputstream, object_hook=load_object)

def save(outputstream, obj):
    """Save objects as JSON to outputstream."""
    json.dump(obj=obj, fp=outputstream, default=save_object)

def save_object(obj):
    """JSON conversion method for outputting objects.
       Appends an class type metadata so we can recreate."""

    if hasattr(obj, '__dict__'):

        # Make dictionary of all class items and recurse if needed
        result = dict((k, save_object(v)) for k, v in obj.__dict__.iteritems() if v)
    
        # Append class name and module so we can recreate
        result[CLASSID] = get_class_id(obj)

        # return this
        return result

    elif isinstance(obj, datetime):

        # Special case for datetime objects
        return { CLASSID: 'datetime.datetime', DATETIMEVALUEKEY: datetime_numeric.build(obj) }

    else:

        # Not a class type - just return
        return obj


def load_object(obj):
    """JSON conversion method for reading catalogue from file.
       If there is module and class metadata, recreate the class instance.
    """

    if isinstance(obj, dict) and (CLASSID in obj):

        # get type of new class ready to build and split into module and class names
        [ modulename, classname ] = obj[CLASSID].rsplit('.', 1)
        
        if (modulename == 'datetime') and (classname == 'datetime'):

            # special case for datetime objects
            return datetime_numeric.parse(obj[DATETIMEVALUEKEY])

        else:

            # retrieve class type to build
            class_to_make = getattr(importlib.import_module(modulename), classname)

            # build parameter dictionary and recurse if needed
            parameters = dict((k, load_object(v)) for k, v in obj.iteritems())

            # remove the unwanted metadata
            del parameters[CLASSID]

            # build it
            return class_to_make(**parameters)

    elif isinstance(obj, list):

        # recurse into list too
        return [ load_object(v) for v in obj ]

    elif isinstance(obj, unicode):

        # ensure local encoding is used
        return str(obj)

    else:

        # not dealing with a class type - just return it
        return obj

def get_class_id(obj):
    """Build class ID string for given object."""

    return obj.__class__.__module__ + '.' + obj.__class__.__name__

    
