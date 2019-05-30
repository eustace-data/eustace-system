"""Output variable definition."""

__version__ = "$Revision: 1326 $"
__author__ = "Joel R. Mitchelson"

import re

class OutputVariableTemplate(object):
    """An output variable template specifies any number of meta-data items for an output variable, excluding any specific short name and long name."""

    # pylint: disable=too-many-instance-attributes,too-many-arguments,too-many-locals

    # reserved items not to be included when getting all attributes
    RESERVED = ['name', 'dtype', 'fill_value', 'template_long_name' ]

    @staticmethod
    def extend(base, **kwargs):
        """ Create a new template by extending an existing one."""

        parameters =  { key: value for (key, value) in base.__dict__.iteritems() }
        for key, value in kwargs.iteritems():
            parameters[key] = value
        return OutputVariableTemplate(**parameters)

    def __init__(self, dtype, fill_value, template_long_name=None, standard_name=None, units=None, cell_methods=None, cf_role=None, flag_masks=None, flag_values=None, flag_meanings=None, valid_range=None, scale_factor=None, add_offset=None, calendar=None, bounds=None, ancillary_variables=None, axis=None, length_scale=None,length_scale_units=None,time_scale=None,time_scale_units=None):

        self.dtype = dtype
        self.fill_value = fill_value
        self.template_long_name = template_long_name
        self.standard_name = standard_name
        self.units = units
        self.cell_methods = cell_methods
        self.cf_role = cf_role
        self.flag_masks = flag_masks
        self.flag_values = flag_values
        self.flag_meanings = flag_meanings
        self.valid_range = valid_range
        self.scale_factor = scale_factor
        self.add_offset = add_offset
        self.calendar = calendar
        self.bounds = bounds
        self.ancillary_variables = ancillary_variables
        self.axis = axis
        self.length_scale=length_scale
        self.length_scale_units=length_scale_units
        self.time_scale=time_scale
        self.time_scale_units=time_scale_units


    def get_template_long_name_keywords(self):
        """If there is a template long name return any keywords in it."""

        if self.template_long_name is not None:

            # Find start and end markers for keywords
            return re.findall('\{([A-Za-z0-9_]+)\}', self.template_long_name)

        else:

            return [ ]
                
    def get_metadata(self):
        """Get meta-data members as a dictionary, excluding any reserved keywords or values which are set to None."""

        return {key: value for (key, value) in self.__dict__.iteritems() if (key not in OutputVariable.RESERVED) and (value is not None)}


class OutputVariable(OutputVariableTemplate):
    """Class to represent output variable.  This extends variable templates to include specific short name and long name."""

    @staticmethod
    def from_template(base, name, long_name=None, **kwargs):
        """Create a specifiec output variable by associating the meta-data from a template with specific short name and long name."""

        # Meta-data parameters
        parameters = base.get_metadata()

        # Might have a long name template
        keys_for_long_name = base.get_template_long_name_keywords()

        # Format long name if no override specified
        if long_name is None:
            if base.template_long_name is not None:
                long_name = base.template_long_name.format(**kwargs)
                if len(long_name) >= 2:
                    long_name = long_name[0].upper() + long_name[1:]

        # Any keyword not already used should be given to constructor
        for key, value in kwargs.iteritems():
            if key not in keys_for_long_name:
                parameters[key] = value

        return OutputVariable(name, base.dtype, base.fill_value, long_name=long_name, **parameters)
    
    def __init__(self, name, dtype, fill_value, long_name=None, **kwargs):

        super(OutputVariable, self).__init__(dtype, fill_value, **kwargs)
        self.name = name
        self.long_name = long_name
