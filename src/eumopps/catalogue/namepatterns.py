"""Utilities for searching and parsing of file names and paths according to specified patterns."""

# pylint: disable=anomalous-backslash-in-string,invalid-name,bad-whitespace

__version__ = "$Revision: 785 $"
__author__ = "Joel R. Mitchelson"

import re
from eumopps.timeutils import datetime_dictionary

TAG_FIELD_A = 'A'
TAG_FIELD_B = 'B'
TAG_FIELD_C = 'C'
TAG_FIELDS = [TAG_FIELD_A, TAG_FIELD_B, TAG_FIELD_C]

EXPRESSION_SUBSTITUTES = {
    'Y': (datetime_dictionary.YEAR,   '[0-9][0-9][0-9][0-9]'),
    'm': (datetime_dictionary.MONTH,  '[0-9][0-9]'),
    'd': (datetime_dictionary.DAY,    '[0-9][0-9]'),
    'H': (datetime_dictionary.HOUR,   '[0-9][0-9]'),
    'M': (datetime_dictionary.MINUTE, '[0-9][0-9]'),
    'S': (datetime_dictionary.SECOND, '[0-9][0-9]'),
    'A': (TAG_FIELD_A, '.*'),
    'B': (TAG_FIELD_B, '.*'),
    'C': (TAG_FIELD_C, '.*'),
}


def make_search_expression(pattern, knownvalues):
    """
      Build regular expression suitable for use by the re package.

      The following sequences have special meaning:

      %Y 4-digit year
      %m 2-digit month
      %d 2-digit day of month
      %H 2-digit hour
      %M 2-digit minute
      %S 2-digit seconds
      %A Dataset category (any string)
      %B Dataset sub-category (any string)
      %C Dataset sub-sub-category (any string)
      %a if %a is in the child pattern and %A has been resolved using the parent pattern, then %a is substituted by a lower-case version of the value of %A.
      %b Similar to %a
      %c Similar to %a
    """

    result = pattern
    result = result.replace('.', '\.')  # dots are taken literally
    result = result.replace('?', '.')  # question marks match any character
    for key, subs in EXPRESSION_SUBSTITUTES.iteritems():
        lookupkey = '%' + key
        result = result.replace(lookupkey, knownvalues[subs[0]]
                                if subs[0] in knownvalues else '(?P<' + subs[0] + '>' + subs[1] + ')')
    for subs in TAG_FIELDS:
        lookup_lowercase_subs = '%' + subs.lower()
        if subs in knownvalues:
            result = result.replace(lookup_lowercase_subs, knownvalues[subs].lower())
    # print '\n\npattern: ', result, '\n\n'
    return result


def hash_string_list(values):
    """Create one hashkey corresponding to an entire list of strings."""
    concat = ''
    for value in values:
        concat += value
    return concat


def hash_string_dictionary(d):
    """Create one hashkey corresponding to an entire dictionary."""
    return hash_string_list(d.values() if d else [])


class NameCandidate(object):
    """Input to pattern matching: a name and an associated flag indicating whether it has been successfully matched."""

    def __init__(self, name):
        """Initialise with name and set processing flag to unprocessed."""
        self.name = name
        self.matched = False

    def get_name(self):
        """Get associated filename."""
        return self.name

    def set_matched(self):
        """Set matched flag to True."""
        self.matched = True

    def is_matched(self):
        """Get processed flag."""
        return self.matched


class ExtractedTextFields(object):
    """A set of text fields which occur in the pathname of one file, corresponding to elements in a pathname pattern."""

    def __init__(self, name, fields):
        self.name = name
        self.fields = fields


def nested_search_patterns(patterns, candidates):
    """Extract fields corresponding to patterns in list of items,
       and return as list of ExtractedTextFields objects."""

    # support only 1 or 2 patterns
    if not patterns:
        raise ValueError('missing patterns')
    elif len(patterns) > 2:
        raise ValueError(
            'only a single pattern or two nested patterns are supported')

    # This will be populated with key-value pairs as:
    # k: hashkey of dictionary of values
    # v: a tuple giving dictionary of values, list of tuples giving item and
    #    part of item remaining
    parentvalues = {}

    # Build dictionary of known distinct values using parent pattern
    if len(patterns) > 1:
        matcher = re.compile(make_search_expression(patterns[0], {}))
        for candidate in candidates:
            outcome = matcher.match(candidate.get_name())
            if outcome:
                values = outcome.groupdict()
                key = hash_string_dictionary(values)
                remainingpart = candidate.get_name()[outcome.end():]
                entry = (candidate, remainingpart)
                if key not in parentvalues:
                    parentvalues[key] = (values, [entry])
                else:
                    parentvalues[key][1].append(entry)
    else:
        parentvalues['1'] = ({}, [(candidate, candidate.get_name()) for candidate in candidates])

    # For each distinct set of values from the parent pattern, match the child pattern
    childpattern = patterns[1] if (len(patterns) > 1) else patterns[0]
    extracts = []
    for key, keydetails in parentvalues.iteritems():
        matcher = re.compile(make_search_expression(
            childpattern, keydetails[0]))
        for entry in keydetails[1]:
            outcome = matcher.match(entry[1])
            if outcome:
                combined_outcome = dict(keydetails[0], **outcome.groupdict())
                extracts.append(ExtractedTextFields(entry[0].get_name(), combined_outcome))
                entry[0].set_matched()

    # return as a search results item
    return extracts
