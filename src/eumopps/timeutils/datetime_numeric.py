"""Expression of datetime as string of numerals."""

from datetime import datetime

FORMAT = "%Y%m%d%H%M%S"

def parse(inputstring):
    """Convert a string in the form YYYYmmddHHMMSS into the corresponding datetime object."""
    return datetime.strptime(inputstring, FORMAT)

def build(inputdatetime):
    """Convert a datetime object into the corresponding string in the form YYYYmmddHHMMSS."""
    iso = inputdatetime.isoformat(' ')
    return iso[0:4] + iso[5:7] + iso[8:10] + iso[11:13] + iso[14:16] + iso[17:19]

def build_from_pattern(pattern, inputdatetime):
    """Convert datetime object according to the given pattern.
       Recognises %Y, %m, %d, %H, %M, %S."""

    iso = inputdatetime.isoformat(' ')
    result = pattern
    result = result.replace('%Y', iso[0:4])
    result = result.replace('%m', iso[5:7])
    result = result.replace('%d', iso[8:10])
    result = result.replace('%H', iso[11:13])
    result = result.replace('%M', iso[14:16])
    result = result.replace('%S', iso[17:19])
    return result

