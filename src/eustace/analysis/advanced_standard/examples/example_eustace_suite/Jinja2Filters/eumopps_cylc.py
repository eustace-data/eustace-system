"""Call EUMOPPS cylc implementation."""

import eumopps.catalogue.jinja.eumopps_cylc

def eumopps_cylc(config, methodname, modulename=None):
    """Call EUMOPPS method."""

    return eumopps.catalogue.jinja.eumopps_cylc.eumopps_cylc(config, methodname, modulename)
