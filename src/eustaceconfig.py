"""Global configuration of EUSTACE workspace"""

__version__ = "$Revision: 1335 $"

import os

CODE_PATH = os.path.normpath(os.path.dirname(os.path.abspath(__file__)))
"""Path to the EUSTACE src folder."""

SYSTEM_PATH = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..'))
"""Path to the EUSTACE system folder."""

WORKSPACE_PATH = '/gws/nopw/j04/eustace'
"""Path to EUSTACE workspace for retrieving and storing data."""

DOCS_PATH = '/gws/nopw/j04/eustace/public/developerguide'
"""Path for output of developer guide HTML."""
