"""Iris CubeList together with utility functions for defining global fields."""

from eustace.satgrid_iris.gridbin.gridfieldlist import GridFieldList
from eustace.satgrid_iris.gridbin.gridfieldlist import GridFieldAxis
from definitions import GLOBAL_FIELD_SHAPE
from definitions import GLOBAL_FIELD_RESOLUTION
from definitions import LATITUDE
from definitions import LONGITUDE
from definitions import TIME
import numpy

class GlobalFieldAxes2D(object):
    """Define coordinate axes for global field."""

    def __init__(self):
        """Construct with latitude and longitude properties."""

        self.latitude = GridFieldAxis(zeroth= -90.0-0.5*GLOBAL_FIELD_RESOLUTION, step=GLOBAL_FIELD_RESOLUTION, count=GLOBAL_FIELD_SHAPE[1], var_name=LATITUDE.name, circular=True, units=LATITUDE.units, standard_name=LATITUDE.standard_name, long_name=LATITUDE.long_name)

        self.longitude = GridFieldAxis(zeroth=-180.0-0.5*GLOBAL_FIELD_RESOLUTION, step=GLOBAL_FIELD_RESOLUTION, count=GLOBAL_FIELD_SHAPE[2], var_name=LONGITUDE.name, circular=True, units=LONGITUDE.units, standard_name=LONGITUDE.standard_name, long_name=LONGITUDE.long_name)

    def aslist(self):
        """Return as list of latitude and longitude objects."""

        return [ self.latitude, self.longitude ]

class GlobalFieldAxes2DWithTime(GlobalFieldAxes2D):
    """Define 2D coordinate axes for global field together with a single time point."""

    def __init__(self, fieldtime):

        # define time coord
        self.time = GridFieldAxis(zeroth=fieldtime-1.0, step=1.0, count=1, var_name=TIME.name, circular=False, units=TIME.units, standard_name=TIME.standard_name, long_name=TIME.long_name)

        # define lat/lon as in base class
        super(GlobalFieldAxes2DWithTime, self).__init__()


    def aslist(self):
        """Return as list of time, latitude, longitude objects."""

        return [ self.time, self.latitude, self.longitude ]

class GlobalFieldCubeList(GridFieldList):

    def __init__(self, fieldtime):
        """Create empty."""

        super(GlobalFieldCubeList, self).__init__(GlobalFieldAxes2DWithTime(fieldtime).aslist())
