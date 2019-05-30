"""
Surface-air ocean data
------------------------
Read surface-air model output for ocean surface as observation source.
Read WP1 NetCDF output and remap suitable for consistent output.
"""

from cubemap import ObservationSourceCubeMapGrid
from cubemap import ObservationMap
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from iris import Constraint
from iris.coords import DimCoord
from iris.coords import AuxCoord
from remap import RemapNetCDF
from remap import RemapNetCDFSpecCopy
from remap import RemapNetCDFSpecDayNumber
from remap import RemapNetCDFSpecTotalUncertainty

from eustace.outputformats.definitions import TAS, TASMIN, TASMAX, TASUNCERTAINTY, TASMINUNCERTAINTY, TASMAXUNCERTAINTY
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_RANDOM
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED2
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC2
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER0
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER1
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER2
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER3
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER4
from eustace.outputformats.definitions import OLD_OCEAN_ANCILLARY_NAMES
from eustace.surfaceairmodel.correlation_ranges import CorrelationRanges

#OCEAN details for correlated uncertainties ranges
OCEAN_CORRELATED_UNCERTAINTIES=['tas_unc_corr_sat','tas_unc_corr_mod']
OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS=[100.,1000.]
OCEAN_CORRELATED_UNCERTAINTIES_TIMES=[1.,3.]
UNITS=['km','day']

class ObservationSourceOceanMetOfficeWP1(ObservationSourceCubeMapGrid):
    """Provide ObservationSource interface to ocean surface data loaded from Iris."""

    OBSERVATIONMAPS = { 
        ObservationSource.TMEAN: ObservationMap('tas', 'unc_rand_tas', [ 'unc_corr_tas' ], [ 100000.0 ]),
    }

    def __init__(self, filename):
      
        # Default loading
        super(ObservationSourceOceanMetOfficeWP1, self).__init__(ObservationSourceOceanMetOfficeWP1.OBSERVATIONMAPS, filename)

        # And also put time coordinate onto aux coords
        time = self.extract(Constraint(cube_func=lambda cube: cube.var_name=='time'))[0]
        latitude = self.extract(Constraint(cube_func=lambda cube: cube.var_name=='latitude'))[0]
        longitude = self.extract(Constraint(cube_func=lambda cube: cube.var_name=='longitude'))[0]
        for obsmap in self.obsmaps.values():
            measurement = self.extract(Constraint(cube_func=lambda cube: cube.var_name==obsmap.measurement))[0]
            measurement.data += 273.15
            measurement.add_aux_coord(AuxCoord(points=time.data, var_name='time'))
            measurement.add_dim_coord(DimCoord(points=latitude.data, standard_name='latitude'), 0)
            measurement.add_dim_coord(DimCoord(points=longitude.data, standard_name='longitude'), 1)

class RemapNetCDFOceanWP1(RemapNetCDF):
    """Remapping of intermediate WP1-format output from ocean model suitable for consistent output."""
    
    # Adding information about correlation ranges
    ocean_correlation_ranges=CorrelationRanges(OCEAN_CORRELATED_UNCERTAINTIES,UNITS,OCEAN_CORRELATED_UNCERTAINTIES_LENGTHS,OCEAN_CORRELATED_UNCERTAINTIES_TIMES)
    ocean_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED)
    ocean_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED2)

    ANCILLARY_VARIABLES = [
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_RANDOM,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_LOCALLYCORRELATED2,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_SYSTEMATIC2,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER0,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER1,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER2,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER3,
        SURFACEAIRMODEL_ANCILLARY_OCEAN_TAS_UNCERTAINTY_PARAMETER4 
    ]

    def __init__(self):

        # Primary fields for ocean output
        mapping = [
            RemapNetCDFSpecDayNumber(),
            RemapNetCDFSpecCopy(outputname=TAS.name, inputname=TAS.name, offset=273.15),
            RemapNetCDFSpecTotalUncertainty(
                outputname=TASUNCERTAINTY.name,
                inputcomponents=OLD_OCEAN_ANCILLARY_NAMES)
        ]

        # Also map uncertainty components directly so that we can store in ancillary files
        mapping.extend([ RemapNetCDFSpecCopy(outputname=v.name, inputname=old_name) for v, old_name in zip(RemapNetCDFOceanWP1.ANCILLARY_VARIABLES, OLD_OCEAN_ANCILLARY_NAMES) ])
        # Call constructor to build class for this mapping
        super(RemapNetCDFOceanWP1, self).__init__(mapping)
