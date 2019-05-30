"""
Surface-air land source
------------------------
Read surface-air model output for land surface.
Read WP1 NetCDF output and remap suitable for consistent output.
"""

from cubemap import ObservationSourceCubeMapGrid
from cubemap import ObservationMap
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from remap import RemapNetCDF
from remap import RemapNetCDFSpecDayNumber
from remap import RemapNetCDFSpecCopy
from remap import RemapNetCDFSpecTotalUncertainty
from eustace.outputformats.definitions import TASMIN, TASMAX, TASMINUNCERTAINTY, TASMAXUNCERTAINTY
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMIN
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMIN 
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMIN 
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCSYSTASMIN
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMAX
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMAX 
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMAX 
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_UNCSYSTASMAX
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_TMINMODELNUMBER  
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_LAND_TMAXMODELNUMBER
from eustace.outputformats.definitions import OLD_LAND_ANCILLARY_NAMES
from eustace.surfaceairmodel.correlation_ranges import CorrelationRanges                     
                                                                                             
#LAND details for correlated uncertainties ranges                                            
LABELS=['tasmin','tasmax']                                                                   
SUFFIXES=['_unc_corr_sfc','_unc_corr_atm']                                                   
UNITS=['km','day']                                                                           
LAND_CORRELATED_UNCERTAINTIES=[label+suffix for label in LABELS for suffix in SUFFIXES]     
LAND_CORRELATED_UNCERTAINTIES_LENGTHS=['unknown','unknown']*len(LABELS)
LAND_CORRELATED_UNCERTAINTIES_TIMES=[30., 'unknown']*len(LABELS)

class ObservationSourceLandMetOfficeWP1(ObservationSourceCubeMapGrid):
    """Provide ObservationSource interface to land surface data loaded from Iris."""

    OBSERVATIONMAPS = { 
        ObservationSource.TMAX: ObservationMap('tasmax', 'unc_rand_tasmax', [ 'unc_corr_tasmax' ], [ 1.0 ]),
        ObservationSource.TMIN: ObservationMap('tasmin', 'unc_rand_tasmin', [ 'unc_corr_tasmin' ], [ 1.0 ])
    }

    def __init__(self, filename):
        super(ObservationSourceLandMetOfficeWP1, self).__init__(ObservationSourceLandMetOfficeWP1.OBSERVATIONMAPS, filename)

    def observables(self):
        """The names of variables estimated from this source."""

        return [ ObservationSource.TMEAN, ObservationSource.TMIN, ObservationSource.TMAX ]

    def observations(self, observable):
        """Override observations method to compute Tmean from Tmax and Tmin when requested."""

        if observable == ObservationSource.TMEAN:

            # get max and min
            tmax = self.observations(ObservationSource.TMAX)
            tmin = self.observations(ObservationSource.TMIN)

            # compute the mean
            tmean = Observations.mean(tmax, tmin)

            return tmean

        else:

            return super(ObservationSourceLandMetOfficeWP1, self).observations(observable)

    def local_correlation_length_scale(self, observable):
        """Also override to use common length scale fot Tmean."""

        # Override to take info from max
        if observable == ObservationSource.TMEAN:
            observable = ObservationSource.TMAX

        # Call base class
        return super(ObservationSourceLandMetOfficeWP1, self).local_correlation_length_scale(observable)


class RemapNetCDFLandWP1(RemapNetCDF):
    """Remapping of intermediate WP1-format output from ocean model suitable for consistent output."""

    # Adding information about correlation ranges
    land_correlation_ranges=CorrelationRanges(LAND_CORRELATED_UNCERTAINTIES,UNITS,LAND_CORRELATED_UNCERTAINTIES_LENGTHS,LAND_CORRELATED_UNCERTAINTIES_TIMES)
    land_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMIN)
    land_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMIN)
    land_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMAX)
    land_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMAX)

    VARIABLES_ANCILLARY = [
        SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMIN,
        SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMIN,
	SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMIN,
        SURFACEAIRMODEL_ANCILLARY_LAND_UNCSYSTASMIN,
        SURFACEAIRMODEL_ANCILLARY_LAND_UNCRANDTASMAX,
        SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRATMTASMAX,
	SURFACEAIRMODEL_ANCILLARY_LAND_UNCCORRSFCTASMAX,
        SURFACEAIRMODEL_ANCILLARY_LAND_UNCSYSTASMAX,
        SURFACEAIRMODEL_ANCILLARY_LAND_TMINMODELNUMBER,
        SURFACEAIRMODEL_ANCILLARY_LAND_TMAXMODELNUMBER ]

    def __init__(self):

        # Map primary variables
        mapping = [

            RemapNetCDFSpecDayNumber(),

            RemapNetCDFSpecCopy(outputname=TASMIN.name, inputname=TASMIN.name),

            RemapNetCDFSpecTotalUncertainty(
                outputname=TASMINUNCERTAINTY.name,
                inputcomponents=OLD_LAND_ANCILLARY_NAMES[:len(OLD_LAND_ANCILLARY_NAMES)/2]),

            RemapNetCDFSpecCopy(outputname=TASMAX.name, inputname=TASMAX.name),

            RemapNetCDFSpecTotalUncertainty(
                outputname=TASMAXUNCERTAINTY.name,
                inputcomponents=OLD_LAND_ANCILLARY_NAMES[len(OLD_LAND_ANCILLARY_NAMES)/2:])
        ]

        # And ancillaries too
        mapping.extend([ RemapNetCDFSpecCopy(outputname=v.name, inputname=old_name) for v, old_name in zip(RemapNetCDFLandWP1.VARIABLES_ANCILLARY, OLD_LAND_ANCILLARY_NAMES+['tasmin_model_number', 'tasmax_model_number']) ])

        # Call constructor to build class for this mapping
        super(RemapNetCDFLandWP1, self).__init__(mapping)
