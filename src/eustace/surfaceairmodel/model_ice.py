"""
Surface-air ice model
----------------------
Ported from DMI Matlab code
"""

import numpy
import scipy.signal
import scipy.io
import sys

# Standard names
from eustace.outputformats.definitions import TAS, TASMIN, TASMAX

# Ancillary variables
from eustace.outputformats.definitions import  SURFACEAIRMODEL_ANCILLARY_ICE_TEMPLATE_UNCERTAINTY
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_EXCLUDINGCLOUD
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_UNCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_SYNOPTICALLYCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_SYSTEMATIC
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_CLOUD
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_EXCLUDINGCLOUD
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_UNCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_SYNOPTICALLYCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_SYSTEMATIC
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_CLOUD
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_EXCLUDINGCLOUD
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_UNCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_SYNOPTICALLYCORRELATED
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_SYSTEMATIC
from eustace.outputformats.definitions import SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_CLOUD
from eustace.surfaceairmodel.correlation_ranges import CorrelationRanges

#ICE details for correlated uncertainties ranges
LABELS=['tas','tasmin','tasmax']
SUFFIXES=['_unc_corr_local']
ICE_CORRELATED_UNCERTAINTIES=[label+suffix for suffix in SUFFIXES for label in LABELS ]
ICE_CORRELATED_UNCERTAINTIES_LENGTHS=[500.]*len(LABELS)
ICE_CORRELATED_UNCERTAINTIES_TIMES=[5]*len(LABELS)
UNITS=['km','day']

class ModelIce(object):
    """The ice surface-air model."""

    # Constants used by model
    LAND_ICE = 'LandIce'
    SEA_ICE = 'SeaIce'
    NORTHERN_HEMISPHERE = 'nh'
    SOUTHERN_HEMISPHERE = 'sh'

    # Names from global output definitions
    TAS = TAS.name
    TASMIN = TASMIN.name
    TASMAX = TASMAX.name

    # List of observables to calculate
    OBSERVABLES = [ TAS, TASMIN, TASMAX ]

    # Offsets about which to compute things
    # (0 Celsius, and 1/1/1900 expressed as days since epoch 1/1/1850)
    TEMPERATURE_OFFSET = 273.15
    DAYNUMBER_OFFSET = 18262.0
    DAYNUMBER_SCALE = 365.2425

    # Global large scale error excluding cloud effects
    glob_LSU = 0.2

    # Global extra cloud uncertainty on minimum temperatures
    glob_CU_offset = { TAS: 0.0, TASMIN: 2.0, TASMAX: 0.0 }

    # Sampling uncertainty determined from validation against insitu matchups
    # for each hemisphere and surface type (SH sea ice == NH sea ice)
    # Taken from original Matlab code:
    #     RU_sampling(1,1,:) = [1.6 2.2 2.1] # NH LandIce avg, min & max
    #     RU_sampling(1,2,:) = [0.0 2.9 2.0] # NH SeaIce
    #     RU_sampling(2,1,:) = [1.6 1.5 3.8] # SH LandIce
    #     RU_sampling(2,2,:) = [1.7 5.7 2.9] # SH SeaIce
    RU_sampling = {
        NORTHERN_HEMISPHERE: {
            LAND_ICE: { TAS: 1.6, TASMIN: 2.2, TASMAX: 2.1 },
            SEA_ICE : { TAS: 0.0, TASMIN: 2.9, TASMAX: 2.0 }
        },
        SOUTHERN_HEMISPHERE: {
            LAND_ICE: { TAS: 1.6, TASMIN: 1.5, TASMAX: 3.8 },
            SEA_ICE : { TAS: 1.7, TASMIN: 5.7, TASMAX: 2.9 }
        }
    }

    # Relationship uncertainty determined from regression analysis of
    # insitu-only observations, for each surface type
    # Taken from original Matlab code:
    #     SSU_relation(1,:) = [1.5 1.8 2.0]  # LandIce avg, min & max
    #     SSU_relation(2,:) = [1.7 1.8 2.1]  # SeaIce avg, min & max
    SSU_relation = {
        LAND_ICE: { TAS: 1.5, TASMIN: 1.8, TASMAX: 2.0 },
        SEA_ICE : { TAS: 1.7, TASMIN: 1.8, TASMAX: 2.1 }
    }

    COMPONENT_UNCERTAINTY_TOTAL = 'uncertainty'
    COMPONENT_UNCERTAINTY_EXCLUDING_CLOUD = '_unc_no_cloud'
    COMPONENT_UNCERTAINTY_UNCORRELATED = '_unc_rand'
    COMPONENT_UNCERTAINTY_SYNOPTICALLY_CORRELATED = '_unc_corr_local'
    COMPONENT_UNCERTAINTY_LARGE_SCALE_CORRELATED = '_unc_sys'
    COMPONENT_UNCERTAINTY_CLOUD = '_unc_cloud'

    # Adding information about correlation ranges
    ice_correlation_ranges=CorrelationRanges(ICE_CORRELATED_UNCERTAINTIES,UNITS,ICE_CORRELATED_UNCERTAINTIES_LENGTHS,ICE_CORRELATED_UNCERTAINTIES_TIMES)
    ice_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_SYNOPTICALLYCORRELATED)
    ice_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_SYNOPTICALLYCORRELATED)
    ice_correlation_ranges.update_correlated_uncertainty_ranges(SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_SYNOPTICALLYCORRELATED)

    ANCILLARY_VARIABLES = [
        SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_EXCLUDINGCLOUD,
        SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_UNCORRELATED,
        SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_SYNOPTICALLYCORRELATED,
        SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_SYSTEMATIC,
        SURFACEAIRMODEL_ANCILLARY_ICE_TAS_UNCERTAINTY_CLOUD,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_EXCLUDINGCLOUD,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_UNCORRELATED,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_SYNOPTICALLYCORRELATED,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_SYSTEMATIC,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMIN_UNCERTAINTY_CLOUD,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_EXCLUDINGCLOUD,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_UNCORRELATED,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_SYNOPTICALLYCORRELATED,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_SYSTEMATIC,
        SURFACEAIRMODEL_ANCILLARY_ICE_TASMAX_UNCERTAINTY_CLOUD
    ]

    def __init__(self, filenames=None):
        """Initialise and optionally load model parameters."""

        if filenames is not None:
            self.load_model_parameters(filenames)

    def load_model_parameters(self, filenames):
        """Filenames should be dictionary of names for each hemisphere."""

        self.components = { hemisphere: ModelIce.load_hemisphere_model_parameters(filename) for hemisphere, filename in filenames.iteritems() }

    @staticmethod
    def load_hemisphere_model_parameters(filename):
        """Load from the .mat file output by DMI model creation."""

        # Load (Matlab format)
        parameters = scipy.io.loadmat(filename)
        
        # Should be list containing 'LandIce' and 'SeaIce'
        namelist = [ str(sf_name[0]) for sf_name in parameters['sf_name'][0] ]

        # Check these are present
        if ModelIce.LAND_ICE not in namelist:
            raise ValueError('Missing parameters: ' + ModelIce.LAND_ICE)
        if ModelIce.SEA_ICE not in namelist:
            raise ValueError('Missing parameters: ' + ModelIce.SEA_ICE)

        # Dictionary of components
        components = { }

        # Should have three vectors of four numbers for every surface type
        for nameindex, name in enumerate(namelist):
            
            # retrieve components
            component_list = parameters['M'][nameindex]

            # should be three of them
            if len(component_list) != 3:
                raise ValueError('Expected 3 components for {0}, got {1}'.format(name, len(component_list)))

            # Retrieve vectors (one per output variable)
            vectors = [ component[:,0] for component in component_list ]

            # They are in order of tas, tasmin, tasmax
            components[name] = { observable: vectors[index] for (index, observable) in enumerate(ModelIce.OBSERVABLES) }

        # Return dictionary of components
        return components
            
    def process_hemisphere(self, inputs):
        """Apply model to the inputs (expect an IceSurfaceTemperatureQualityControl object for inputs)."""

        # Processor object
        processor = ModelIceProcessor(self, inputs)

        # Do the processing
        results = { }
        for observable in ModelIce.OBSERVABLES:

            # Compute for this observable and append to results dictionary
            results.update( processor.process_observable(observable) )
            
        # Return dictionary of all results
        return results

    def process_global(self, inputs_southern, inputs_northern):
        """Apply model to the two hemisphere input sets to create global output fields."""

        # Check hemisphere info
        if inputs_southern.hemisphere() != ModelIce.SOUTHERN_HEMISPHERE:
            raise ValueError('first file should be southern hemisphere')
        if inputs_northern.hemisphere() != ModelIce.NORTHERN_HEMISPHERE:
            raise ValueError('second file should be northern hemisphere')
        
        # Check day numbers compatible
        if (inputs_southern.daynumber != inputs_northern.daynumber):
            raise ValueError('mismatch in day numbers')

        # Check same number of longitudes
        if inputs_southern.longitude.shape[0] != inputs_northern.longitude.shape[0]:
            raise ValueError('mismatch in number of longitudes')

        # Retrieve latitude resolution
        latitude_resolution = inputs_southern.latitude[1] - inputs_southern.latitude[0]

        # Should be > 0 (increasing)
        if latitude_resolution <= 0.0:
            raise ValueError('latitudes expected to increase but are descreasing')

        # Check sense of sourtern hemisphere latitudes < northern hemisphere latitudes
        if inputs_northern.latitude[0] <= inputs_southern.latitude[-1]:
            raise ValueError('second file to have latitudes greater than first (but does not)')

        # Latitudes should be at same resolution for both input files (and same sign of direction of increase)
        diff_latitude0 = numpy.max(numpy.abs((inputs_southern.latitude[1:] - inputs_southern.latitude[0:-1]) - latitude_resolution))
        if diff_latitude0 > 0.000001:
            raise ValueError('latitude resolution non-constant')
        diff_latitude1 = numpy.max(numpy.abs((inputs_northern.latitude[1:] - inputs_northern.latitude[0:-1]) - latitude_resolution))
        if diff_latitude1 > 0.000001:
            raise ValueError('latitude resolution mismatch')

        # Get the first
        first_latitude = inputs_southern.latitude[0]

        # Total latitudes computed from inputs
        total_latitudes = int((inputs_northern.latitude[-1] - first_latitude + 0.000001) / latitude_resolution) + 1

        # Total longitudes from first input
        total_longitudes = inputs_southern.longitude.shape[0]

        # Results array initially just day number
        results = { 'daynumber': inputs_southern.daynumber }

        # Process each hemisphere and append to results
        for hemisphere_inputs in [ inputs_southern, inputs_northern ]:

            # Compute latitude info for this hemisphere
            num_latitudes = hemisphere_inputs.latitude.shape[0]
            latitude_index_start = int((hemisphere_inputs.latitude[0] - first_latitude + 0.000001) / latitude_resolution)
            latitude_index_end = latitude_index_start + num_latitudes

            # Compute for each hemisphere
            hemisphere_results = self.process_hemisphere(hemisphere_inputs)

            # Combine into results
            for fieldname, fielddata in hemisphere_results.iteritems():

                # Attept get from dictionary
                combined_field = results.get(fieldname, None)

                # Insert if not there
                if combined_field is None:
                    combined_field = numpy.ma.masked_all((total_latitudes, total_longitudes))
                    results[fieldname] = combined_field

                # Put hemisphere info in
                combined_field[latitude_index_start:latitude_index_end,:] = fielddata

        return results


    @staticmethod
    def compute_air_temperature_TskinSeason(coeffs, daynumber, surface_temperature):
        """Computations from TskinSeason algorithm (called by ModelIceObservableProcessor)."""

        # Offset temperature
        temperature_normalised = surface_temperature - ModelIce.TEMPERATURE_OFFSET

        # Offset daynumber
        time_normalised = 2.0 * numpy.pi * (daynumber - ModelIce.DAYNUMBER_OFFSET) / ModelIce.DAYNUMBER_SCALE

        # Taken from 'TskinSeason' algorithm
        result_normalised = coeffs[1]*temperature_normalised + (coeffs[0] + coeffs[2]*numpy.cos(time_normalised) + coeffs[3]*numpy.sin(time_normalised))

        # And put offset back to Kelvin
        result_kelvin = result_normalised + ModelIce.TEMPERATURE_OFFSET

        return result_kelvin

class ModelIceProcessor(object):
    """Processing of a given set of model inputs."""

    def __init__(self, model, inputs):

        self.model = model
        self.inputs = inputs

    def process_observable(self, observable):
        """Apply model to one input (one of ts | tsmin | tsmax)."""

        # Corresponding methods to get surface temperature input for each output observable
        corresponding_surface_temperature = {
            ModelIce.TAS   : self.inputs.quality_controlled_surface_temperature,
            ModelIce.TASMIN: self.inputs.quality_controlled_surface_temperature_min,
            ModelIce.TASMAX: self.inputs.quality_controlled_surface_temperature_max
        }

        # Masks for regions to use land ice vs sea ice model
        corresponding_part_mask = { 
            ModelIce.LAND_ICE: self.inputs.land_ice_mask, 
            ModelIce.SEA_ICE: self.inputs.sea_ice_mask
        }

        # Processor for this observable
        observable_processor = ModelIceObservableProcessor(observable)

        # The hemisphere being processed
        hemisphere = self.inputs.hemisphere()

        # Process each part
        for part in [  ModelIce.LAND_ICE, ModelIce.SEA_ICE ]:

            # Process
            observable_processor.process_observable_part(
                part,
                corresponding_part_mask[part](),
                self.model.components[hemisphere][part][observable],
                self.inputs.daynumber,
                hemisphere,
                corresponding_surface_temperature[observable](),
                self.inputs.uncertainty_uncorrelated(),
                self.inputs.uncertainty_synoptically_correlated(),
                self.inputs.uncertainty_large_scale_correlated())

        # And total
        uncertainty_total_excluding_cloud = numpy.sqrt(numpy.square(observable_processor.uncertainty_uncorrelated) + numpy.square(observable_processor.uncertainty_synoptically_correlated) + numpy.square(observable_processor.uncertainty_large_scale_correlated))

        # Total with cloud
        uncertainty_total_with_cloud = numpy.sqrt(numpy.square(observable_processor.uncertainty_uncorrelated) + numpy.square(observable_processor.uncertainty_synoptically_correlated) + numpy.square(observable_processor.uncertainty_large_scale_correlated) + numpy.square(observable_processor.uncertainty_cloud))

        # Put results into dictionary
        results = { }
        results[observable] = observable_processor.air_temperature
        results[observable+ModelIce.COMPONENT_UNCERTAINTY_TOTAL] = uncertainty_total_with_cloud
        results[observable+ModelIce.COMPONENT_UNCERTAINTY_EXCLUDING_CLOUD] = uncertainty_total_excluding_cloud
        results[observable+ModelIce.COMPONENT_UNCERTAINTY_UNCORRELATED] = observable_processor.uncertainty_uncorrelated
        results[observable+ModelIce.COMPONENT_UNCERTAINTY_SYNOPTICALLY_CORRELATED] = observable_processor.uncertainty_synoptically_correlated
        results[observable+ModelIce.COMPONENT_UNCERTAINTY_LARGE_SCALE_CORRELATED] = observable_processor.uncertainty_large_scale_correlated
        results[observable+ModelIce.COMPONENT_UNCERTAINTY_CLOUD] = observable_processor.uncertainty_cloud

        # Return results dictionary
        return results

class ModelIceObservableProcessor(object):

    def __init__(self, observable):
        """Initialise for given observable, daynumber, hemisphere, and input data."""

        # Keep parameter
        self.observable = observable

        # Arrays will be created here when first part is processed


    def process_observable_part(self, part, part_mask, coeffs, daynumber, hemisphere, surface_temperature, input_uncorrelated_uncertainty, input_synoptically_correlated_uncertainty, input_large_scale_correlated_uncertainty):
        """Apply model to one part."""

        # Allocate output arrays the first time this is called
        if not hasattr(self, 'air_temperature'):

            # Shape of output
            shape = surface_temperature.shape

            # Matrices to hold result
            self.air_temperature = numpy.ma.masked_all(shape)
            self.uncertainty_uncorrelated = numpy.ma.masked_all(shape)
            self.uncertainty_synoptically_correlated = numpy.ma.masked_all(shape)
            self.uncertainty_large_scale_correlated = numpy.ma.masked_all(shape)
            self.uncertainty_cloud = numpy.ma.masked_all(shape)

        # Valid mask is from part (ocean or sea ice) and input observations
        combined_mask = surface_temperature.mask | part_mask

        # Indices of valid temperature over this part
        indices = numpy.nonzero(combined_mask.ravel() == False)

        # Check necessary uncertainty components present
        missing = numpy.nonzero(input_uncorrelated_uncertainty.mask.ravel()[indices])[0].tolist() + \
                  numpy.nonzero(input_synoptically_correlated_uncertainty.mask.ravel()[indices])[0].tolist() + \
                  numpy.nonzero(input_large_scale_correlated_uncertainty.mask.ravel()[indices])[0].tolist()

        if len(missing):
            sys.stderr.write('WARNING: not all observations have uncertainty information (ignoring those without it).\n')
            indices = numpy.delete(indices, missing)
            # raise ValueError('ERROR: not all outputs have uncorrelated uncertainty information.')

        # Temperature over valid region
        surface_temperature_valid = surface_temperature.data.ravel()[indices]

        # Do computation
        part_air_temperature = ModelIce.compute_air_temperature_TskinSeason(coeffs, daynumber, surface_temperature_valid)

        # Uncertainty components (masked arrays)
        satellite_uncertainty_uncorrelated = input_uncorrelated_uncertainty.data.ravel()[indices]
        satellite_uncertainty_synoptically_correlated = input_synoptically_correlated_uncertainty.data.ravel()[indices]
        satellite_uncertainty_large_scale_correlated = input_large_scale_correlated_uncertainty.data.ravel()[indices]

        # Uncertainty model
        part_RU = numpy.sqrt(numpy.square(satellite_uncertainty_uncorrelated) + numpy.square(ModelIce.RU_sampling[hemisphere][part][self.observable])) * coeffs[1]
        part_SSU = numpy.sqrt(numpy.square(satellite_uncertainty_synoptically_correlated * coeffs[1]) + numpy.square(ModelIce.SSU_relation[part][self.observable]))  

        # Also constants for this part
        part_LSU = coeffs[1] * ModelIce.glob_LSU
        part_CU = coeffs[1] * (satellite_uncertainty_large_scale_correlated - ModelIce.glob_LSU + ModelIce.glob_CU_offset[self.observable])

        # Check that parts are not overlapping
        if numpy.count_nonzero(~self.air_temperature.mask.ravel()[indices]):
            raise ValueError('ERROR: overlap')

        # Put into result
        self.air_temperature.mask.ravel()[indices] = False
        self.air_temperature.data.ravel()[indices] = part_air_temperature
        self.uncertainty_uncorrelated.mask.ravel()[indices] = False
        self.uncertainty_uncorrelated.data.ravel()[indices] = part_RU
        self.uncertainty_synoptically_correlated.mask.ravel()[indices] = False
        self.uncertainty_synoptically_correlated.data.ravel()[indices] = part_SSU        
        self.uncertainty_large_scale_correlated.mask.ravel()[indices] = False
        self.uncertainty_large_scale_correlated.data.ravel()[indices] = part_LSU        
        self.uncertainty_cloud.mask.ravel()[indices] = False
        self.uncertainty_cloud.data.ravel()[indices] = part_CU
