"""Output format for consistent air temperature fields."""

from eustace.outputformats.globalfield_filebuilder import DatasetAttributesGlobalField, FileBuilderGlobalField
from eustace.outputformats.definitions import TAS, TASMIN, TASMAX, TASUNCERTAINTY, TASMINUNCERTAINTY, TASMAXUNCERTAINTY
from eustace.analysis.observationsource import ObservationSource
import eumopps.version.svn
import numpy
import time
from post_processing import apply_land_post_processing
import eustace
from iris import load
from iris import Constraint
from iris.cube import CubeList
import numpy
from eustace.analysis.observationsource import ObservationSource
from eustace.analysis.observationsource import Observations
from cubemap import ObservationSourceCubeMapGrid
from cubemap import ObservationMap

class DatasetAttributesConsistentModelOutput(DatasetAttributesGlobalField):

    DATASET = 'global'

    def __init__(self, modulename, institution, catalogue_id):

        # Attempt to get version (may fail if SVN not in clean state)
        try:
            version = eumopps.version.svn.get_revision_id_for_module(eustace)
        except eumopps.version.svn.RepositoryStateException:
            version = eumopps.version.svn.PROVENANCE_BYPASS_WARNING

        # Construct using these parameters
        super(DatasetAttributesConsistentModelOutput, self).__init__(
            dataset=DatasetAttributesConsistentModelOutput.DATASET,
            version=version,
            mainvariable=TAS.name,
            institution=institution,
            comment='',
            history=time.strftime('%c') + ' ' + modulename,
            source='EUSTACE Catalogue ' + catalogue_id)

class ConsistentModelOutputNetCDF(object):
    """Write output from surface-air models."""
          
    PRIMARY_VARIABLES = [ TAS, TASMIN, TASMAX, TASUNCERTAINTY, TASMINUNCERTAINTY, TASMAXUNCERTAINTY ]
    FIELDNAME_DAYNUMBER = 'daynumber'

    def __init__(self, dataset_attributes, surface_dependent_variables):
        """Construct with dataset attributes ( DatasetAttributesGlobalField object ) and any surface-dependent variables to write to ancillary."""
        
        self.dataset_attributes = dataset_attributes
        self.surface_dependent_variables = surface_dependent_variables

    def write_element(self, filename, results, variables):
        """Write a single output file using given variable list."""

        # Make filebuilder
        filebuilder = FileBuilderGlobalField(filename, results[ConsistentModelOutputNetCDF.FIELDNAME_DAYNUMBER], **self.dataset_attributes.__dict__)
        
        # Add any fields defined
        for variable in variables:

            if variable.name in results:

                # Get these values
                values = results[variable.name]

                # If expressed as 2D matrix add a time dimension
                if values.ndim == 2:
                   values = numpy.ma.masked_array(data=numpy.expand_dims(values.data, axis=0), mask=numpy.expand_dims(values.mask, axis=0))

                # Add it
                filebuilder.add_global_field(variable, values)

        # Store
        filebuilder.save_and_close()

    def write_primary(self, tasfilename, results):
        """Write the results to primary consistently formatted file."""

        self.write_element(tasfilename, results, ConsistentModelOutputNetCDF.PRIMARY_VARIABLES)

    def write_ancillary(self, ancillaryfilename, results):
        """Write ancillary information according to the given surface-dependent variables."""

        self.write_element(ancillaryfilename, results, self.surface_dependent_variables)

    def post_process(self, results, surface):
	"""Apply post-processing operations depending on processed surface"""
  
	if surface == 'LAND':
	    apply_land_post_processing(results, self.surface_dependent_variables)
	else:
	    raise ValueError('Invalid surface argument. Valid values are \"ICE\", \"LAND\" and \"OCEAN\"')
        

class ObservationSourceSatstace(ObservationSourceCubeMapGrid):
    """Provide ObservationSource interface to satstace output."""

    MAPLOOKUP = {
        ObservationSource.TMEAN: ObservationMap(TAS.name, TASUNCERTAINTY.name, [ ], [ ]),
        ObservationSource.TMIN: ObservationMap(TASMIN.name, TASMINUNCERTAINTY.name, [ ], [ ]),
        ObservationSource.TMAX: ObservationMap(TASMAX.name, TASMAXUNCERTAINTY.name, [ ], [ ])
    }

    def __init__(self, observables, filename):

        # Make map of observables
        obsmap = { observable: ObservationSourceSatstace.MAPLOOKUP[observable] for observable in observables }

        # Call base class
        super(ObservationSourceSatstace, self).__init__(obsmap, filename)
