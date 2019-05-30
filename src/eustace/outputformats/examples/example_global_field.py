"""Generation of example files containing global field data."""

__version__ = "$Revision: 821 $"
__author__ = "Joel R. Mitchelson"

from eustace.outputformats import definitions
from eustace.outputformats.globalfield_filebuilder import FileBuilderGlobalField
from eustace.outputformats.globalfield_filebuilder import DatasetAttributesGlobalField
from eustace.timeutils.epoch import days_since_epoch
import numpy
import time
from datetime import datetime

def write_example_global_field(source, version, outputdirectory, institution):
    """Make an example output format for a global field."""

    # set fixed seed so psuedo-random noise is actually repeatable
    numpy.random.seed(1)

    # Build attributes for example
    attributes = DatasetAttributesGlobalField(
        dataset='Example',
        version=version,
        mainvariable='tas',
        source=source,
        institution=institution,
        comment='EUSTACE project example file format for global field',
        history='Created ' + time.strftime('%c'))

    # Day number for the given date
    daynumber = int(days_since_epoch(datetime(2015, 11, 5)))

    # Make a pathname
    pathname = attributes.build_pathname(outputdirectory, daynumber)

    # object to build global field file at current time
    builder = FileBuilderGlobalField(pathname, daynumber, **attributes.__dict__)

    # get some global field data
    field_data = get_global_field_data()

    # add examples to output
    builder.add_global_field(definitions.TAS, field_data)
    builder.add_global_field(definitions.TASMIN, field_data - float(10.0))
    builder.add_global_field(definitions.TASMAX, field_data + float(10.0))

    # also get some uncertainty data
    uncertainty_data = get_uncertainty_example_data()

    # add to output
    builder.add_uncertainty_parameter('uncertainty_example', 'An example of an uncertainty variable (K)', uncertainty_data)

    # store the result
    builder.save_and_close()


def get_global_field_data():
    """Create some demonstration data."""

    # import an image of the world
    import worldmap
    global_field = (numpy.flipud(numpy.array(
        worldmap.image, numpy.float32)) * float(60.0 / 255.0)) + float(223.15)

    # superimpose logo data on map as perturbation at logo points
    global_field += numpy.multiply(numpy.random.normal(0.0, 15.0, global_field.shape), get_scaled_eustace_logo())

    # add time dimension
    return numpy.expand_dims(global_field, axis=0)


def get_uncertainty_example_data():
    """Create some demonstration data for uncertainty."""

    # retrieve logo
    logo = get_scaled_eustace_logo()

    # random noise with logo superimposed
    shape = definitions.GLOBAL_FIELD_SHAPE[1:]
    uncertainty = numpy.multiply(numpy.random.normal(0.0, 30.0, shape), logo) + \
        numpy.random.normal(0.0, 3.0, shape)

    # logo with time dimension added
    return numpy.expand_dims(uncertainty, axis=0)


def get_scaled_eustace_logo(shape=definitions.GLOBAL_FIELD_SHAPE[1:]):
    """Read EUSTACE logo file and convert to global field size."""

    # import eustace logo
    import eustacelogo
    logo = numpy.flipud(numpy.array(eustacelogo.image, numpy.float32)) / float(255.0)

    # remap to global field size
    indices_onecol = numpy.arange(0, shape[0], dtype=numpy.int32) / numpy.int32(shape[0] / logo.shape[0])
    indices_onerow = numpy.arange(0, shape[1], dtype=numpy.int32) / numpy.int32(shape[1] / logo.shape[1])
    indices_cols = numpy.kron(indices_onerow, numpy.ones((shape[0], 1), numpy.int32))
    indices_rows = numpy.kron(indices_onecol, numpy.ones((shape[1], 1), numpy.int32)).transpose()
    return logo[indices_rows, indices_cols]
