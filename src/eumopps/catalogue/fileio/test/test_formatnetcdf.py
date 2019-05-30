"""Tests for Catalogue input/output in NetCDF format."""

# pylint: disable=missing-docstring,invalid-name

import unittest
from netCDF4 import Dataset
from netCDF4 import chartostring
from netCDF4 import stringtoarr
import tempfile
import shutil
from ..formatnetcdf import CatalogueReaderNetCDF
from ..formatnetcdf import CatalogueWriterNetCDF
from ...dataset import CatalogueDataSet
from ...dataset import CatalogueDataSubset
from ...dataset import CatalogueFileEntry
from ...storage import DataStorageFiles
from ...catalogue import Catalogue
from eumopps.netcdfobjects import netcdfobjects
from datetime import datetime
import numpy
import os
import subprocess

class TestCatalogueWriterNetCDF(unittest.TestCase):

    def test_write_empty(self):

        # temporary file to use for output
        testfile = tempfile.NamedTemporaryFile(suffix='.nc')

        # instantiate writer and do output
        writer = CatalogueWriterNetCDF()
        writer.save(testfile.name, Catalogue())

        # check we can load it with netCDF library
        netcdf = Dataset(testfile.name, 'r')   # pylint: disable=unused-variable

        # the objects library should recognise it as a catalogue too
        result = netcdfobjects.load_object(netcdf)
        self.assertIsInstance(result, Catalogue)

    def test_write_two_datasets(self):

        # Dictionaries representing data sets
        catalogue = Catalogue(

            identifier='bob',

            datasets=[
                
                CatalogueDataSet(
                    name='SomeData',
                    path='/where/is/it',
                    subsets = [ 
                        CatalogueDataSubset(
                            layout=DataStorageFiles(patterns=['something_fancy-with % signs and everything.nc', 'and another array member']),
                            matches=[
                                CatalogueFileEntry(name='myname.nc', time=datetime(2016, 05, 20, 13, 05), tags=['onetag', 'twotags'], size=1356, checksum='ABC129549'),
                                CatalogueFileEntry(name='othername.nc', time=datetime(2016, 05, 20, 14, 05), tags=['hello'], size=None, checksum=None),
                                CatalogueFileEntry(name='bob.nc', time=datetime(2016, 05, 20, 15, 05), tags=None, size=55552, checksum='checkitout')
                                ],
                            value_errors=['whatwentwrong?', 'we have a problem'],
                            archive_unused=['why is this even here :-/']
                            )
                        ],
                    non_matching=['/bob/is/not/here', '/neither/is/ermintrude']
                    ),
             
                CatalogueDataSet(
                    name='Other Stuff Here',
                    path='/somewhere/here',
                    subsets=[],
                    non_matching=['harrumph']),
                ]

            )

        # temporary file to use for output
        testfile = tempfile.NamedTemporaryFile(suffix='.nc')

        # instantiate writer and do output
        writer = CatalogueWriterNetCDF()
        writer.save(testfile.name, catalogue)

        # TEMP: show netcdf contents
        # print '\n' + subprocess.Popen(['ncdump', testfile.name], stdout=subprocess.PIPE).communicate()[0]

        # TEMP: show cfchecks output
        # print '\n' + subprocess.Popen(['cfchecks', testfile.name],
        # stdout=subprocess.PIPE).communicate()[0]

        # load with netCDF library
        netcdf = Dataset(testfile.name)

        # Check top level
        self.assertEqual('bob', netcdf.identifier)
        self.assertTrue('datasets' in netcdf.dimensions)
        self.assertEqual(2, len(netcdf.dimensions['datasets']))
        self.assertFalse(netcdf.dimensions['datasets'].isunlimited())
        self.assertEqual(2, len(netcdf.groups))

        # Groups for each data set
        group0 = netcdf.groups['datasets_00000000']
        group1 = netcdf.groups['datasets_00000001']

        # Check dataset group attributes
        self.assertEqual('SomeData', group0.getncattr('name'))
        self.assertEqual('/where/is/it', group0.getncattr('path'))

        # Check numbered subsets
        self.assertEqual(1, len(group0.dimensions['subsets']))
        self.assertFalse(group0.dimensions['subsets'].isunlimited())
        subset0 = group0.groups['subsets_00000000']

        # Storage spec
        # self.assertEqual('eumopps.storage.DataStorageFiles', subset0.storage)

        # Layout subgroup
        layout = subset0.groups['layout']
        self.assertEqual(2, len(layout.dimensions['patterns']))
        self.assertFalse(layout.dimensions['patterns'].isunlimited())
        self.assertEqual(2, len(layout.variables['patterns']))
        self.assertEqual('something_fancy-with % signs and everything.nc',
                         chartostring(layout.variables['patterns'][0]))
        self.assertEqual('and another array member',
                         chartostring(layout.variables['patterns'][1]))

        # Matches arrays
        matches = subset0.groups['matches']
        self.assertEqual(3, matches.dimensions['list_count'].size)
        self.assertEqual('seconds since 1850-01-01 00:00:00 UTC', matches.variables['time'].units)
        self.assertEqual(5250575100, matches.variables['time'][0])
        self.assertEqual('myname.nc', chartostring(matches.variables['name'][0]))
        self.assertEqual('onetag', chartostring(matches.variables['tags'][0][0]))
        self.assertEqual('twotags', chartostring(matches.variables['tags'][0][1]))
        self.assertIsInstance(matches.variables['size'][0], numpy.int64)
        self.assertEqual(1356, matches.variables['size'][0])
        self.assertEqual('ABC129549', chartostring(matches.variables['checksum'][0]))
        self.assertEqual('othername.nc', chartostring(matches.variables['name'][1]))
        self.assertEqual(5250578700, matches.variables['time'][1])
        self.assertEqual('hello', chartostring(matches.variables['tags'][1][0]))
        self.assertEqual('', chartostring(matches.variables['tags'][2][1]))
        self.assertTrue(matches.variables['size'][1].mask)
        self.assertEqual('', chartostring(matches.variables['checksum'][1]))
        self.assertEqual('bob.nc', chartostring(matches.variables['name'][2]))
        self.assertEqual(5250582300, matches.variables['time'][2])
        self.assertEqual('', chartostring(matches.variables['tags'][2][0]))
        self.assertEqual('', chartostring(matches.variables['tags'][2][1]))
        self.assertTrue(isinstance(matches.variables['size'][2], numpy.int64))
        self.assertEqual(55552, matches.variables['size'][2])
        self.assertEqual('checkitout', chartostring(matches.variables['checksum'][2]))

        # Value errors array of strings
        self.assertEqual(2, len(subset0.variables['value_errors']))
        self.assertEqual('whatwentwrong?', chartostring(subset0.variables['value_errors'][0]))
        self.assertEqual('we have a problem', chartostring(subset0.variables['value_errors'][1]))

        # 'archive_unused' array of strings
        self.assertEqual(1, len(subset0.variables['archive_unused']))
        self.assertEqual('why is this even here :-/', chartostring(subset0.variables['archive_unused'][0]))

        # Second data set example
        self.assertEqual(2, len(group0.dimensions['non_matching']))
        self.assertFalse(group0.dimensions['non_matching'].isunlimited())
        self.assertEqual('/bob/is/not/here', chartostring(group0.variables['non_matching'][0]))
        self.assertEqual('/neither/is/ermintrude', chartostring(group0.variables['non_matching'][1]))
        self.assertEqual('Other Stuff Here', group1.getncattr('name'))
        self.assertEqual('/somewhere/here', group1.getncattr('path'))
        self.assertEqual(1, len(group1.dimensions['non_matching']))
        self.assertFalse(group1.dimensions['non_matching'].isunlimited())
        self.assertEqual('harrumph', chartostring(group1.variables['non_matching'][0]))


class TestCatalogueReaderNetCDF(unittest.TestCase):

    def test_read_empty(self):
        """Read a file that has zero data sets."""

        # Empty file with unique temporary name
        testfile = tempfile.NamedTemporaryFile(suffix='.nc')

        # Build NetCDF manually
        nc = Dataset(testfile.name, 'w')
        nc.setncattr('Conventions', 'CF-1.6')
        nc.setncattr('identifier', 'someidentity-here')
        nc.setncattr('python_class', 'eumopps.catalogue.catalogue.Catalogue')
        nc.createDimension('datasets', 0)
        nc.createDimension('operations', 0)
        nc.close()

        # Read
        reader = CatalogueReaderNetCDF()
        result = reader.load(testfile.name)

        # Check result (should be empty dataset list and have the specified
        # identifier)
        self.assertTrue(isinstance(result.datasets, list))
        self.assertEqual(0, len(result.datasets))
        self.assertEqual('someidentity-here', result.identifier)

    def test_read_two_datasets(self):
        """Read a catalogue file for two data sets."""

        # Empty file with unique temporary name
        testfile = tempfile.NamedTemporaryFile(suffix='.nc')

        # Build NetCDF manually
        nc = Dataset(testfile.name, 'w')
        nc.setncattr('Conventions', 'CF-1.6')
        nc.setncattr('identifier', 'citizen-of-the-world')
        nc.setncattr('python_class', 'eumopps.catalogue.catalogue.Catalogue')
        nc.createDimension('datasets', 2)
        nc.createDimension('default_strlen', 23)
        group0 = nc.createGroup('datasets_00000000')
        group0.setncattr('name', 'Ffflip')
        group0.setncattr('path', '/find/it/here')
        group0.setncattr('python_class', 'eumopps.catalogue.dataset.CatalogueDataSet')
        group1 = nc.createGroup('datasets_00000001')
        group1.setncattr('name', 'Ffflop')
        group1.setncattr('path', '/or/here')
        group1.setncattr('python_class', 'eumopps.catalogue.dataset.CatalogueDataSet')
        group1.createDimension('subsets', 1)
        subset = group1.createGroup('subsets_00000000')
        subset.setncattr('python_class', 'eumopps.catalogue.dataset.CatalogueDataSubset')
        layout = subset.createGroup('layout')
        layout.setncattr('python_class', 'eumopps.catalogue.storage.DataStorageFiles')
        layout.createDimension('patterns', 2)
        layout.createVariable('patterns', 'S1', ['patterns', 'default_strlen'])
        layout.variables['patterns'][0] = stringtoarr('splendid', 23)
        layout.variables['patterns'][1] = stringtoarr('pretty', 23)
        matches = subset.createGroup('matches')
        matches.setncattr('python_class', 'eumopps.catalogue.dataset.CatalogueFileEntry')
        matches.createDimension('list_count', 0)
        matches.createDimension('tags', 2)
        matches_name = matches.createVariable('name', 'S1', ['list_count', 'default_strlen'])
        matches_time = matches.createVariable('time', 'i8', ['list_count'])
        matches_time.units = 'seconds since 1850-01-01 00:00:00 UTC'
        matches_size = matches.createVariable('size', 'i8', ['list_count'])
        matches_tags = matches.createVariable('tags', 'S1', ['list_count', 'tags', 'default_strlen'])
        matches_name[0] = stringtoarr('bob', 23)
        # 2016-05-30 17:29:33 = 60780 days * 86400  + 17 hours * 3600 + 29 minutes * 60 + 33
        matches_time[0] = 5251454973
        matches_size[0] = 39877123421
        tags = numpy.zeros((2, 23), 'S1')
        tags[0] = stringtoarr('onetag', 23)
        tags[1] = stringtoarr('twotags', 23)
        matches_tags[0] = tags
        subset.createDimension('archive_unused', 2)
        subset.createVariable('archive_unused', 'S1', ['archive_unused', 'default_strlen'])
        subset.variables['archive_unused'][0] = stringtoarr('nothing', 23)
        subset.variables['archive_unused'][1] = stringtoarr('notmuch', 23)
        group1.createDimension('non_matching', 3)
        group1.createVariable('non_matching', 'S1', ['non_matching', 'default_strlen'])
        group1.variables['non_matching'][0] = stringtoarr('floop', 23)
        group1.variables['non_matching'][1] = stringtoarr('sloop', 23)
        group1.variables['non_matching'][2] = stringtoarr('kaput', 23)
        nc.close()

        # TEMP: show netcdf contents
        # print '\n' + subprocess.Popen(['ncdump', testfile.name], stdout=subprocess.PIPE).communicate()[0]

        # Read
        reader = CatalogueReaderNetCDF()
        result = reader.load(testfile.name)

        # Check results
        self.assertTrue(isinstance(result, Catalogue))
        self.assertEqual('citizen-of-the-world', result.identifier)
        self.assertTrue(isinstance(result.datasets, list))
        self.assertEqual(2, len(result.datasets))
        self.assertTrue(isinstance(result.datasets[0], CatalogueDataSet))
        self.assertEqual('Ffflip', result.datasets[0].name)
        self.assertEqual('/find/it/here', result.datasets[0].path)
        self.assertTrue(isinstance(result.datasets[0].name, basestring))
        self.assertTrue(isinstance(result.datasets[0].path, basestring))
        self.assertTrue(isinstance(result.datasets[1], CatalogueDataSet))
        self.assertEqual('Ffflop', result.datasets[1].name)
        self.assertEqual('/or/here', result.datasets[1].path)
        self.assertTrue(isinstance(result.datasets[1].name, basestring))
        self.assertTrue(isinstance(result.datasets[1].path, basestring))
        self.assertEqual(['floop', 'sloop', 'kaput'], result.datasets[1].non_matching)
        self.assertTrue(isinstance(result.datasets[1].non_matching[0], basestring))
        self.assertTrue(isinstance(result.datasets[1].non_matching[1], basestring))
        self.assertTrue(isinstance(result.datasets[1].non_matching[2], basestring))
        self.assertEqual(1, len(result.datasets[1].subsets))
        self.assertTrue(isinstance(result.datasets[1].subsets[0].layout, DataStorageFiles))
        self.assertEqual(1, isinstance(result.datasets[1].subsets[0].matches, list))
        self.assertEqual(1, len(result.datasets[1].subsets[0].matches))
        self.assertEqual('bob', result.datasets[1].subsets[0].matches[0].name)
        self.assertTrue(isinstance(result.datasets[1].subsets[0].matches[0].name, basestring))
        self.assertEqual(datetime(2016, 05, 30, 17, 29, 33), result.datasets[1].subsets[0].matches[0].time)
        self.assertEqual(39877123421, result.datasets[1].subsets[0].matches[0].size)
        self.assertEqual(2, len(result.datasets[1].subsets[0].matches[0].tags))
        self.assertEqual('onetag', result.datasets[1].subsets[0].matches[0].tags[0])
        self.assertEqual('twotags', result.datasets[1].subsets[0].matches[0].tags[1])
        self.assertEqual(['nothing', 'notmuch'], result.datasets[1].subsets[0].archive_unused)

 
class NamedTemporaryDirectory(object):
    """Temporary directory removed when instance goes out of scope."""

    def __init__(self):
        self.pathname = tempfile.mkdtemp()

    def __del__(self):
        shutil.rmtree(self.pathname)


class TestNamedTemporaryDirectory(unittest.TestCase):
    """Test for named temporary directory class (which is itself used primarily for testing)."""

    def test_create_delete(self):

        # Make the temp directory
        testdir = NamedTemporaryDirectory()

        # Write some file inside it
        filename = os.path.join(testdir.pathname, 'somefile')
        open(filename, 'w').write('sometext')

        pathname = testdir.pathname
        self.assertTrue(os.path.exists(pathname))
        self.assertTrue(os.path.isdir(pathname))
        self.assertTrue(os.access(filename, os.R_OK))
        del testdir
        self.assertFalse(os.path.exists(pathname))
        self.assertFalse(os.access(filename, os.R_OK))
