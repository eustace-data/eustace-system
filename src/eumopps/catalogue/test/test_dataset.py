"""Unit tests for dataset catalogue container classes."""

__version__ = "$Revision: 1019 $"
__author__ = "Joel R. Mitchelson"

import unittest
import os
from ..dataset import ChecksumThread
from ..dataset import CatalogueFileEntry
from ..dataset import CatalogueDataSubset

class TestCatalogueDataSubset(unittest.TestCase):
    """Test container class that holds information about subset and does checksum calculations."""

    def test_compute_checksum_thread_indices(self):
        """Check method to compute which indices to run on which threads for multi-threaded checksum."""

        # Normal case: number of matches is greater than number of threads and
        # doesn't exactly divide
        self.assertEqual([(0, 2), (2, 5)], CatalogueDataSubset.compute_checksum_thread_indices(
            total_matches=5, numthreads=2))

        # Normal case: number of matches greater than number of threads and
        # divides exactly
        self.assertEqual([(0, 2), (2, 4), (4, 6)], CatalogueDataSubset.compute_checksum_thread_indices(
            total_matches=6, numthreads=3))

        # Small data set: number of matches is less than number of threads
        # available
        self.assertEqual([(0, 1), (1, 2), (2, 3), (3, 4), (4, 5)],
                         CatalogueDataSubset.compute_checksum_thread_indices(total_matches=5, numthreads=64))

        # Edge case: number of matches exactly equals number of threads
        # available
        self.assertEqual([(0, 1), (1, 2), (2, 3)], CatalogueDataSubset.compute_checksum_thread_indices(
            total_matches=3, numthreads=3))

    def test_compute_checksum(self):
        """Test multi-threaded checksum."""

        # Just use same example file multiple times for test
        subset = CatalogueDataSubset(
            layout=None,
            matches=[
                CatalogueFileEntry(name='checksumexample.txt'),
                CatalogueFileEntry(name='nonsense'),
                CatalogueFileEntry(name='checksumexample.txt')])

        # Multi-threaded computation of checksum
        subset.compute_checksum(datasetpath=os.path.dirname(__file__))

        # Middle one should be None due to non-existent file
        self.assertEqual('3469475163', subset.matches[0].checksum)
        self.assertEqual(1431, subset.matches[0].size)
        self.assertEqual(None, subset.matches[1].checksum)
        self.assertEqual(None, subset.matches[1].size)
        self.assertEqual('3469475163', subset.matches[2].checksum)
        self.assertEqual(1431, subset.matches[2].size)


class TestChecksumThread(unittest.TestCase):
    """Tests for helper that does checksum on separate threads."""

    def test_constructor(self):
        """Ensure class constructs with correct member variables set."""

        # Run constructor
        example = ChecksumThread(datasetpath='/some/path', matches=[CatalogueFileEntry(
            'firstone'), CatalogueFileEntry('secondone'), CatalogueFileEntry('andthis')], indexstart=1, indexstop=3)

        # Check results
        self.assertEqual('/some/path', example.datasetpath)
        self.assertEqual(3, len(example.matches))
        self.assertEqual('firstone', example.matches[0].name)
        self.assertEqual('secondone', example.matches[1].name)
        self.assertEqual('andthis', example.matches[2].name)
        self.assertEqual(1, example.indexstart)
        self.assertEqual(3, example.indexstop)

    def test_run(self):
        """Run without actually spawning separate thread."""

        # Just use same example file multiple times for test
        matches = [
            CatalogueFileEntry('checksumexample.txt'),
            CatalogueFileEntry('checksumexample.txt'),
            CatalogueFileEntry('checksumexample.txt')]

        # Construct thread object
        example = ChecksumThread(
            datasetpath=os.path.dirname(__file__),
            matches=matches,
            indexstart=1,
            indexstop=3)

        # Should be flagged as not run initially
        self.assertFalse(example.iscomplete.is_set())

        # Now run
        example.run()

        # Should be flagged as run
        self.assertTrue(example.iscomplete.is_set())

        # Member variables of 2nd and 3rd array elements set accordingly
        # 1st array element should still be None
        self.assertEqual(None, matches[0].checksum)
        self.assertEqual(None, matches[0].size)
        self.assertEqual('3469475163', matches[1].checksum)
        self.assertEqual(1431, matches[1].size)
        self.assertEqual('3469475163', matches[2].checksum)
        self.assertEqual(1431, matches[2].size)

    def test_start_and_wait(self):
        """Run by creating child thread and waiting for result."""

        # Just use same example file multiple times for test
        matches = [
            CatalogueFileEntry('checksumexample.txt'),
            CatalogueFileEntry('checksumexample.txt'),
            CatalogueFileEntry('checksumexample.txt')]

        # Construct thread object
        example = ChecksumThread(
            datasetpath=os.path.dirname(__file__),
            matches=matches,
            indexstart=1,
            indexstop=3)

        # Now run and wait for answer
        example.start()
        example.iscomplete.wait()

        # Member variables of 2nd and 3rd array elements set accordingly
        # 1st array element should still be None
        self.assertEqual(None, matches[0].checksum)
        self.assertEqual(None, matches[0].size)
        self.assertEqual('3469475163', matches[1].checksum)
        self.assertEqual(1431, matches[1].size)
        self.assertEqual('3469475163', matches[2].checksum)
        self.assertEqual(1431, matches[2].size)

    def test_run_with_bad_file(self):
        """Check that run still completes with bad input."""

        # Just use same example file multiple times for test
        matches = [
            CatalogueFileEntry('checksumexample.txt'),
            CatalogueFileEntry('thisdoesnotexist'),
            CatalogueFileEntry('checksumexample.txt')]

        # Construct thread object
        example = ChecksumThread(
            datasetpath=os.path.dirname(__file__),
            matches=matches,
            indexstart=1,
            indexstop=3)

        # Should be flagged as not run initially
        self.assertFalse(example.iscomplete.is_set())

        # Now run
        example.run()

        # Should be flagged as run
        self.assertTrue(example.iscomplete.is_set())

        # Member variables of 2nd and 3rd array elements set accordingly
        # 1st array element should still be None
        self.assertEqual(None, matches[0].checksum)
        self.assertEqual(None, matches[0].size)
        self.assertEqual(None, matches[1].checksum)
        self.assertEqual(None, matches[1].size)
        self.assertEqual('3469475163', matches[2].checksum)
        self.assertEqual(1431, matches[2].size)
