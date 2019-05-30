
import unittest
from ..irisconnector import MeshToIrisCubeList
from ..irisconnector import IrisCubeListToMesh
from ..mesh import MeshIcosahedronSubdivision
import iris
import numpy

class TestIrisConnector(unittest.TestCase):

    def test_round_trip(self):

        testmesh = MeshIcosahedronSubdivision.build(0)
        cubelist = MeshToIrisCubeList(testmesh)
        result = IrisCubeListToMesh(cubelist)
        numpy.testing.assert_almost_equal(result.points, testmesh.points)
        numpy.testing.assert_equal(result.triangles, testmesh.triangles)
        numpy.testing.assert_equal(result.edges, testmesh.edges)
