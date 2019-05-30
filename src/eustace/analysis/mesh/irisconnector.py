"""Express mesh as iris cube suitable for storing."""

import numpy
from iris.coords import DimCoord
from iris.cube import Cube
from iris.cube import CubeList
import netCDF4
from geometry import cartesian_to_polar2d
from geometry import polar2d_to_cartesian
from mesh import Mesh

class MeshToIrisCubeList(CubeList):

    def __init__(self, mesh):

        # Call iris base constructor
        super(MeshToIrisCubeList, self).__init__()

        # Mesh identifier
        mesh2identifier = Cube(
            numpy.ma.array(data=numpy.int32(-1), mask=True, fill_value=numpy.int32(-1)),
            var_name='Mesh2',
            long_name='Topology data of 2D unstructured mesh',
            attributes = {
                'topology_dimension': numpy.int32(2),
                'node_coordinates': 'Mesh2_node_x Mesh2_node_y',
                'face_node_connectivity': 'Mesh2_face_nodes',
                'edge_node_connectivity': 'Mesh2_edge_nodes',
                'EUSTACE_subdivison_level': numpy.int32(mesh.level),
                }
            )

        # Convert data
        polarcoords = cartesian_to_polar2d(mesh.points)

        longitudes = Cube(
            polarcoords[:,1],
            var_name='Mesh2_node_x',
            # standard_name='node_longitude',
            long_name='Longitude of 2D mesh nodes.',
            units='degrees_east')

        latitudes = Cube(
            polarcoords[:,0],
            var_name='Mesh2_node_y',
            # standard_name='node_latitude',
            long_name='Latitude of 2D mesh nodes.',
            units='degrees_north')

        # Triangle numbers - a EUSTACE-specific concept
        trianglenumbers = Cube(
            numpy.array(mesh.triangles[:,0], numpy.uint64),
            var_name='Mesh2_EUSTACE_trianglenumbers')

        # Indices of faces
        face_nodes = Cube(
            numpy.array(mesh.triangles[:,1:], numpy.int32),
            var_name='Mesh2_face_nodes',
            long_name = 'Maps every triangle face to its corner nodes.',
            attributes={'cf_role': 'face_node_connectivity', 'start_index': 0})

        edge_nodes = Cube(
            numpy.array(mesh.edges[:,0:2], numpy.int32),
            var_name='Mesh2_edge_nodes',
            long_name = 'Maps every edge to the two nodes that it connects.',
            attributes={'cf_role': 'edge_node_connectivity', 'start_index': 0})

        edge_detail = Cube(
            numpy.array(mesh.edges[:,2:], numpy.uint64),
            var_name='Mesh2_EUSTACE_edge_triangle_map',
            long_name = 'Four numbers for each edge containing triangle 0 id, triangle 0 edge number, triangle 1 id, triangle 1 edge number.')
        
        self.append(mesh2identifier)
        self.append(longitudes)
        self.append(latitudes)
        self.append(trianglenumbers)
        self.append(face_nodes)
        self.append(edge_nodes)
        self.append(edge_detail)

class IrisCubeListToMesh(Mesh):

    def __init__(self, cubelist):

        # retrieve locations
        latitudes = next(field for field in cubelist if field.var_name == 'Mesh2_node_y')
        longitudes = next(field for field in cubelist if field.var_name == 'Mesh2_node_x')

        # convert to cartesian
        cartesian = polar2d_to_cartesian( numpy.vstack((latitudes.data, longitudes.data)).T )

        # load EUSTACE-specific triangle number info
        trianglenumbers = next(field for field in cubelist if field.var_name == 'Mesh2_EUSTACE_trianglenumbers').data

        # load face info
        faces = numpy.array(next(field for field in cubelist if field.var_name == 'Mesh2_face_nodes').data, numpy.uint64)

        # insert first column of triangle numbers
        triangles = numpy.vstack((trianglenumbers.T, faces.T)).T

        # load edge info
        faces = numpy.array(next(field for field in cubelist if field.var_name == 'Mesh2_face_nodes').data, numpy.uint64)

        # rebuild edge information
        edge_nodes = numpy.array(next(field for field in cubelist if field.var_name == 'Mesh2_edge_nodes').data, numpy.uint64)
        edge_triangle_map = next(field for field in cubelist if field.var_name == 'Mesh2_EUSTACE_edge_triangle_map').data
        edges = numpy.vstack((edge_nodes.T, edge_triangle_map.T)).T

        # build mesh object
        super(IrisCubeListToMesh, self).__init__(points=cartesian, triangles=triangles, edges=edges)
