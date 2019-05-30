"""Save mesh as NetCDF."""

from mesh import MeshIcosahedronSubdivision
from irisconnector import MeshToIrisCubeList
from eustace.outputformats.iriscubesaver import NetCDFAttributesEUSTACE
import iris
import argparse

def main():

    # Parse command line to get file path
    parser = argparse.ArgumentParser('eustace.analysis.mesh.writenetcdf')
    parser.add_argument('filename', help='Output filename')
    parser.add_argument('--subdivision', type=int, default=0, help='Levels of subdivision')
    parser.add_argument('--face_centre_at_poles', type=bool, default=True, help='Rotate to put face centre at poles')
    args = parser.parse_args()

    # Make mesh
    mesh = MeshIcosahedronSubdivision.build(args.subdivision, face_centre_at_poles=args.face_centre_at_poles)

    # Convert to iris object
    cubelist = MeshToIrisCubeList(mesh)

    # apply attributes so that they appear in final NetCDF
    NetCDFAttributesEUSTACE(title='EUSTACE Mesh', institution='Met Office', source=MeshIcosahedronSubdivision.__module__, comment='').apply(cubelist)

    # Write
    iris.FUTURE.netcdf_no_unlimited = True
    iris.save(cubelist, args.filename)


# Call main entry point
if __name__ == "__main__":
    main()
