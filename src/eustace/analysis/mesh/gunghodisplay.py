"""Display of cube using modified GungHo code."""

import cartopy.crs as ccrs
import iris
import matplotlib.colors as mcolors
import matplotlib.cm
import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry.polygon import Polygon
import argparse

def main():

    parser = argparse.ArgumentParser('eustace.analysis.mesh.gunghodisplay')
    parser.add_argument('filepath')
    args = parser.parse_args()

    # Get path to NetCDF file from command line
    filepath = args.filepath

    # Load in the cube that defines the face-node connectivity.
    face_node = iris.load_cube(filepath, 'Maps every triangle face to its corner nodes.')

    # Load in cubes of lat-lon locations of the nodes.
    lats = iris.load_cube(filepath, 'Latitude of 2D mesh nodes.')
    lons = iris.load_cube(filepath, 'Longitude of 2D mesh nodes.')

    # To avoid problem with wrapping longitudes at GMT, use longitude
    # range -180 to 180.
    lons.data[lons.data>180] -= 360

    # Mock up some data for each face.
    data = np.random.rand(face_node.shape[0])

    # Create a colourmap.
    norm = mcolors.Normalize(vmin=data.min(), vmax=data.max())

    cmap = matplotlib.cm.get_cmap('plasma')

    projections = [ccrs.Orthographic(), ccrs.PlateCarree(), ccrs.Robinson()]

    for projection in projections:

        fig, ax = plt.subplots(figsize=(19, 11), subplot_kw=dict(projection=projection))

        # For each face create a polygon and assign it a colour.
        for face_index, face in enumerate(face_node.data):

            # Determine the lat-lon locations of the node that makes up
            # each face.
            locs = [ (lons.data[index], lats.data[index]) for index in reversed(face) ]

            # Add start node to the end to complete the loop.
            locs.append(locs[0])

            # Create polygon to represent the face
            poly = Polygon(locs)

            # Get data value at that face.
            data_val = data[face_index]

            # Choose colour of polygon based on the data value.
            colour = cmap(norm(data_val))

            ax.add_geometries([poly], ccrs.Geodetic(), facecolor=colour)
        
    ax.coastlines()
    plt.show()

# Call main entry point
if __name__ == "__main__":
    main()
