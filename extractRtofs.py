## This script must be run inside a conda environment. 
## RTOFS is a high resolution current forecast model developed jointly by the US Navy and NOAA.
## The model covers the entire world, but is only available in a particularly nasty NetCDF configuration.
## All this script does is look locally for any .nc files, and it assumes any .nc files found are 
## RTOFS files. This script then 1 by 1 reads in each file,
## and writes it to two NEW NetCDF files (one for north/south direction, one for east/west direction),
## that can be understood by other weather utilities for later conversion to Grib. 

from numpy.testing._private.utils import IgnoreException

# These libraries are why we have to do this in Conda Python instead of Node.
import iris
import iris_grib
import cartopy
import numpy as np
import numpy.ma as ma
import os
import datetime

# For timing how long this takes. On average it takes about 1.33 hours to do 24 hours of RTOFS forecasts. Not great.
start = datetime.datetime.now()

files = []
for file in os.listdir("./"):
    if file.endswith(".nc"):
        files.append(file)


def formatCube(cube, new_name):
    
    cube.remove_coord(cube.coord(var_name='Layer'))
    # These lines should be commented out for the daily file, and should be restored for the hourly files. (RTOFS gives a 24 hour file)
    cube.remove_coord(cube.coord(var_name='Y'))
    cube.remove_coord(cube.coord(var_name='X'))

    # In order to convert to grib, we need to add some dimensions.
    cube.add_aux_coord(iris.coords.DimCoord(0, standard_name='forecast_period', units='hours'))
    cube.add_aux_coord(iris.coords.DimCoord(0, "height", units="m"))

    # These are magic values from Grib documentation.
    cube.rename(new_name)
    cube.units = 'm s-1'
    
    # TODO: Can we get away with spliting the file by the two projections and just assigning this?
    WGS_84_CS = iris.coord_systems.GeogCS(semi_major_axis=6378137,
    inverse_flattening=298.257223563)

    new_cube, extent = iris.analysis.cartography.project(cube, WGS_84_CS)

    new_cube.remove_coord('latitude')
    new_cube.remove_coord('longitude')
    new_cube.coord(standard_name='projection_y_coordinate').rename('latitude')
    new_cube.coord(standard_name='projection_x_coordinate').rename('longitude')
    new_cube.coord(standard_name='latitude').units = 'degrees'
    new_cube.coord(standard_name='longitude').units = 'degrees'
 
    return new_cube

## RTOFS NetCDF contains 2 projections in 1 file. Above 47 is an Arctic Bi-Polar, below is Mercator. 
# [leftLon, lowerLat, rightLon, upperLat]
arctic_bbox = [-180., 47., 180., 90.]
mercator_bbox = [-180., -90., 180., 47.]

minmax = lambda x: (np.min(x), np.max(x))

def bbox_extract_2Dcoords(cube, bbox):
    """
    Extract a sub-set of a cube inside a lon, lat bounding box
    bbox=[lon_min lon_max lat_min lat_max].
    NOTE: This is a work around to subset an iris cube that has
    2D lon, lat coords.
    """
    lons = cube.coord('longitude').points
    lats = cube.coord('latitude').points

    lons_inregion = np.logical_and(lons > bbox[0], lons < bbox[2])
    lats_inregion = np.logical_and(lats > bbox[1], lats < bbox[3])
    inregion = np.logical_and(lons_inregion, lats_inregion)

    region_inds = np.where(inregion)
    imin, imax = minmax(region_inds[0])
    jmin, jmax = minmax(region_inds[1])
    return cube[..., imin:imax+1, jmin:jmax+1]

# Leave this because we may want to try to split the files so we don't need to re-project
# east_arctic = bbox_extract_2Dcoords(new_east, arctic_bbox)
# east_mercator = bbox_extract_2Dcoords(new_east, mercator_bbox)

# north_arctic = bbox_extract_2Dcoords(new_north, arctic_bbox)
# north_mercator = bbox_extract_2Dcoords(new_north, mercator_bbox)

# cube.aux_coords[1].coord_system=iris.coord_systems.GeogCS(654321)
# cube.aux_coords[2].coord_system=iris.coord_systems.GeogCS(654321)

num_files = 0

for file in files:
    file_start = datetime.datetime.now()
    print('Starting file: ')
    print(file)
    cubes = iris.load(file)
    east_cube = cubes.extract("eastward_sea_water_velocity")[0]
    east_cube.attributes = {'GRIB_PARAM': 'GRIB2:d010c001n002'}
    north_cube = cubes.extract("northward_sea_water_velocity")[0]
    north_cube.attributes = {'GRIB_PARAM': 'GRIB2:d010c001n003'}

    # Magic variables for CDO to understand what we have in the NetCDF.
    eastward = formatCube(east_cube, "ucurr")
    northward = formatCube(north_cube, "vcurr")

    iris.save(northward, file.split('.nc')[0] + '_northward.nc', fill_value=0)
    iris.save(eastward, file.split('.nc')[0] + '_eastward.nc', fill_value=0)

    file_end = datetime.datetime.now()
    file_delta = file_end - file_start
    print('File took seconds: ')
    print(file_delta.seconds)
    num_files += 1

finish = datetime.datetime.now()
delta = finish - start
print('Seconds: ')
print(delta.seconds)

print('File average process time in seconds: ')
print(delta.seconds / num_files)


