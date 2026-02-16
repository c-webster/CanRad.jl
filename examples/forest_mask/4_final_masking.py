from shapely.geometry import shape
from shapely.ops import unary_union

"""

Part of workflow to calculate forest mask for running CanRad-c2r over large domains

STEP 4 (optional):

Takes the buffered binary forest mask at canrad resolution:
    1. (optional) masks large water bodies to 0
    2. (optional) clips raster to domain polygon 

Uses multi-processing to speed up processing of large datasets by breaking the input raster into blocks

- input: binary buffered forest mask and different resolutions
- output: final forest masks with 0 = nodata, 1 = forest, 2 = open 

NOTE: unlike the other files, each individual mask must be run individually


# .gitignore: note to CW: run on a SC-canrad machine in the environment gis-env in ~/windows/projects/eurad/eu_data


"""


import os
import math
import rasterio
from rasterio.features import geometry_mask
from affine import Affine
import fiona
from rasterio.mask import mask
import numpy as np

inraster_fn = "forest_data/ForestMask/calculations/eualps_ForestMask_010m_40mbuffer_agg0100m.tif"
outstr = "oshd_alps"
outraster_fn = f"forest_data/ForestMask/{outstr}_ForestMask_100m.tif"


mask_land = True
land_mask_poly = "shapefiles/land_polygon/land_polygon.shp"

clip_to_domain = True
# domain_polygon = "../oshd_alps/shapefiles/oshd_alps_domain_03035_polygon.shp"
domain_polygon = "shapefiles/alps_polygon/alps_polygon_03035.shp"


# ---------------
# Main
# ---------------

with rasterio.open(inraster_fn) as src:
    profile = src.profile
    data = src.read(1).astype(np.int8)

    # Set NaN values to 2
    # data[np.isnan(data)] = 2
    # data[data == 0] = 2  # Assuming original open values are 0

    if mask_land:
        # Mask large water bodies to 0 using land polygon
        with fiona.open(land_mask_poly, "r") as shapefile:
            shapes = [feature["geometry"] for feature in shapefile]
            land_mask = geometry_mask(
                shapes,
                transform=src.transform,
                invert=True,
                out_shape=(src.height, src.width)
            )

        data[~land_mask] = -1
    
    if clip_to_domain:
        # Clip raster to domain polygon - manually crop to preserve existing data modifications
        with fiona.open(domain_polygon, "r") as shapefile:
            geoms = [shape(feature["geometry"]) for feature in shapefile]
            domain_geom = unary_union(geoms)
            bounds = domain_geom.bounds  # (minx, miny, maxx, maxy)
            
            # Round bounds to nearest km
            x_min_rounded = math.floor(bounds[0] / 1000) * 1000
            y_min_rounded = math.floor(bounds[1] / 1000) * 1000
            x_max_rounded = math.ceil(bounds[2] / 1000) * 1000
            y_max_rounded = math.ceil(bounds[3] / 1000) * 1000
            
            # Convert bounds to pixel coordinates
            inv_transform = ~src.transform
            col_min, row_max = inv_transform * (x_min_rounded, y_min_rounded)
            col_max, row_min = inv_transform * (x_max_rounded, y_max_rounded)
            
            col_min = int(math.floor(col_min))
            col_max = int(math.ceil(col_max))
            row_min = int(math.floor(row_min))
            row_max = int(math.ceil(row_max))
            
            # Clip data array
            data = data[row_min:row_max, col_min:col_max]
            
            # Create new transform for clipped extent
            out_transform = Affine(
                src.transform.a, src.transform.b, x_min_rounded,
                src.transform.d, src.transform.e, y_max_rounded
            )
            profile.update(transform=out_transform, height=data.shape[0], width=data.shape[1])
            
            # Mask pixels outside polygon geometry to -1
            domain_mask = geometry_mask(
                geoms,
                transform=out_transform,
                invert=True,
                out_shape=(data.shape[0], data.shape[1])
            )
            data[~domain_mask] = -1

    # Write output raster
    profile.update(dtype=rasterio.int8, nodata=-1)

    with rasterio.open(outraster_fn, 'w', **profile) as dst:
        dst.write(data.astype(rasterio.int8), 1)