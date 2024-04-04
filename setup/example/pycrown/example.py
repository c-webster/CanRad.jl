"""
PyCrown - Fast raster-based individual tree segmentation for LiDAR data
-----------------------------------------------------------------------
Copyright: 2018, Jan Zörner
Licence: GNU GPLv3
"""

from datetime import datetime

from pycrown import PyCrown

import sys

if __name__ == '__main__':

    TSTART = datetime.now()

    F_CHM = sys.argv[1]+'/data/CHM.tif'
    F_DTM = sys.argv[1]+'/data/DTM.tif'
    F_DSM = sys.argv[1]+'/data/DSM.tif'
    F_LAS = sys.argv[1]+'/data/POINTS.laz'

    PC = PyCrown(F_CHM, F_DTM, F_DSM, F_LAS, outpath=sys.argv[1]+'/result')

    # Smooth CHM with 5m median filter
    PC.filter_chm(3, ws_in_pixels=False)
    
    # Tree Detection with local maximum filter
    PC.tree_detection(PC.chm, ws=3, ws_in_pixels=True, hmin=8.)

    # Crown Delineation
	# ['dalponte_cython', 'dalponte_numba',
    # 	'dalponteCIRC_numba', 'watershed_skimage']
    PC.crown_delineation(algorithm='watershed_skimage', th_tree=8.,
                         th_seed=0.7, th_crown=0.55, max_crown=10.)

    # Correct tree tops on steep terrain
    PC.correct_tree_tops()

    # Calculate tree height and elevation
    PC.get_tree_height_elevation(loc='top')
    PC.get_tree_height_elevation(loc='top_cor')

    # Screen small trees
    PC.screen_small_trees(hmin=5., loc='top')

    # Convert raster crowns to polygons
    PC.crowns_to_polys_raster()
    PC.crowns_to_polys_smooth(store_las=True)

    # Check that all geometries are valid
    PC.quality_control()

    # Export results
    PC.export_tree_locations(loc='top')
    PC.export_tree_locations(loc='top_cor')
    PC.export_tree_crowns(crowntype='crown_poly_raster')
    PC.export_raster(PC.chm, PC.outpath / 'chm.tif', 'CHM')

    # PC.export_tree_crowns(crowntype='crown_poly_smooth')

    TEND = datetime.now()

    print(f"Number of trees detected: {len(PC.trees)}")
    print(f'Processing time: {TEND-TSTART} [HH:MM:SS]')
