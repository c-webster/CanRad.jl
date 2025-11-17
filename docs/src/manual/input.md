**variable name [relevant model]**

## Terrain Input
**hrdtmf [ L2R, C2R, T2R ]**\
high-resolution digital terrain model\
used to calculate the near horizon-line within 300m\
format: GeoTIFF

**lrdtmf [ L2R, C2R, T2R ]**\
low-resolution digital terrain model\
used to calculate the distant topographic horizon-line\
format: GeoTIFF

**termf [ L2R, C2R, T2R ]**\
terrain mask (calculated in CanRad if terrain_precalc = true)\
format: NetCDF

**hlmf [ L2R, C2R, T2R ]**\
horizon line matrix file of zenith angle of topography above the horizon\
(calculated externally, or in CanRad if save_horizon = true)\
format: NetCDF

## Canopy Input
**chmf [ C2R ]**\
canopy height model\
format: GeoTIFF

**lasf [ L2R ]**\
digital surface model\
lidar point cloud, ideally classified\
format: LAS/LAZ

**lavdf [ C2R ]**\
leaf-area volume density\
an output from canrad_precalc.jl\
format: GeoTIFF

**dbhf [ L2R, C2R ]**\
diameter at breast height file\
an output from canrad_precalc.jl\
contains dimensions for including trunks in the synthetic images (dbh)\
format: Plain Text as *_treeinfo.txt file

**ltcf [ L2R ]**\
lidar tree classifcation file
an output from pycrown and canrad_precalc.jl\
specifyies which points belong to the tree crowns in treeinfo\
format: LAZ

**cbhf [ C2R ]**\
canopy base height\
an output from canrad_precalc.jl\
used if crown base height varies over domain (e.g. 60% of total tree height)\
format: GeoTIFF at same resolution as chmf

## Meteorological Input
**swrf [ L2R, C2R, T2R ]**\
shortwave radiation file\
only needed if calc_swr = 2 in the settings file\
the model assumes the start and end date are the same as what you have specified in the settings file\
format: Plain Text tab-delimited files of three columns (date,time,value)
