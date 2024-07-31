
**setting [relevant model]**

Settings with no default values are required in the settings file


## Terrain Settings

**terrain_highres [ L2R, C2R, T2R ]**\
*boolean*
include high high resolution local terrain (1-5 m)

**highres_peri [ L2R, C2R, T2R ]**\
*integer*\
radius for including high-res terrain 
value in metres
default = 300


**terrain_lowres [ L2R, C2R, T2R ]**\
*boolean*
include low resolution regional terrain (> 25 m)

**lowres_peri [ L2R, C2R, T2R ]**\
*integer*\
radius for low-res terrain (calculating regional horizon line)
default = 25000

**terrainmask_precalc [ L2R, C2R ]**\
*boolean*\
if the terrain mask has been pre-calculated by CanRad-T2R (can be useful when calibrating L2R and needing to run multiple times)
default = false

**calc_terrain [ L2R, C2R ]**\
*boolean*\
default = false
calculate variables for terrain-only as well as forest

**buildings [ L2R, C2R, T2R ]**\
*boolean*\
default = false
include buildings in terrain calculation - requires blmf file

**hlm_precalc [ L2R, C2R ]**\
*boolean*\
default = false
use pre-calculated horizon line matrix (e.g. for compatibility with topo-downscaled swr data)


## Forest settings [ L2R, C2R ]

**keep_points [ L2R ]**
*String*
options:
"canopy" (default) : only loads lidar points classified as vegetation (classes 3,4,5)
"ground+canopy" : loads ground and canopy classified points (classes 2-5)
"all" : loads all points, including those unclassified and those classified as noise

**forest_peri  [ L2R, C2R ]**
*integer*
radius within which forest is included (note all model points must be a minimum of this distance from the edge of the forest canopy datasets)


**forest_type [ C2R ]**\
*String*\
options:
"evergreen" 100% evergreen forest; no seasonality in calculations. note: only developed/tested for coniferous evergreen species.\
"deciduous" 100% deciduous forest. Requires *seasonality* to also be specified.\
"mixed" forest is a mix of deciduous and evergreen species. R


**season [ L2R, C2R ]**
*String*

In C2R:
Options are **summer**, **winter**, **both**\
Model calculates output for summer (leaf-on) and/or winter (leaf-off) scenarios. 
Only applies if "deciduous" or "mixed" are indicated in the forest type setting.

In L2R:
Options are **winter**, **none**
Used in generation of the hemispheric images to represent trunk and canopy points differently under winter canopy conditions.

**tree_species [ C2R ]**\
*String*\
**needleleaf** or **broadleaf**\
Used in 
1. the generation of the hemispheric images to distinguish between the different crown shapes of needleleaf and broadleaf forests
2. distinguishing between larch (deciduous needleleaf) and decidious broadleaf forests.  


**trunks [ L2R ]**
*boolean*\
default = false

**branches [ L2R ]**
*boolean*\
default = false

**branch_spacing [ L2R ]**
*Float*
spacing interval for creating branches/
units = metres/
default = 0.1/

**cbh [ C2R ]**
*Float*
height of canopy base above ground in m across field area/
defual = 2/

## Calculation settings

**calc_trans**\
*boolean*

**calc_swr**\
*boolean*
0 = off; 1 = potential swr (atm_trans = 1); 2 = real swr (needs swrf in dat_in)


## Time step settings

**t1** and **t2** \start and end time steps for calculation
*String* 
format: "dd.mm.yyyy HH:MM:SS"

**tstep** \time step in minutes for the output data 


## Location settings

**time_zone** relative to UTC (e.g. CET = UTC+1, time_zone = 1)
*integer*

**coor_system** (UTM)
*String*
e.g. "UTM 32 N" (spaces required)\
note: all input datasets must be in this coordinate system\
other possible options are "CH1903", "CH1903+" for Swiss datasets


## Image settings

**image_height** height of virtual camera position above ground
*Float*

**point_size [ L2R ]**  size of points in hemispheric images
*2-element Vector{Int64}*
First number is point size are closest to the camera, second number is point size futherest away from camera
Values are in pixels (image is 10000x10000 pixels)


## Run settings

**batch**
*boolean*


**save_images**
*boolean*
default = false
save the calculated images in a netcdf file and generates .png files

**save_horizon [ C2R, T2R ]**\
*boolean*\
default = false

**progress**
*boolean*
default = false
reports individual step progress in /ProgressLastPoint


## Special implementations [ C2R ]

Written for a specific application and uses specific combinations of datasets

Current options:\
**swissrad** used to calculate the Swiss nationwide dataset. 
**oshd** used to calculate across domains relevant to oshd (uses a pre-calculated horizon line to ensure compatibility with downscaled SWR)


