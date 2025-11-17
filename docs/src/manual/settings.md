
**setting [relevant model]**

Settings with no default values are required in the settings file


## Terrain Settings

**terrain_highres [ L2R, C2R, T2R ]**\
*Boolean*\
*default = false*\
include high resolution (1-25 m) local terrain 

**highres_peri [ L2R, C2R, T2R ]**\
*Integer*\
*default = 300*\
radius (in metres) for including high-res terrain 
value in metres

**terrain_lowres [ L2R, C2R, T2R ]**\
*Boolean*\
*default = false (ground assumed flat)*\
include lower resolution (> 25 m) regional terrain 

**lowres_peri [ L2R, C2R, T2R ]**\
*integer*\
*default = 25000*\
radius (in metres) for low-res terrain (calculating regional horizon line)

**terrainmask_precalc [ L2R, C2R ]**\
*Boolean*\
*default = false*\
used only if the terrain mask has been pre-calculated and saved by T2R (can be useful when calibrating L2R and needing to run multiple times)

**calc_terrain [ L2R, C2R ]**\
*Boolean*\
*default = false*\
calculate variables for terrain-only as well as forest

**buildings [ L2R, C2R, T2R ]**\
*Boolean*\
*default = false*\
include buildings in terrain calculation - requires blmf file\
buildings within the same radius as highres_peri are included in the calculations

**hlm_precalc [ L2R, C2R ]**\
*Boolean*\
*default = false*\
use pre-calculated horizon line matrix (e.g. for compatibility with topo-downscaled swr data)

## Forest settings [ L2R, C2R ]

**keep_points [ L2R ]**\
*String*\
options:\
"canopy" (*default*) : only loads lidar points classified as vegetation (classes 3,4,5)\
"ground+canopy" : loads ground and canopy classified points (classes 2-5)\
"all" : loads all points, including those unclassified and those classified as noise

**forest_peri  [ L2R, C2R ]**\
*Integer*\
*default = 100*\
radius within which forest is included (note all model points must be a minimum of this distance from the edge of the forest canopy datasets)

**forest_type [ C2R ]**\
*String*\
options:\
"evergreen" (*default*) : 100% evergreen forest; no phenology in calculations (note: only developed/tested for coniferous evergreen species)\
"deciduous"* : 100% deciduous forest\
"mixed"* : a mix of deciduous and evergreen species (only available with one of the "special_implementation" options \
*Requires *phenology* to also be specified\

**phenology [ L2R, C2R ]**\
*String*\
*default = "none"*\
In C2R:\
&emsp;"leafon" : calculates for leaf-on/summertime canopy conditions \
&emsp;"leafoff" : calculates for leaf-off/wintertime canopy conditions  \
&emsp;"both" : calculates output data for both leaf-on and leaf-off conditions \
&emsp;"none" (*default*) : phenology is not accounted for and evergreen forest_forest is assumed \
In L2R:\
&emsp;"leafoff" : calculates for leaf-off/wintertime canopy conditions (calculates trunk and canopy points differently to reduce canopy density in winter)  \
&emsp;"leafon"* : doesn't change calculations but used to name variables in output files
&emsp;"none"* (*default*) :  doesn't change calculations but used to name variables in output files (i.e forest_type = "evergreen")
*uses only user input *point_size* to create hemispheric image

**leaf_type [ C2R ]** \
<formerly tree_species>
*String*\
*default = "needleleaf"*\
"needleleaf" or "broadleaf"\
Used in \
a) the generation of the hemispheric images to distinguish between the different crown shapes of needleleaf and broadleaf forests\
b) distinguishing between larch (deciduous needleleaf) and decidious broadleaf forests in *special_implementation*

**trunks [ L2R, C2R ]**\
*Boolean*\
*default = false*\
adds opaque trunks to the hemispheric images - requires that the canrad_precalc.jl steps are run to create additional datasets. 

**trunks_peri [ L2R, C2R ]**\
*Integer*\
*default = forest_peri*\
radius within which trunks are included - in dense canopies, this can be reduced to lower than forest_peri to save computation time. 


**branches [ L2R ]**\
*Boolean*\
*default = false*\
adds additional branch points to the image to increase canopy density\
method: adds additional points between existing points and the centre of the tree crown - requires that the canrad_precalc.jl steps are run to create additional datasets. 

**branch_spacing [ L2R ]**\
*Float*\
*default = 0.1*\
spacing interval (m) for creating branches\

**cbh [ C2R ]**\
*Float*\
*default = 2.0*\
height of canopy base above ground (in m)\


## Calculation settings

**calc_trans [ L2R, C2R, T2R ]**\
*Boolean*\
*default = false*\
false = only sky-view fraction is calculated and saved\
true  = time-varying direct-beam transmissivity is calculated and saved


**calc_swr [ L2R, C2R, T2R ]**\
*Integer*\
*default = 0*\
0  = shortwave radiation not calculated \
1* = calculates potential swr using atm_trans = 1 \
2* = calculates true swr, requiring measured or modelled above-canopy swr data suppled in the file 'swrf' in dat_in\
*requires that calc_trans=true


## Time settings

**t1 and t2 [ L2R, C2R, T2R ]**\
*String* \
*default =  "22.06.2020 00:00:00" and "22.06.2020 23:00:00"*\
start and end time steps for calculation\
format: "dd.mm.yyyy HH:MM:SS"

**tstep [ L2R, C2R, T2R ]**\
*Integer*\
*default = 60*\
time step in minutes for the output transmissivity and shortwave radiation data \
all data is calculated at a 2 minute interval and aggregated to tstep to saving, therefore the minimum value for this setting is 2 \



## Location settings

**time_zone [ L2R, C2R, T2R ]**\
*Integer*\
*required if calc_trans=true*\
local time zone relative to UTC (e.g. CET = UTC+1, time_zone = 1)


**epsg_code [ L2R, C2R, T2R ]**\
*Integer*\
*required if calc_trans=true*\
relevant epsg code for the coordinate system of the input and output datasets\
note: all input datasets must be in this coordinate system\


## Image settings

**image_height [ L2R, C2R, T2R ]**\
*Float*\
*default = 0.5*\
height (m) of virtual camera position above ground

**point_size [ L2R ]**\
*2-element Vector{Float64}*\
*default = [ 30.0 , 10.0 ]*\
size of points in synthetic hemispheric images
First number is the size of points closest to the camera, second number is point size futherest away from camera\
Values are in pixels and image is 10000x10000 pixels


## Run settings

**batch [ L2R, C2R, T2R ]**\
*boolean*\
*default = false*\
used if running multiple tiles to manage sub-folders in the output folder structure

**save_images [ L2R, C2R, T2R ]**\
*boolean*\
*default = false*\
save the calculated images in a netcdf file 

**make_pngs [ L2R, C2R, T2R ]**\
*boolean*\
*default = false*\
generates .png files from the saved images

**save_horizon [ C2R, T2R ]**\
*boolean*\
*default = false*\
saves the calculated topographic horizon line (terrain only)

**progress [ L2R, C2R, T2R ]**\
*boolean*\
*default = false*\
reports individual step progress in /ProgressLastPoint\
(good for checking run-time in L2R, both otherwise doesn't need to be used)


## Special implementations [ C2R ]

Written for a specific application and uses specific combinations of datasets

Current options:\
**swissrad** used to calculate the Swiss nationwide dataset \
**oshd** used to calculate across domains relevant to oshd (uses a pre-calculated horizon line to ensure compatibility with downscaled SWR) \
**oshd-alps** used to calculate across the European Alps\
**none** used for all other C2R use cases and requires lavdf in dat_in

