# Configuration Settings 
 
This page provides a complete reference for all configuration settings in CanRad. Settings are defined in Julia script files (e.g., `C2R_settings.jl`, `L2R_settings.jl`, `T2R_settings.jl`). 
 
## Overview 
 
**Required vs Optional Settings:** 
- All configuration settings have default values, but the default configuration is designed for the code to run without error rather than produce a default result for an example forest.  
- Some settings are required only if certain features are enabled (e.g., `t1`, `t2`, and `time_zone` are required if `calc_trans = true`). 
 
**Model Applicability:** 
- Each setting indicates which models it applies to: [L2R], [C2R], [T2R] 
- You need only include settings relevant to your chosen model (settings not relevant to your chosen model will be ignored) 
 
 
## Terrain Settings 
 
Configure topographic shading effects from local and regional terrain. 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|-------------| 
| `terrain_highres` | Boolean | `false` | L2R, C2R, T2R | Include high-resolution (1-25 m) local terrain for near-field shading calculations | 
| `highres_peri` | Integer | `300` | L2R, C2R, T2R | Radius (meters) for including high-resolution terrain around each calculation point | 
| `terrain_lowres` | Boolean | `false` | L2R, C2R, T2R | Include low-resolution (>25 m) regional terrain for distant horizon calculations. If false, ground is assumed flat. | 
| `lowres_peri` | Integer | `25000` | L2R, C2R, T2R | Radius (meters) for low-resolution terrain (calculating regional horizon line). Default = 25 km. | 
| `terrainmask_precalc` | Boolean | `false` | L2R, C2R | Use pre-calculated terrain mask (useful for calibration runs or when reusing T2R outputs) | 
| `calc_terrain` | Boolean | `false` | L2R, C2R | Also calculate and save terrain-only (above-canopy) variables for comparison | 
| `buildings` | Boolean | `false` | L2R, C2R, T2R | Include buildings in terrain shading calculations (requires `blmf` input file). Buildings within `highres_peri` are included. | 
| `hlm_precalc` | Boolean | `false` | L2R, C2R | Use pre-calculated horizon line matrix instead of computing from terrain (ensures compatibility with topographically-downscaled radiation data) | 
 
**Typical configurations:** 
- **Flat terrain:** Keep both `terrain_highres` and `terrain_lowres` as `false` 
- **Forest gaps:** Set `terrain_highres = true`, `highres_peri = 300` 
- **Mountain valleys:** Set both to `true`, `lowres_peri = 25000` or larger 
 
<br> 
 
## Forest settings [ L2R, C2R ] 
 
Configure how forest canopy is represented in the model. 
 
### Spatial Extent Settings 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|-------------| 
| `forest_peri` | Integer | `100` | L2R, C2R | Radius (meters) within which forest is included in calculations. All calculation points must be at least this distance inside the extent of all forest canopy datasets. | 
| `trunks_peri` | Integer | `forest_peri` | L2R, C2R | Radius (meters) within which trunks are included. Can be reduced below `forest_peri` in dense canopies to save computation time. | 
 
**Important:** Ensure input data extends at least `forest_peri` beyond all calculation points to avoid edge effects. 
 
 
 
 
### Forest Type and Phenology 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|   ----------| 
| `forest_type` | String | `"evergreen"` | C2R |  Forest type of the model domain | 
 
 
> | Option | Description | 
> |--------|-------------| 
> | `"evergreen"` | Coniferous or evergreen forest (i.e. no phenology); currently only developed/tested for coniferous evergreen forests | 
> | `"deciduous"` | Deciduous forest: requires `phenology` setting | 
> | `"mixed"` | Mixed forest: requires `special_implementation` and dataset of forest type | 
 
### Phenology 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|-------------| 
| `phenology` | String | `"none"` | L2R, C2R | Seasonal canopy state (see table below) | 
 
**Phenology Options for C2R:** 
 
> | Options | Description | 
> |-------|-------------| 
> | `"none"` | No phenology; assumes evergreen forest | 
> | `"leafon"` | Leaf-on (summer) conditions only | 
> | `"leafoff"` | Leaf-off (winter) conditions only | 
> | `"both"` | Calculates both leaf-on and leaf-off outputs | 
 
**Phenology Options for L2R:** 
> | Options | Description | 
> |-------|-------------| 
> | `"none"` | No phenology adjustments (evergreen) | 
> | `"leafoff"` | Reduces canopy density for winter; treats trunks and canopy separately | 
> | `"leafon"` | No calculation changes; used for output variable naming | 
 
### Leaf Characteristics (C2R only) 
 
| Setting | Type | Default | Description | 
|---------|------|---------|-------------| 
| `leaf_type` | String | `"needleleaf"` | Leaf morphology.<br>Used to: (a) distinguish crown shapes in image generation, (b) differentiate larch from broadleaf deciduous species in special implementations. | 
 
> | Options | Description | 
> |---------|-------------| 
> | `"needleleaf"` | typically used for coniferous evergreen forests and larch deciduous forests; assumes conical tree crown shape 
> | `"broadleaf"` | typically used for mid-latitude/altitude deciduous broadleaf forests; assumes spherical tree crown shape 
 
### Point Cloud Settings (L2R only) 
 
| Setting | Type | Default | Description | 
|---------|------|---------|-------------| 
| `keep_points` | String | `"canopy"` | Sets lidar point classifications to load. Requires the input point cloud is classified. | 
 
> | Option | Description | 
> |--------|-------------|   
> | `"canopy"` | Recommended for most users; minimizes memory usage and focuses on canopy structure | 
> | `"ground+canopy"` | Not recommended; including ground points can significantly increase memory usage and may not improve transmissivity estimates if a DTM is used for terrain | 
> | `"all"` | Not recommended; includes unclassified points which may be noise and further increases memory usage | 
 
 
### Canopy Structure Enhancement 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|-------------| 
| `trunks` | Boolean | `false` | L2R, C2R | Add opaque tree trunks to synthetic images. Requires preprocessing with `canrad_precalc.jl` to generate tree information files. | 
| `branches` | Boolean | `false` | L2R | Add synthetic branch points to increase canopy density. Adds points between existing points and tree crown centers. Requires preprocessing. | 
| `branch_spacing` | Float | `0.1` | L2R | Spacing interval (meters) for creating synthetic branch points | 
| `cbh` | Float | `2.0` | C2R | Height of canopy base above ground (meters). Used if not providing a `cbhf` (canopy base height file). | 
 
**When to enable trunks:** Important for winter/leaf-off conditions and/or sparse canopies. Setting has a minimal effect in dense evergreen forests (can reduce `trunks_peri` to save time) 
 
**When to enable branches:** Sparse lidar point clouds (<20 pts/m²) to increase tree-crown density 
 
<br> 
 
## Calculation Settings 
 
Control what outputs are calculated and saved. 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|-------------| 
| `calc_trans` | Boolean | `false` | L2R, C2R, T2R | Calculate and save time-varying direct-beam transmissivity.<br>• `false`: Only sky-view fraction is calculated<br>• `true`: Transmissivity time series is calculated | 
| `calc_swr` | Integer | `0` | L2R, C2R, T2R | Shortwave radiation calculation mode:<br>• `0`: No radiation calculation<br>• `1`: Potential radiation (atmospheric transmissivity = 1)<br>• `2`: Actual radiation using measured data (requires `swrf` input file)<br>**Note:** Options 1 and 2 require `calc_trans = true` | 
 
 
<br> 
 
## Time Settings 
 
Define the temporal extent and resolution of calculations. 
 
| Setting | Type | Default | Models | Required | Description | 
|---------|------|---------|--------|----------|-------------| 
| `t1` | String | `"22.06.2020 00:00:00"` | L2R, C2R, T2R | If `calc_trans=true` | Start date and time for calculations | 
| `t2` | String | `"22.06.2020 23:00:00"` | L2R, C2R, T2R | If `calc_trans=true` | End date and time for calculations | 
| `tstep` | Integer | `60` | L2R, C2R, T2R | No | Output time step in minutes (transmissivity and radiation are aggregated to this interval) | 
 
**Format:** `"dd.mm.yyyy HH:MM:SS"` 
 
**Examples:** 
```julia 
# Single day in summer 
t1 = "21.06.2023 00:00:00" 
t2 = "21.06.2023 23:00:00" 
tstep = 30  # 30-minute output 
 
# Full year 
t1 = "01.01.2023 00:00:00" 
t2 = "31.12.2023 23:00:00" 
tstep = 60  # hourly output 
``` 
 
**Important notes:** 
- All calculations are performed at 2-minute intervals internally 
- Results are then averaged to `tstep` for output 
- Minimum `tstep` value is 2 minutes 
- For `calc_swr = 2`, measurement timestamps must match output `tstep` exactly 
 
## Location Settings 
 
Specify geographic location and coordinate system. 
 
| Setting | Type | Models | Required | Description | 
|---------|------|--------|----------|-------------| 
| `time_zone` | Integer | L2R, C2R, T2R | If `calc_trans=true` | Local time zone offset from UTC (e.g., CET = UTC+1 → `time_zone = 1`) | 
| `epsg_code` | Integer | L2R, C2R, T2R | If `calc_trans=true` | EPSG code for coordinate reference system of all spatial inputs and outputs | 
 
**Examples:** 
```julia 
# Switzerland (LV95) 
epsg_code = 2056 
time_zone = 1  # CET 
 
# Western US (UTM Zone 10N) 
epsg_code = 32610 
time_zone = -8  # PST 
``` 
 
**Critical:** All input datasets (rasters, point clouds, calculation points) must use this coordinate system. 
 
 
## Image Settings 
 
Configure synthetic hemispheric image generation. 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|-------------| 
| `image_height` | Float | `0.5` | L2R, C2R, T2R | Height (meters) of virtual camera position above ground | 
| `point_size` | Vector{Float64} | `[30.0, 10.0]` | L2R | Size of points in synthetic hemispheric images (pixels). First value = points closest to camera, second value = points farthest from camera. Image resolution is 10000×10000 pixels internally. | 
 
**Adjusting image_height:** 
- `0.5`: Standard for radiation sensors or snow surface 
- `1.5-2.0`: For above-ground vegetation sensors 
- Lower values emphasize nearby canopy; higher values give wider perspective 
 
**Adjusting point_size (L2R):** 
- Larger values: Denser canopy representation, higher transmissivity sensitivity 
- Smaller values: May miss canopy gaps in sparse point clouds 
- Depends on point cloud density 
 
## Output and Run Settings 
 
Control what files are saved and how the model runs. 
 
| Setting | Type | Default | Models | Description | 
|---------|------|---------|--------|-------------| 
| `batch` | Boolean | `false` | L2R, C2R, T2R | Enable batch mode for processing multiple tiles with organized sub-folder structure | 
| `save_images` | Boolean | `false` | L2R, C2R, T2R | Save synthetic hemispheric images in NetCDF file | 
| `make_pngs` | Boolean | `false` | L2R, C2R, T2R | Export saved images as PNG files to a separate folder. Note that if `phenology` is enabled, the output is a .gif for each location. | 
| `save_horizon` | Boolean | `false` | C2R, T2R | Save calculated topographic horizon line matrix for reuse | 
| `progress` | Boolean | `false` | L2R, C2R, T2R | Report detailed per-point progress to `/ProgressLastPoint/` folder (useful for debugging or timing L2R runs) | 
 
**Considerations:** 
- `save_images = true` and `make_pngs = true`: Creates many files; useful for validation but storage-intensive 
- `save_horizon = true`: Used for topographically downscaling radiation data in external software
 
## Special Implementations (C2R only) 
 
Pre-configured specialised configurations for specific large-scale applications and datasets. Use `"none"` unless working with these specific projects. 
 
| Option | Description | Use Case | 
|--------|-------------|----------| 
| `"none"` | Standard C2R configuration | General applications (requires `lavdf` input) | 
| `"swissrad"` | Swiss nationwide dataset configuration | Switzerland-wide radiation modeling | 
| `"oshd"` | Uses pre-calculated horizon lines for compatibility | Integration with OSHD downscaled radiation | 
| `"oshd-alps"` | Configuration for European Alps domain | Large-scale Alpine applications | 

 
**Usage:** 
```julia 
special_implementation = "none"  # For most users 
``` 
 
 
 
## Settings File Templates 
 
See files in `testset/` for complete working examples of each model.  
 
 
## Important Notes 
 
- All input datasets but be in the same coordinate reference system specified by `epsg_code` 
- Ensure canopy input data extent covers at least `forest_peri` beyond all calculation points and terrain data covers `highres_peri` and `lowres_peri` as applicable 
- For `calc_swr = 2`, ensure `swrf` measurement timestamps exactly match output `tstep` 
 
## Troubleshooting Settings 
 
**Memory issues:** 
- Reduce `forest_peri` or `trunks_peri` 
- For L2R: Use `keep_points = "canopy"` instead of `"all"` 
- Reduce spatial or temporal extent 
 
**Slow computation:** 
- Disable `branches` if not needed 
- Reduce `trunks_peri` in dense forests 
- Increase `tstep` if high temporal resolution isn't required