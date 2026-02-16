# Model Outputs

CanRad produces several types of output variables that characterize radiation transmission through forest canopies and terrain. All outputs are saved in NetCDF format.

## Output Variables

### Sky-View Factor (SVF)

Sky-view factor quantifies the proportion of the sky hemisphere that is visible from a given location, accounting for obstructions from canopy and terrain.

| Variable | Description | Calculation Method |
|----------|-------------|-------------------|
| `svf_planar` | Planar sky-view factor | Calculated for a horizontal flat surface (e.g., radiation sensor, snow surface). Zenith rings in the synthetic hemispheric image are weighted by their surface area projected onto a flat horizontal plane. |
| `svf_hemi` | Hemispherical sky-view factor | Calculated for a hemispherically shaped surface (e.g., plant canopy). Zenith rings are weighted by their surface area on the hemisphere itself. |

**Use cases:**
- `svf_planar`: Appropriate for modeling radiation on flat surfaces (snow, soil, horizontal sensors)
- `svf_hemi`: Appropriate for modeling radiation on curved surfaces (plant leaves, spherical sensors)

### Time-Varying Transmissivity (TVT)

Time-varying direct beam shortwave transmissivity represents the fraction of direct solar radiation that passes through the canopy at each moment in time.

**Calculation:**
- The model determines the ratio of sky pixels to total pixels in front of the projected position of the solar disc
- Computed at **2-minute intervals** to ensure complete temporal coverage of the solar track
- Averaged to the user-defined time step (`tstep` in settings)
- Output timestamps represent the **beginning** of each averaged period

**Variable naming:** `tvt_<environment>` (e.g., `tvt_evergreen`, `tvt_leafoff`)

### Shortwave Radiation (SWR)

CanRad can output shortwave radiation fluxes in two ways: calculated internally or post-processed from transmissivity data.

#### Output Variables

| Variable | Description |
|----------|-------------|
| `SWR_total` | Total incoming shortwave radiation (direct + diffuse) |
| `SWR_direct` | Direct beam shortwave radiation along the path of the sun |
| `SWR_diffuse` | Diffuse shortwave radiation from the sky hemisphere |

#### Calculation Methods

**Option 1: Calculated within CanRad (T2R/L2R/C2R)**

When `calc_swr = 1` in settings:
- Calculated at 2-minute intervals (matching transmissivity), then averaged to `tstep`
- Computes **potential radiation** using atmospheric transmissivity = 1 (clear-sky maximum)
- `SWR_diffuse` is not saved separately (calculate as: `SWR_total - SWR_direct`)

When `calc_swr = 2` in settings:
- Uses input measured above-canopy radiation to compute compute **actual radiation** 
- Calculated at the temporal resolution of the input data
- ⚠️ **Important:** Above-canopy data time steps must exactly match the transmissivity data in the NetCDF file

**Option 2: Post-processed with CalcSWRFromNC.jl**\
*TO DO* 

### Synthetic Hemispheric Images

CanRad generates synthetic hemispheric (fisheye) images that represent the view looking upward from each calculation point.

**Saving images:**
- Enable in settings file with `save_shi = true`
- Images stored in NetCDF file and can be exported as PNG files (`make_pngs = true`)
- Useful for visualization and validation against field photographs

**Image characteristics:**
- Circular fisheye projection (180° field of view)
- Black pixels = obstructions (canopy, terrain)
- White pixels = visible sky
- Image resolution is 1000 x 1000 pixels

## Environment Suffixes

All output variables include a suffix indicating the environment type, characterised by the user input in the settings file. 

| Suffix | Environment | Description |
|--------|-------------|-------------|
| `_terrain` | Terrain only | Above-canopy conditions with no overlying vegetation. Used for open areas or as a baseline for comparison. |
| `_evergreen` | Evergreen | Coniferous or other evergreen forests with no seasonal phenology. Constant canopy density year-round. |
| `_leafoff` | Leaf-off | Deciduous forests during winter/dormant season. Reduced canopy density, visible branch structure. |
| `_leafon` | Leaf-on | Deciduous forests during growing season. Full canopy density with leaves present. |

**Example variable names:**
- `svf_planar_evergreen`
- `tvt_leafoff`
- `SWR_total_terrain`

## Output File Structure

### NetCDF File Organization

**Primary output file:** `Output_<sitename>.nc`

**Dimensions:**
- `x`, `y`: Spatial coordinates
- `time`: Temporal dimension for time-varying outputs

**Variables by type:**
- Static: `svf_planar_*`, `svf_hemi_*` (no time dimension)
- Time-varying: `tvt_*`, `SWR_total_*`, `SWR_direct_*` (x, y, time)

### Additional Output Files

| File | Content | When Created |
|------|---------|--------------|
| `HLM_<sitename>.nc` | Horizon line matrix (zenith angles) | When `save_horizon = true` |
| `SHIs_<sitename>.nc` | Synthetic hemispheric images | When `save_shi = true` |
| `Settings_<sitename>.txt` | Copy of settings file used for run | Always |
| `SHIs_<sitename>/` | Exported PNG images | When `make_pngs = true` |

### Progress Tracking

Every model run creates a progress file e.g. `Processing 10% ... 15 of 150.txt` file in the output directory, which is updated in real-time to reflect the current progress of the model run as it iterates over each point in the task. This file is used to monitor long-running jobs, especially on HPC clusters. Because it is updated as each point is calculated, it can also be used as an indication of how fast/slow the model is running.

A `Progress_<sitename>/` folder can be created (`step_progress = true`) during model runs containing:
- Step-by-step progress updates
- Timing information for performance monitoring/benchmarking
- It is *not recommended* to enable this option for large runs due to the large number of files created. It is designed for development and benchmarking rather than routine use.

## Units and Conventions

| Variable | Units | Range |
|----------|-------|-------|
| Sky-view fraction | Dimensionless | 0–1 (0 = fully obstructed, 1 = open sky) |
| Transmissivity | Dimensionless | 0–1 (0 = fully blocked, 1 = fully transmitted) |
| Shortwave radiation | W m⁻² | 0–1400 (approximately) |
| Time | Hours or custom | Defined by `tstep` setting |

