
# !!! WORKING VERSION !!! DO NOT USE!!!




# Data Preparation and Setup

Before running CanRad models (L2R or C2R), you might need to prepare specific input datasets. This guide walks through the preprocessing workflow using the provided `canrad_precalc.jl` script.


# Basic setup

## Running model at point scale (basic setup)

To run CanRad at specific points without preprocessing large datasets, follow these steps:





## Using Julia's Distributed Computing 

CanRad can leverage Julia's `Distributed` package to run simulations in parallel across multiple CPU cores. This is particularly useful for large datasets or when running multiple scenarios.



## Running on HPC Clusters with SLURM





## Overview of precalc workflow

The precalc workflow consists of two main parts:

1. **Part 1: Tree Detection** — Uses [pycrown](https://github.com/manaakiwhenua/pycrown) to detect individual trees and delineate crown boundaries
2. **Part 2: Dataset Generation** — Calculates CanRad-specific datasets from tree information

**Output datasets:**
- For **L2R**: Tree information file (`*_treeinfo.txt`), tree-classified point cloud (`*_lastreeclass.laz`)
- For **C2R**: Tree information file (`*_treeinfo.txt`), leaf area volume density (`lavd.tif`), optional canopy base height

## Prerequisites

### Required Software

| Tool | Purpose | Installation |
|------|---------|--------------|
| **Python with pycrown** | Tree crown detection | `conda env create -f environment.yml` from [pycrown repository](https://github.com/manaakiwhenua/pycrown) |
| **Julia** | Data processing | [julialang.org](https://julialang.org/downloads/) |
| **LAStools** (L2R only) | Point cloud classification | [lastools.github.io](https://lastools.github.io/) |

### Required Julia Packages

Install these packages in Julia (press `]` to enter package mode):

```julia
]add Conda
using Conda
Conda.add("numpy")

]add https://github.com/c-webster/SpatialFileIO.jl
]add PolygonOps, StaticArrays, DelimitedFiles, Statistics, SpatialFileIO, PyCall
]add GeoJSON@0.5.1  # Important: Must be version 0.5.1
]add ArchGDAL
```

⚠️ **Important:** GeoJSON must be version 0.5.1 or the script will fail.

### Required Input Data

For tree detection (Part 1), you need four datasets with **identical extent and resolution**:

| File | Description | How to Generate |
|------|-------------|-----------------|
| `CHM.tif` | Canopy Height Model | From lidar using PDAL, lidR, or LAStools |
| `DTM.tif` | Digital Terrain Model | From lidar ground points |
| `DSM.tif` | Digital Surface Model | From lidar first returns |
| `POINTS.laz` | Lidar point cloud | Original data or subset |

**Best practice:** Generate CHM/DTM/DSM directly from your lidar point cloud to ensure perfect alignment. Use PDAL, lidR (R package), or LAStools.

## Directory Structure

Organize your data as follows:

```
your_project/
├── pycrown/
│   ├── data/
│   │   ├── CHM.tif
│   │   ├── DTM.tif
│   │   ├── DSM.tif
│   │   └── POINTS.laz
│   ├── result/          # Created by pycrown
│   └── your_settings.py # Pycrown configuration
└── your_data.laz        # Full point cloud (for L2R)
```

See `setup/example/` in the CanRad repository for a working example.

## Step-by-Step Workflow

### Step 1: Configure Pycrown Settings

Create a Python script (e.g., `your_settings.py`) to configure pycrown for your data. Key parameters to adjust:

```python
from pycrown import PyCrown

PC = PyCrown(F_CHM, F_DTM, F_DSM, F_LAS, outpath='./result')

# Smooth CHM (adjust window size for your resolution)
PC.filter_chm(3, ws_in_pixels=False)  # 3m median filter

# Tree detection (adjust parameters for your forest)
PC.tree_detection(PC.chm, ws=3, ws_in_pixels=True, hmin=8.0)  # hmin = minimum tree height

# Crown delineation (tune thresholds for your site)
PC.crown_delineation(algorithm='watershed_skimage', 
                     th_tree=8.0,    # Tree height threshold
                     th_seed=0.7,    # Seed threshold
                     th_crown=0.55,  # Crown threshold
                     max_crown=10.0) # Max crown diameter (m)

# Screen small trees
PC.screen_small_trees(hmin=5.0, loc='top')
```

**Critical parameters to test:**
- `hmin`: Minimum tree height (m) — depends on your forest
- `ws`: Window size for local maximum filter — smaller for dense forests
- `th_seed`, `th_crown`: Control crown boundary detection — test visually
- `max_crown`: Maximum crown diameter — prevents over-segmentation

**Testing pycrown:** Run pycrown on a small test area and visually inspect the crown polygons. Only proceed to Part 2 once you're satisfied with tree detection quality.

See `setup/example/pycrown/example.py` for a complete configuration example.

### Step 2: Configure canrad_precalc.jl

Edit the user settings section in `setup/canrad_precalc.jl`:

```julia
# Working directory containing your pycrown folder
wrkdir = "/path/to/your_project"

# Path to shp2json.py helper script
shp2json = "../setup/shp2json.py"

# Path to LAStools (for L2R only)
lastoolspath = "/path/to/LAStools/bin"

# Model type: "L2R" or "C2R"
model = "C2R"

# Pycrown settings script name
pycrownscript = "your_settings.py"

# Forest type for DBH estimation (see below)
forest_type = "TCG"

# Tree species for leaf area calculation (C2R only)
tree_species = "Spruce"
```

### Step 3: Set Forest Type for DBH Estimation

The script estimates diameter at breast height (DBH) from tree height and crown diameter using allometric equations from [Jucker et al. (2017)](https://doi.org/10.1111/gcb.13388):

| Code | Forest Type | Region | Example Species | Equation |
|------|-------------|--------|-----------------|----------|
| `"TCG"` | Temperate coniferous gymnosperm | Palearctic | Norway spruce | DBH = 0.974 × (H × CD)^0.748 |
| `"BG"` | Boreal gymnosperm | Palearctic | Scots pine | DBH = 1.43 × (H × CD)^0.649 |
| `"TMG"` | Temperate mixed gymnosperm | Palearctic | Mixed beech/spruce | DBH = 1.001 × (H × CD)^0.73 |

Where H = tree height (m), CD = crown diameter (m), DBH in cm.

### Step 4: Set Tree Species for Leaf Area (C2R Only)

For C2R, the script calculates leaf area (LA) for each tree to generate the leaf area volume density (LAVD) file:

| Species | Equation | Reference |
|---------|----------|-----------|
| `"Spruce"` | LA = exp(-8.31 + 2.61×ln(DBH×10) - 0.07×H) | Goude et al. (2019) |
| `"Scot Pine"` | LA = 0.0091 × DBH^2.363 × H^0.106 | — |
| `"Beech"` | LA = 8.56 + 0.0286 × DBH^2.623 | DOI 10.1007/s10342-009-0345-8 |

**LAVD calculation:**
- Crown volume (CV) is estimated as 80% of a cone (spruce) or 75% (pine)
- LAVD = LA / CV (m²/m³)
- Values are assigned to all CHM pixels within each crown polygon

### Step 5: Run canrad_precalc.jl

Activate the pycrown conda environment and run Julia:

```bash
conda activate pycrown
julia
```

In Julia, run the preprocessing script:

```julia
include("setup/canrad_precalc.jl")
```

**What happens:**

**Part 1 — Tree Detection:**
1. Runs your pycrown Python script
2. Detects tree tops and delineates crown boundaries
3. Exports crown polygons as shapefiles
4. Converts shapefiles to GeoJSON (using `shp2json.py`)

**Part 2 — Dataset Generation:**

**For both L2R and C2R:**
- Reads tree crown GeoJSON files
- Matches tree tops to crown polygons
- Calculates crown area and diameter
- Estimates DBH using forest-specific allometric equations
- Exports tree information file: `*_treeinfo.txt`

**For L2R additionally:**
- Re-classifies point cloud with ground points using LAStools
- Exports: `*_lastreeclass.laz`

**For C2R additionally:**
- Calculates leaf area for each tree (species-specific)
- Estimates crown volume from crown dimensions
- Computes LAVD = leaf area / crown volume
- Assigns LAVD values to CHM pixels within crowns
- Fills gaps: pixels above 2m with no LAVD get the mean LAVD value
- Exports: `lavd.tif` (same extent/resolution as CHM)

### Step 6: Verify Outputs

Check that the following files were created:

**L2R outputs:**
- `*_treeinfo.txt` — Tree information (x, y, height, DBH, tree ID)
- `*_lastreeclass.laz` — Point cloud with ground classification

**C2R outputs:**
- `*_treeinfo.txt` — Tree information
- `pycrown/lavd.tif` — Leaf area volume density raster

**Tree information file format:**
```
treeptsX    treeptsY    treeptsH_m    treeptsdbh_cm    tree_num
2647123.45  1185678.90  15.2          25.3             1
2647128.12  1185682.34  18.6          32.1             2
```

## Troubleshooting

### Pycrown Issues

**No trees detected:**
- Lower `hmin` parameter
- Decrease `ws` (window size) for denser forests
- Smooth CHM less aggressively

**Over-segmentation (too many small trees):**
- Increase `hmin` and `ws`
- Adjust `th_seed` and `th_crown` thresholds
- Increase minimum tree screening height

**Poor crown boundaries:**
- Try different algorithms: `'watershed_skimage'`, `'dalponte_cython'`
- Adjust `th_crown` parameter
- Check CHM quality (gaps, noise)

### canrad_precalc.jl Issues

**"forest type unknown" error:**
- Check `forest_type` is one of: `"TCG"`, `"BG"`, `"TMG"`

**"tree species not defined" (C2R):**
- Check `tree_species` is one of: `"Spruce"`, `"Scot Pine"`, `"Beech"`

**GeoJSON version error:**
- Ensure GeoJSON.jl version 0.5.1: `]add GeoJSON@0.5.1`

**Missing lastreeclass.laz (L2R):**
- Verify `lastoolspath` is correct
- Check LAStools is installed and accessible
- On Linux, ensure `wine` is installed for running .exe files

**CHM/DSM/DTM extent mismatch:**
- Regenerate all three rasters with identical parameters
- Use same resolution, extent, and CRS

### Data Quality Issues

**Sparse LAVD values:**
- Check pycrown crown polygons cover the canopy
- Verify CHM has no large data gaps
- Consider adjusting gap-filling threshold (currently 2m)

**Unrealistic DBH values:**
- Verify correct `forest_type` for your region
- Check tree height and crown diameter values are reasonable
- Consider calibrating with field measurements

## Example Test Data

A complete example is provided in `setup/example/` with the same data used in `testset/input/`. Use this to:
- Learn the directory structure
- Test your installation
- See correctly formatted inputs/outputs
- Understand pycrown parameter effects

## Additional Notes

### For L2R: Why Re-classify Points?

Pycrown creates a new point cloud with tree ID classifications but doesn't preserve original ground/vegetation classifications. The script re-runs the point cloud through `lasground.exe` to restore ground classification. This is critical because L2R discards ground points when loading data, significantly reducing memory requirements.

### For C2R: Crown Height Assumptions

The script assumes:
- **Spruce:** Crown base at 20% of tree height (crown = 80% of tree)
- **Scots Pine:** Crown base at 25% of tree height (crown = 75% of tree)

If your forest has different crown proportions, edit these values in the script (lines dealing with crown volume calculation).

### Custom Allometric Equations

To add your own DBH or leaf area equations:

1. Add a new `forest_type` or `tree_species` option
2. Insert the equation in the appropriate section
3. Cite the source in comments

**Example:**
```julia
elseif forest_type == "MY_FOREST"
    dbh_t = [your equation using cz_t and cd_m]
```

### Coordinate Systems

All input rasters and point clouds must share the same projected coordinate reference system (CRS). The script reads CRS from the CHM file and applies it to outputs.

## Next Steps

After successfully generating preprocessed datasets:

1. Review the [Input Data](input.md) page to understand what each file contains
2. Create a settings file for your model — see [Configuration Settings](settings.md)
3. Specify file paths in your settings file
4. Run CanRad — see [Models](models.md) for execution instructions

## References

- **Pycrown:** Zörner, J., Dymond, J. R., Shepherd, J. D., & Jolly, B. (2018). PyCrown - Fast raster-based individual tree segmentation for LiDAR data. DOI: [10.7931/M0SR-DN55](https://doi.org/10.7931/M0SR-DN55)
- **Allometric equations:** Jucker, T., et al. (2017). Allometric equations for integrating remote sensing imagery into forest monitoring programmes. *Global Change Biology*, 23(1), 177-190. DOI: [10.1111/gcb.13388](https://doi.org/10.1111/gcb.13388) 