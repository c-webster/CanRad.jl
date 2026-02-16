# CanRad.jl

The **Canopy Radiation Model** (CanRad) is a Julia package for calculating direct and diffuse shortwave radiation transmission in forest canopies and complex terrain.

## Overview

CanRad provides three complementary models to simulate radiation transmission in forested landscapes:

- **L2R (Lidar to Radiation)**: Uses airborne lidar point clouds for high-fidelity 3D forest structure representation
- **C2R (Canopy Height Model to Radiation)**: Achieves lidar-level accuracy using the only canopy height model
- **T2R (Terrain to Radiation)**: Calculates radiation transmission for terrain-only (above-canopy) conditions

These models generate synthetic hemispheric images to calculate sky-view factors, time-varying transmissivity, and incoming shortwave radiation at any point in the landscape.

## Key Features

- **High spatial resolution**: Calculate radiation at up to 1 m resolution over large domains
- **Time-varying transmissivity**: Calculate direct beam transmission throughout daily and annual solar cycles
- **Terrain integration**: Include local (1-25 m) and regional (>25 km) topographic shading
- **Phenology support**: Model seasonal changes in deciduous and mixed forests 
- **Scalable computing**: Run large domains using Julia's distributed package or use slurm on high-performance computing clusters
- **Flexible outputs**: Generate sky-view factors, transmissivity time series, and synthetic hemispheric images with customizable time steps

## When to Use Each Model

| Model | Best For | Key Requirement |
|-------|----------|-----------------|
| **L2R** | Highest accuracy; research applications | Classified lidar point clouds |
| **C2R** | Large-scale applications; operational use | Canopy height model + forest type information |
| **T2R** | Above-canopy conditions; model benchmarking | Terrain data only |

## Quick Installation

CanRad requires [SpatialFileIO.jl](https://github.com/c-webster/SpatialFileIO.jl) and the Python package `scipy`.

```julia
using Pkg

# Install SpatialFileIO
Pkg.add(url="https://github.com/c-webster/SpatialFileIO.jl")

# Install scipy via Conda
Pkg.add("Conda")
using Conda
Conda.add("scipy")

# Install CanRad
Pkg.add(url="https://github.com/c-webster/CanRad.jl")
```

## Basic Usage (applies to both L2R and C2R)

```julia
using CanRad, DelimitedFiles

# Load settings file (see manual/settings for details)
include("C2R_settings.jl")

# Load points at which to run the model
pts = readdlm("pts.txt")

# Generate the settings dictionaries
dat_in, par_in = chm2rad_settings()

# Specify output directory
exdir = "C2R_output"

# (optional) give the model run a taskID
taskID = "C2R_test"

# Run the C2R model
chm2rad!(pts,dat_in,par_in,exdir,taskID)
```

See the [Setup Guide](manual/setup.md) for detailed preprocessing steps and the [Settings Reference](manual/settings.md) for configuration options.

## Documentation Contents

```@contents
Pages = [
    "manual/setup.md",
    "manual/models.md",
    "manual/input.md",
    "manual/settings.md",
    "manual/output.md"
]
Depth = 2
```

## Scientific References

### L2R Model
Webster, C., Mazzotti, G., Essery, R., & Jonas, T. (2020). [Enhancing airborne lidar data for improved forest structure representation in shortwave transmission models](https://doi.org/10.1016/j.rse.2020.112017). *Remote Sensing of Environment*, 249, 112017.

### C2R Model
Webster, C., Essery, R., Mazzotti, G., & Jonas, T. (2023). [Using just a canopy height model to obtain lidar-level accuracy in 3D forest canopy shortwave transmissivity estimates](https://doi.org/10.1016/j.agrformet.2023.109429). *Agricultural and Forest Meteorology*, 333, 109429.

### Larger-scale implementations (C2R + T2R)
Webster, C., Ginzler, C., Marty, M., Nussbaumer, A., Mazzotti, G. and Jonas, T. (2025). [Hourly potential light availability maps at 10 m resolution over Switzerland](https://doi.org/10.1038/s41597-025-06152-9). *Scientific Data*, 12(1), p.1882.

### Model Usage



## Related Software

CanRad builds upon and improves:
- [Lidar2HemiEval](https://github.com/c-webster/Lidar2HemiEval) - An earlier version of L2R implemented in MATLAB that enhances the lidar point cloud with trunk and branch points, and evalauates the resulting hemispheric images for radiation modeling
- [Lidar2Hemi](https://github.com/Tobias-Jonas-SLF/Lidar2Hemi) - The original MATLAB implementation generating synthetic hemispheric images from lidar data

## Getting Help

- Check the [Settings Reference](manual/settings.md) for configuration options
- Review test cases in `testset/` for working examples
- Contact the package maintainer for questions about implementation

## System Requirements

- Julia 1.12 or later
- For L2R setup: [LAStools](https://lastools.github.io/) for point cloud processing (see manual/setup.md for details)
- For C2R setup: Python with [pycrown](https://github.com/manaakiwhenua/pycrown) (see manual/setup.md for details)
