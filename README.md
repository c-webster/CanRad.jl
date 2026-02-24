# The Canopy Radiation Model

[![Build Status](https://travis-ci.com/c-webster/CanRad.jl.svg?branch=master)](https://travis-ci.com/c-webster/CanRad.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/c-webster/CanRad.jl?svg=true)](https://ci.appveyor.com/project/c-webster/CanRad-jl)
[![Codecov](https://codecov.io/gh/c-webster/CanRad.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/c-webster/CanRad.jl)
[![Coveralls](https://coveralls.io/repos/github/c-webster/CanRad.jl/badge.svg?branch=master)](https://coveralls.io/github/c-webster/CanRad.jl?branch=master)


## Description

The Canopy Radiation Model (CanRad) calculates direct and diffuse shortwave radiation transmission in forest canopies using either airborne lidar data (L2R) or a canopy height model (C2R). An additional module is available that calculates radiation transmission using only terrain data (T2R).

For a detailed description of L2R, see:
Webster C, Mazzotti G, Essery E and Jonas T 2020, [Enhancing airborne lidar data for improved forest structure representation in shortwave transmission models](https://doi.org/10.1016/j.rse.2020.112017), Remote Sensing of Environment 

L2R is an optimised and improved version of [Lidar2HemiEval](https://github.com/c-webster/Lidar2HemiEval), which was based on [Lidar2Hemi](https://github.com/Tobias-Jonas-SLF/Lidar2Hemi), written in MATLAB. 

For a detailed description of C2R, see:
Webster, C, Essery E, Mazzotti G and Jonas T, 2023 [Using just a canopy height model to obtain lidar-level accuracy in 3D forest canopy shortwave transmissivity estimates](https://doi.org/10.1016/j.agrformet.2023.109429), Agricultural and Forest Meteorology

Both models require some additional preperatory steps. Scripts to easily carry out these steps are currently under development to facilitate easy implementation for different sites/users and will be added ~/examples/prep/ as they are completed. In the meantime, please contact the author if you wish to have earlier access. 

Both models can be run on an HPC to cover large domains. An example implementation will be added to ~/examples/cluster/. Feel free to contact the author for more information.


## Installation

CanRad requires [SpatialFileIO.jl](https://github.com/c-webster/SpatialFileIO.jl) and the python package `scipy`. 

`SpatialFileIO` can be added by
```
using Pkg
Pkg.add(url="https://github.com/c-webster/SpatialFileIO.jl")
```

`scipy` can be added by
```
Pkg.add("Conda")
using Conda
Conda.add("scipy")
```

Then add `CanRad`
```
Pkg.add(url="https://github.com/c-webster/CanRad.jl")
```

## Tests

Tests compilation and basic usage of `CanRad`. The test can take several minutes to run as they include a full L2R run on a small test site, although testing LAS2Rad is skipped if there is not enough RAM available for loading in the test point cloud.

The output sky-view factor from the test runs of L2R and C2R are benchmarked against values calculated from hemispherical photographs taken at the same locations using [HPEVal](https://github.com/Tobias-Jonas-SLF/HPEval) and available [here](https://github.com/c-webster/CanRad.jl/tree/main/testset/real_HPs). The tests are considered passed if the code runs without error and the sky-view factor values from L2R and C2R are within 0.05 of the values from HPEVal. 

Add required packages
```
Pkg.add(["DelimitedFiles", "NCDatasets"])
```

Run the test
```
using CanRad
Pkg.test("CanRad")
```

<!-- ## Usage

See specific documnentation for L2R and C2R in the [docs](https://c-webster.github.io/CanRad.jl/stable/). -->



<!-- If *phenology* was enabled in L2R, the output is a .gif for each location.  -->





