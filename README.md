# CanRad

[![Build Status](https://travis-ci.com/c-webster/CanRad.jl.svg?branch=master)](https://travis-ci.com/c-webster/CanRad.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/c-webster/CanRad.jl?svg=true)](https://ci.appveyor.com/project/c-webster/CanRad-jl)
[![Codecov](https://codecov.io/gh/c-webster/CanRad.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/c-webster/CanRad.jl)
[![Coveralls](https://coveralls.io/repos/github/c-webster/CanRad.jl/badge.svg?branch=master)](https://coveralls.io/github/c-webster/CanRad.jl?branch=master)


## Description

CanRad (CanopyRadiation) calculates direct and diffuse shortwave radiation transmission in forest canopies using either airborne lidar data (L2R) or a canopy height model (C2R). 

For a detailed description of L2R, see:
Webster C, Mazzotti G, Essery E and Jonas T 2020, [Enhancing airborne lidar data for improved forest structure representation in shortwave transmission models](https://doi.org/10.1016/j.rse.2020.112017), Remote Sensing of Environment 

L2R is an optimised and improved version of [Lidar2HemiEval](https://github.com/c-webster/Lidar2HemiEval), which was based on [Lidar2Hemi](https://github.com/Tobias-Jonas-SLF/Lidar2Hemi), written in MATLAB. 

For a detailed description of C2R, see:
Webster, C, Essery E, Mazzotti G and Jonas T, 2023 [Using just a canopy height model to obtain lidar-level accuracy in 3D forest canopy shortwave transmissivity estimate](https://doi.org/10.1016/j.agrformet.2023.109429), Agricultural and Forest Meteorology

Both models require some additional preperatory steps. Scripts to easily carry out these steps are currently under development to facilitate easy implementation for different sites/users and will be added ~/examples/prep/ as they are completed. In the meantime, please contact the author if you wish to have earlier access. 

Both models can be run on an HPC to cover large domains. An example implementation will be added to ~/examples/cluster/. Feel free to contact the author for more information.


## Installation

CanRad requires [SpatialFileIO.jl](https://github.com/c-webster/SpatialFileIO.jl) and the python package `scipy`. 

`SpatialFileIO` can be added by
```
]add https://github.com/c-webster/SpatialFileIO.jl
```

`scipy` can be added by
```
]add Conda
using Conda
Conda.add("scipy")
```

Then add `CanRad`
```
]add https://github.com/c-webster/CanRad.jl
```

## Tests

Add required packages
```
]add DelimitedFiles, NCDatasets
```

Run the test (this will take several minutes)
```
]test CanRad
```

## Usage

Use ~/test/L2R_Settings_test.jl or ~/test/C2R_Settings_test.jl for desired input parameters and file paths. 
Edit ~/test/run_CanRad_tests.jl 

... examaples for points, batch and cluster to be added


## Models

### D2R - Digital Terrain Model to Radiation




### L2R - Lidar to Radiation




### C2R - Canopy Height Model to Radiation





## Outputs

### Time-static variables

**Vf_planar**: Sky-view fraction calcuated for the perspective of a horizontal flat uplooking surface (sensor or ground). Calculated by weighting zenith rings in synthetic hemispheric image according to their surface area projected onto a flat horizontal surface

**Vf_hemi**: Sky-view fraction calculated for the perspective of a hemispherically shaped surface (sensor or e.g. plant). Calculated by weighting zenith rings in synthetic hemispheric image according to their surface area on the hemisphere

### Time-varying variables

All variables are calculated at 2 minute intervals to ensure a temporally complete solar track, then averaged to the user-defined interval. The output time stamps are the beginning of the averged period. 

**Forest_Transmissivity**: Time-varying direct beam shortwave transmissivity. Calculated by determining ratio of canopy/sky pixels in front of the projected position of the solar disc.

**SWR_total**: Maximum potential total incoming shortwave radiation (diffuse + direct). 

**SWR_direct**: Maximum potential direct incoming shortwave radiation. 

Shortwave radiation variables use atmospheric transmissivity = 1.

Diffuse radiation can be calculated using the difference between total and direct radiation. 


### Synthetic hemispheric images

The settings file has an option to save the calculated synthetic images in a netcdf file. Use 
```
make_SHIs(/path/to/SHIs)
```
to create .pngs for each image. 

If *seasonal* was enabled in L2R, the output is a .gif for each location. 