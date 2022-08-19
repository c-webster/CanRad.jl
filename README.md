# CanRad

[![Build Status](https://travis-ci.com/c-webster/CanRad.jl.svg?branch=master)](https://travis-ci.com/c-webster/CanRad.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/c-webster/CanRad.jl?svg=true)](https://ci.appveyor.com/project/c-webster/CanRad-jl)
[![Codecov](https://codecov.io/gh/c-webster/CanRad.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/c-webster/CanRad.jl)
[![Coveralls](https://coveralls.io/repos/github/c-webster/CanRad.jl/badge.svg?branch=master)](https://coveralls.io/github/c-webster/CanRad.jl?branch=master)


## Description

CanRad (CanopyRadiation) calculates direct and diffuse shortwave radiation transmission in forest canopies using either airborne lidar data (LAS2Rad/L2R) or a canopy height model (CHM2Rad/C2R). 

For a detailed description of L2R, see:
Webster C, Mazzotti G, Essery E and Jonas T 2020, [Enhancing airborne lidar data for improved forest structure representation in shortwave transmission models](https://doi.org/10.1016/j.rse.2020.112017), Remote Sensing of Environment 

L2R is an optimised and improved version of [Lidar2HemiEval](https://github.com/c-webster/Lidar2HemiEval), which was based on [Lidar2Hemi](https://github.com/Tobias-Jonas-SLF/Lidar2Hemi).

A detailed description of C2R is currently under review in Agricultural and Forest Meteorology:
Webster, C, Essery E, Mazzotti G and Jonas T, Using just a canopy height model to obtain lidar-level accuracy in 3D forest canopy shortwave transmissivity estimates. Submitted to AfM June 2022

Both models require some additional preperatory steps. Scripts to easily carry out these steps are currently under development to facilitate easy implementation for different sites/users and will be added ~/setup/prep/ as they are completed. In the meantime, please contact the author if you wish to have earlier access. 

Both models can be run on an HPC to cover large domains. An example implementation will be added to ~/setup/cluster/. Feel free to contact the author for more information.


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

Then add CanRad
```
]add https://github.com/c-webster/CanRad.jl
```


## Usage

Edit ~/test/L2R_Settings_test.jl or ~/test/C2R_Settings_test.jl for desired input parameters and file paths. 
Edit ~/test/run_CanRad_tests.jl 

To run the tests, you also need DelimitedFiles
```
]add DelimitedFiles
```
then
```
include("~/run_CanRad_tests.jl")
```
