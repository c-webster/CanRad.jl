# CanRad

[![Build Status](https://travis-ci.com/c-webster/CanRad.jl.svg?branch=master)](https://travis-ci.com/c-webster/CanRad.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/c-webster/CanRad.jl?svg=true)](https://ci.appveyor.com/project/c-webster/CanRad-jl)
[![Codecov](https://codecov.io/gh/c-webster/CanRad.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/c-webster/CanRad.jl)
[![Coveralls](https://coveralls.io/repos/github/c-webster/CanRad.jl/badge.svg?branch=master)](https://coveralls.io/github/c-webster/CanRad.jl?branch=master)


## Usage

Package to calculate shortwave radiation transmission in forest canopies using either airborne lidar data (LAS2Rad) or a canopy height model (CHM2Rad).



## Installation

CanRad requires [SpatialFileIO.jl](https://github.com/c-webster/SpatialFileIO.jl) and the python package `scipy`. CanRad and SpatialFileIO are still private repositories, so must be added using the path to either the private repositories, or the local storage of the packages.

`scipy` can be added by
```
using Pkg
Pkg.add("Conda")
using Conda
Conda.add("scipy")
```
