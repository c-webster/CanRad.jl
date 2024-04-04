
You can use canrad_precalc.jl to calculate the specific datasets required by L2R or C2R. 


The script starts by running the python package [pycrown](https://doi.org/10.7931/M0SR-DN55) (Part 1), and requires four datasets requires four datasets: CHM.tif, DTM.tif, DSM.tif and POINTS.laz\
You can find more details on pycrown in the [github repository](https://github.com/manaakiwhenua/pycrown) but the main point is that the CHM/DTM/DSM must have the same extent and resolution. It's easiest to calculate these files directly from your lidar point cloud using e.g. pdal, lidR or LAStools\
You also need the .py script to run pycrown. \
The settings in example.py are specific to the example dataset. It is recommended that you test the pycrown settings specific to your area. Only once you get the desired tree crown segmentation should you move on to calculating the datasets for CanRad (Part 2)


canrad_precalc.jl is written in julia, but you should run it within a python environment where pycrown is installed:

Example in a terminal: 
```
conda env create -f environment.yml 
conda activate pycrown
julia
```
where environment.yml refers to the file from the [pycrown respository](https://github.com/manaakiwhenua/pycrown)

Then install the relevant packages in julia:
```
]add Conda
using Conda
Conda.add("numpy")

]add https://github.com/c-webster/SpatialFileIO.jl
]add PolygonOps, StaticArrays, DelimitedFiles, Statistics, SpatialFileIO, PyCall
]add GeoJSON@0.5.1
```
where ```]``` activates the julia package manager


The example dataset given in the folder ..setup/example is the same dataset as that in ../testset/input
Using the scripts and datasets in the ```setup``` and ```testset``` folders should give you enough information to run CanRad over your domain. 


Note that running canrad_precalc.jl for L2R requires [LAStools](https://lastools.github.io/). This is because the output point cloud from pycrown does not preserve the point classifications, so the script re-runs the point cloud through lasground.exe. When L2R loads the point cloud, ground points are discarded, drastically reducing RAM requirements. 