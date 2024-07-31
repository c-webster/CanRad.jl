**variable name [relevant model]**


**dtmf [ L2R, C2R, T2R ]**\
high-resolution digital terrain model\
used to calculate the near horizon-line within 300m\
.tif format

**demf [ L2R, C2R, T2R ]**\
low-resolution digital terrain model\
used to calculate the distant topographic horizon-line\
.tif format

**chmf [ C2R ]**\
canopy height model\
.tif format

**dsmf [ L2R ]**\
digital surface model\
lidar point cloud
.laz or .las

**lavd [ C2R ]**\
leaf-area volume density\
an output from canrad_precalc.jl\
.tif format

**dbhf [ L2R, C2R ]**\
diameter at breast height file
an output from canrad_precalc.jl
contains dimensions for including trunks in the synthetic images (dbh)
*_treeinfo.txt file

**ltcf [ L2R ]**\
lidar tree classifcation file
an output from pycrown and canrad_precalc.jl\
.laz file specifying which points belong to the tree crowns in treeinfo

**cbh [ C2R ]**\
canopy base height\
an output from canrad_precalc.jl\
used if crown base height varies over domain (e.g. 60% of total tree height)\
.tif format at same resolution as chmf

**swrf [ L2R, C2R, T2R ]**\
shortwave radiation file\
only needed if calc_swr = 2 in the settings file\
format required is tab-delimited files of three columns (date,time,value)\
the model assumes the start and end date are the same as what you have specified in the settings file