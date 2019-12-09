
# wrkdir     = "C:/Users/v1cwebst/GoogleDrive/Data/LAS2Rad/Scripts/"

setf       = "C:/Users/webster/GDrive/Data/LAS2Rad/Settings/LAS2Rad_Settings_Pts.jl"

ptsf       = "E:/Data/L2HEval/Data_Batches/FluelaSlopes/SegmentA/SegmentA_gridpts_Copy.txt"

basefolder = "E:/Data/L2HEval/Data_Batches/FluelaSlopes/SegmentA/"

exdir      = "C:/Users/webster/GDrive/Data/LAS2Rad/temp/"

###############################################################################
### BEGIN

using CanRad, DelimitedFiles

include(setf)

par_in, dat_in = LAS2Rad_Settings(basefolder)

pts = readdlm(ptsf)

LAS2Rad(pts,dat_in,par_in,exdir)
