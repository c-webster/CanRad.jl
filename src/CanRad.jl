module CanRad

using FileIO, LasIO, DelimitedFiles, MAT, DataStructures,
    Statistics, VectorizedRoutines, Dates,
    Match, DataFrames, Formatting

using Conda, PyCall

function __init__()
    @eval global pyinterp  = pyimport("scipy.interpolate")
    @eval global scipyspat = pyimport("scipy.spatial")
    @eval global netcdf    = pyimport("netCDF4")
    @eval global np        = pyimport("numpy")
end

# using RCall

include("File_Functions.jl")
include("Prepatory_Functions.jl")
include("Can2Hemi_Functions.jl")
include("Solar_Functions.jl")
include("Terrain_Functions.jl")
# include("calcHorizon.jl")

export
    readlas,
    importdtm,
    read_ascii,
    createfiles,
    extract,
    clipdat,
    findelev,
    trunkpoints,
    preallo_trunks,
    calculate_trunks,
    make_branches,
    create_mat,
    dem2pol,
    fillterrain,
    pcd2pol,
    pol2cart,
    prepcrtdat,
    preptercart,
    getimagecentre,
    findpairs,
    fillmat,
    findmincol,
    calc_horizon_lines,
    pyinterp,
    scipyspat,
    netcdf,
    np,
    DTM_create,
    DTM_segment,
    DTM_resample

include("Settings_Files/Settings_DTM_create.jl")
include("Settings_Files/Settings_DTM_segment.jl")


end # module
