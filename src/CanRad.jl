module CanRad

using FileIO, LasIO, DelimitedFiles, MAT, DataStructures,
    Statistics, VectorizedRoutines, Dates,
    Match, DataFrames, Formatting

using PyCall

function __init__()
    @eval global pyinterp  = pyimport("scipy.interpolate")
    @eval global scipyspat = pyimport("scipy.spatial")
    @eval global netcdf    = pyimport("netCDF4")
    @eval global np        = pyimport("numpy")
end

include("File_Functions.jl")
include("calcHorizon.jl")
include("Can2Hemi_Functions.jl")
include("Prepatory_Functions.jl")
include("Can2Hemi_Functions.jl")
include("Solar_Functions.jl")
include("Terrain_Functions.jl")

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
    makeDTM,
    tileDTM,
    Settings_makeDTM

include("Settings_Files/Settings_makeDTM.jl")


end # module
