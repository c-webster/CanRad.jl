module CanRad

using FileIO, LasIO, LazIO, DelimitedFiles, MAT, DataStructures,
    Statistics, VectorizedRoutines, Dates,
    Match, DataFrames, Formatting, Distributed

using Conda, PyCall

Conda.add("scipy")
Conda.add("netCDF4")
Conda.add("numpy")

function __init__()
    @eval global pyinterp  = pyimport("scipy.interpolate")
    @eval global scipyspat = pyimport("scipy.spatial")
    @eval global netcdf    = pyimport("netCDF4")
    @eval global np        = pyimport("numpy")
end


include("File_Functions.jl")
include("Prepatory_Functions.jl")
include("Can2Hemi_Functions.jl")
include("Solar_Functions.jl")
include("CHM2Rad.jl")
include("LAS2Rad.jl")
# include("CanRad_Prep.jl")

export
    # CanRad_Prep,
    extension,
    readlas,
    importdtm,
    read_ascii,
    createfiles,
    create_exmat,
    extract,
    clipdat,
    create_tiles,
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
    prepterdat,
    getimagecentre,
    findpairs,
    fillmat,
    findmincol,
    calc_horizon_lines,
    calcThickness,
    pyinterp,
    scipyspat,
    netcdf,
    np,
    LAS2Rad,
    CHM2Rad



end # module
