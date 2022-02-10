module CanRad

using DelimitedFiles, DataStructures, LasIO, LazIO,
    Statistics, Dates, Interpolations,
    DataFrames, Formatting, Distributed, Distributions,
    SpatialFileIO, NCDatasets

using Conda, PyCall

# Conda.add("scipy")

function __init__()
    @eval global pyinterp  = pyimport("scipy.interpolate")
    @eval global scipyspat = pyimport("scipy.spatial")
end


include("File_Functions.jl")
include("Preparatory_Functions.jl")
include("Can2Hemi_Functions.jl")
include("Solar_Functions.jl")
include("CHM2Rad.jl")
include("LAS2Rad.jl")

export
    # CanRad_Prep,
    aggregate_data,
    calc_horizon_lines,
    calc_solar_track,
    calc_transmissivity,
    calcCHM_Ptrans,
    calcmintht,
    calcmm,
    calcPtrans,
    calcPtrans_dist,
    calcringratios,
    calcThickness,
    calculate_trunks,
    calculateSWR,
    calculateVf,
    cart2sph,
    check_output,
    CHM2Rad,
    clipdat,
    compatability_check,
    correct_sph,
    create_exmat,
    create_mat,
    create_tiles,
    createfiles,
    createVariables,
    dem2pol,
    dist,
    dist3d,
    extension,
    extract,
    fillmat,
    fillterrain,
    filterbyradius,
    findelev,
    findmincol,
    findpairs,
    frbins,
    get_constants,
    getimagecentre,
    getPhiTht,
    getsundxs,
    getsurfdat,
    getsurfdat_lavd,
    LAS2Rad,
    loaddbh,
    load_hlm,
    loadltc_laz,
    loadltc_txt,
    make_branches,
    netcdf,
    normalise,
    np,
    pcd2pol,
    pcd2pol2cart,
    pol2cart,
    preallo_trunks,
    prepcrtdat,
    prepterdat,
    pyinterp,
    scipyspat,
    # trunk_locs,
    trunkpoints,
    utm2deg

end # module
