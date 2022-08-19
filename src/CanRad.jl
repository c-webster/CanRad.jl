module CanRad

using DelimitedFiles, DataStructures, LasIO, LazIO,
    Statistics, Dates, Interpolations, Images,
    DataFrames, Formatting, Distributed, Distributions,
    SpatialFileIO, NCDatasets, Chain, Pkg, MarketTechnicals

using Conda, PyCall

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
include("LAS2Rad_seasonal.jl")

export
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
    get_pkg_version,
    getimagecentre,
    getPhiTht,
    getsundxs,
    getsurfdat,
    getsurfdat_lavd,
    LAS2Rad,
    LAS2Rad_seasonal,
    loaddbh,
    load_hlm,
    loadltc_laz,
    loadltc_txt,
    make_branches,
    make_SHIs,
    netcdf,
    normalise,
    np,
    organise_outf,
    pcd2pol,
    pcd2pol2cart,
    pol2cart,
    preallo_trunks,
    prepcrtdat,
    prepterdat,
    pyinterp,
    scipyspat,
    split_points,
    # trunk_locs,
    trunkpoints,
    utm2deg,
    write_metadata

end # module
