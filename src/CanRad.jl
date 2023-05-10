module CanRad

using DelimitedFiles, DataStructures, LasIO, LazIO,
    Statistics, Dates, Interpolations, Images,
    DataFrames, Formatting, Distributed, Distributions,
    SpatialFileIO, NCDatasets, Chain, Pkg, MarketTechnicals,
    Parameters

using Conda, PyCall

function __init__()
    @eval global pyinterp  = pyimport("scipy.interpolate")
    @eval global scipyspat = pyimport("scipy.spatial")
end

include("constants.jl")
include("fileio.jl")
include("initialise.jl")
include("hemifun.jl")
include("solar.jl")
include("chm2rad.jl")
include("las2rad.jl")
include("ter2rad.jl")

export
    pyinterp,
    scipyspat,
    CANRAD,
    SOLAR,
    RADIATION,
    TER2RAD,
    ter2rad!,
    CHM2RAD,
    chm2rad!,
    LAS2RAD,
    las2rad!,
    aggregate_data,
    calc_horizon_lines,
    calc_solar_track,
    calc_transmissivity!,
    calcCHM_Ptrans!,
    calc_latlon,
    calcmintht!,
    calcPtrans,
    calcPtrans_dist,
    calcringratios,
    calcThickness!,
    calculate_trunks,
    calculateSWR,
    calculateVf,
    cart2sph!,
    check_output,
    clipdat,
    clipdat!,
    compatability_check!,
    create_exmat,
    create_tiles,
    createfiles,
    create_exhlm,
    dist,
    dist3d,
    extract,
    fillmat!,
    fillterrain,
    filterbyradius,
    findelev,
    findelev!,
    findmincol,
    findpairs,
    frbins,
    get_pkg_version,
    getimagecentre,
    getlimits!,
    getPhiTht,
    getsundxs!,
    getsurfdat,
    getsurfdat_chm,
    getterrainmask,
    hlm2cart,
    loaddbh,
    load_hlm,
    load_hlm_oshd,
    loadltc_laz,
    loadltc_txt,
    make_branches,
    make_SHIs,
    normalise!,
    normalise_chmdat!,
    organise_outf,
    pcd2pol2cart!,
    pol2cart,
    pol2cart!,
    preallo_trunks,
    prepsurfdat!,
    prepterdat!,
    prepterdat,
    trunkpoints,
    utm2deg,
    write_metadata

end # module
