function LAS2Rad_Settings(basefolder,segname="empty")

    par_in = Dict(
    "trunks" => true,
    "branches" => true,

    "terrain_peri" => 7000,
    "radius" => 500,

    "tershad" => 1,
    "terrain" => 1,
    "tilt" => false,

    "ch" => 0.5,

    "t1"  => "2016.12.21 00:00:00",
    "t2"  => "2017.06.22 00:00:00",
    "int" => 10,

    "time_zone" => 1,
    "coor_system" => "CH1903",

    "surf_peri" => 100,

    "nrbins" => 5,

    "tolerance" => 1.5,

    "progress" => false
    )

    dat_in = Dict(
    # las
    "dsmf" => basefolder*"/Segments/"*segname*"/"*segname*".las",
    # dtm
    "dtmf" => basefolder*"/Data_Terrain/DTM/DTM_LaretLarge.mat",
    # dem
    "demf" => basefolder*"/Data_Terrain/DEM/dem_100m.txt",
    # dbh
    "dbhf" => basefolder*"/Segments/"*segname*"/"*segname*"_treeinfo.txt",
    # ltc
    "ltcf" => basefolder*"/Segments/"*segname*"/"*segname*"_lastreeclass.txt"
    )

    return par_in, dat_in


end
