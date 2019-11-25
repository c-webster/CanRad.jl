function LAS2Rad_Settings(basefolder,segname)

    par_in = Dict(
    "trunks" => true, # true or false
    "branches" => true, # true or false

    "terrain_peri" => 5000,
    "radius" => 500,

    "tershad" => 1,
    "terrain" => 1,
    "tilt" => true, # true or false

    "ch" => 0.5,

    "t1"  => "2016.12.21 00:00:00",
    "t2"  => "2017.06.22 00:00:00",
    "int" => 10,

    "time_zone" => 1,
    "coor_system" => "CH1903",

    "surf_peri" => 100,

    "tolerance" => 1.5,

    "progress" => true
    )

    dat_in = Dict(
    # las
    "dsmf" => basefolder*"SegmentA.las",
    # dtm
    "dtmf" => "C:/Users/v1cwebst/GoogleDrive/Data/L2HEval/Data_Terrain/DTM/DTM_Fluela.mat",
    # dem
    "demf" => "C:/Users/v1cwebst/GoogleDrive/Data/L2HEval/Data_Terrain/DEM/dem_100m.txt",
    # dbh
    "dbhf" => basefolder*"SegmentA_treeinfo.txt",
    # ltc
    "ltcf" => basefolder*"SegmentA_lastreeclass.txt"
    )

    return par_in, dat_in


end
