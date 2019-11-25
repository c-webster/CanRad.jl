function CHM2Rad_Settings(basefolder,segname)


    par_in = Dict(

    "terrain_peri" => 5000,
    "dtm_peri" => 400, # depends on the local terrain...
    "chm_peri" => 100,
    "radius" => 500,

    "tht_int" => 0.5, # 1, 0.5 or 0.25

    "tershad" => 1,
    "terrain" => 1,

    "ch" => 0.5,

    "t1"  => "2016.12.21 00:00:00",
    "t2"  => "2017.06.22 00:00:00",
    "int" => 10,

    "time_zone" => 1,
    "coor_system" => "CH1903",

    "surf_peri" => 100,

    "progress" => true
    )




    dat_in = Dict(
    # chm
    "chmf" => basefolder*
    # dtm
    "dtmf" => "C:/Users/v1cwebst/GoogleDrive/Data/L2HEval/Data_Terrain/DTM/DTM_Fluela.mat",
    # dem
    "demf" => "C:/Users/v1cwebst/GoogleDrive/Data/L2HEval/Data_Terrain/DEM/dem_100m.txt",
    )









end
