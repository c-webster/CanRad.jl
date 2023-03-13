function chm2rad_settings(basefolder="empty"::String)

    par_in = Dict(
    "trunks" => false,

    # terrain settings
    "terrain_highres" => true, # include high high resolution local terrain (1-5 m)
    "terrain_lowres" => false, # include low resolution regional terrain (> 25 m)
    "terrain_tile" => false, # calculate terrain_lowres once for the centre point of the whole area
    "terrain_precalc" => true, # if the terrain mask has been pre-calculated by CanRad

    "terrain_peri" => 10000, # radius around points to include lowres terrain

    "season" => "nothing", # see wiki for more information about this feature
    "for_type" => 2, # 1=deciduous broadleaf; 2=needleleaf

    "image_height" => 0.25, #FLOAT NUMBER  !!!!!!!!

    "calc_trans" => true, # calculate time-varying transmissivity
    "calc_swr" => 1, # 0 = off; 1 = potential swr (atm_trans = 1); 2 = real swr (needs swrf in dat_in)

    "t1"  => "01.10.2019 00:00:00", # start time stamp; format "dd.mm.yyyy HH:MM:SS"
    "t2"  => "21.06.2020 23:50:00", # end time stamp
    "tstep" => 2, # in minutes

    "time_zone" => 1, # relative to UTC (e.g CET = UTC+1)
    "coor_system" => "CH1903",
    # possible options "CH1903", "CH1903+" or e.g. "UTM 32 N" for UTM systems (spaces required)

    "save_images" => true, # save the calculated images in a netcdf file

    "surf_peri" => 100, # radius around points to include in image

    "pt_corr" => true, # enable pt correction (if false, canopy is 100% opaque)

    "batch" => false,

    "buildings" => false, # include buildings as opague features (requires a building height model, bhmf, like the chm)

    "progress" => true  # report individual step progress in /ProgressLastPoint
    )

    # all spatial input data must be in the same coordinate system
    dat_in = Dict(
    # chm
    "chmf" => basefolder*"/input/test_chm.tif",
    # dtm -> high res terrain
    "dtmf" => basefolder*"/input/test_dtm.tif",
    # dem -> coarse res terrain
    "demf" => basefolder*"/input/dem_50m.tif",
    # terrain mask (if terrain_precalc = true)
    "terf" => basefolder*"/Output_tests_T2R/SHIs_testset.nc",
    # lavdf
    "lavdf" => basefolder*"/input/test_lavd.asc",
    # dbh
    "dbhf" => basefolder*"/input/test_chm_treeinfo.txt",
    # chm_base
    "cbhf" => basefolder*"/input/test_chm_base.asc",
    # swr
    "swrf" => basefolder*"/input/test_swr.txt"
    )

    return dat_in, par_in


end
