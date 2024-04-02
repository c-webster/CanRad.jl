function las2rad_settings(basefolder::String)

    par_in = Dict(

        # terrain settings
        "terrain_highres"     => true, # include high high resolution local terrain (1-5 m)
        "highres_peri"        => 300, # radius for including high-res terrain !! should be float number!
        "terrain_lowres"      => true, # !MUST BE TRUE! include low resolution regional terrain (> 25 m)
        "lowres_peri"         => 10000, # radius around points to include lowres terrain
        "terrainmask_precalc" => true, # if the terrain mask has been pre-calculated by CanRad-T2R
        "calc_terrain"        => true, # calculate variables for terrain-only as well as forest
        "buildings"           => false, # include buildings as opague features (requires a building height model, bhmf, like the chm)
        "horizon_line"        => false, # use pre-calculated horizon line matrix (e.g. for compatibility with topo-downscaled swr data)

        # forest settings
        "trunks"         => true, # true or false
        "branches"       => true, # true or false
        "branch_spacing" => 0.1, # spacing interval for creating branches if left empty the default is 0.1
        "season"         => "winter", # "winter" or "nothing" see docs for more information about this feature
        "forest_peri"    => 100, # radius around points to include in image

        # calculation settings
        "calc_trans" => true,
        "calc_swr" => 1, # 0 = off; 1 = potential swr (atm_trans = 1); 2 = real swr (needs swrf in dat_in)

        # time step settings
        "t1"  => "01.10.2019 00:00:00", # start time stamp; format "dd.mm.yyyy HH:MM:SS"
        "t2"  => "21.06.2020 23:50:00", # end time stamp
        "tstep" => 2, # in minutes; note this is the time step of the saved data; the model is calculated at 2 minute intervals for a complete solar track

        # location settings
        "time_zone"   => 1, # relative to UTC (e.g CET = UTC+1)
        "coor_system" => "CH1903", # possible options "CH1903", "CH1903+" or e.g. "UTM 32 N" for UTM systems (spaces required)

        # image settings
        "image_height" => 0.25, #FLOAT NUMBER  !!!!!!!!
        "point_size"  => [30.0,10.0], # how big to make the points in the output image in pixels (image is 1000x1000) this needs to be calibrated for point cloud density

        # run settings
        "batch"       => false,
        "save_images" => true, # save the calculated images in a netcdf file
        "progress"    => true # report individual step progress in /ProgressLastPoint

    )

    # all spatial input data must be in the same coordinate system
    dat_in = Dict(
    # chm
    "dsmf" => basefolder*"/input/test.laz",
    # dtm
    "dtmf" => basefolder*"/input/test_dtm.tif",
    # dem
    "demf" => basefolder*"/input/dem_50m.tif",
    # terrain mask (if terrain_precalc = true)
    "terf" => basefolder*"/Output_tests_T2R/SHIs_testset.nc",
    # horizon line file
    "hlmf" => basefolder*"/Output_tests_T2R/HLM_testset.nc",
    # dbh
    "dbhf" => basefolder*"/input/test_treeinfo.txt",
    # ltc
    "ltcf" => basefolder*"/input/test_lastreeclass.laz",
    # swr
    "swrf" => basefolder*"/input/test_swr.txt"
    )

    return dat_in, par_in


end
