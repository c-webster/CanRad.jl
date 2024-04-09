function chm2rad_settings(basefolder="empty"::String)

    par_in = Dict(

        "special_implementation" => "none", # "none" or "swissrad" or "oshd"

        # terrain settings
        "terrain_highres"     => true, # include high high resolution local terrain (1-5 m)
        "highres_peri"        => 300, # radius for including high-res terrain 
        "terrain_lowres"      => true, # include low resolution regional terrain (> 25 m)
        "lowres_peri"         => 25000, # radius for low-res terrain (calculating regional horizon line)
        "terrainmask_precalc" => false, # if the terrain mask has been pre-calculated by CanRad
        "calc_terrain"        => true, # calculate variables for terrain-only as well as forest
        "buildings"           => false, # include buildings in terrain calculation - requires blmf file
        "hlm_precalc"         => false, # use pre-calculated horizon line matrix (e.g. for compatibility with topo-downscaled swr data)

        # forest settings
        "season"       => "both", # summer, winter, both
        "forest_type"  => "evergreen", # evergreen, deciduous, mixed
        "tree_species" => "needleleaf", # needleleaf, broadleaf, both
        "forest_peri"  => 100, # radius within which forest is included in SHI calculation
        "cbh"          => 2, # height of canopy base above ground (can also be input as raster for varied base height)

        # calculation settings
        "calc_trans" => true, # calculate time-varying transmissivity
        "calc_swr"   => 1, # 0 = off; 1 = potential swr (atm_trans = 1); 2 = real swr (needs swrf in dat_in)

        # time step settings
        "t1"    => "01.10.2019 00:00:00", # start time stamp; format "dd.mm.yyyy HH:MM:SS"
        "t2"    => "21.06.2020 23:50:00", # end time stamp
        "tstep" => 2, # in minutes

        # location settings
        "time_zone"   => 1, # relative to UTC (e.g CET = UTC+1)
        "coor_system" => "CH1903", # possible options "CH1903", "CH1903+" or e.g. "UTM 32 N" for UTM systems (spaces required)

        # image settings
        "image_height" => 0.25, #FLOAT NUMBER  !!!!!!!!

        # run settings
        "batch"       => false, # running in parallel or single process
        "save_images" => true, # save the calculated images in a netcdf file
        "save_hlm"    => true, # save the calculated terrain horizon line 
        "progress"    => true  # report individual step progress in /ProgressLastPoint

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
        # horizon line file
        "hlmf" => basefolder*"/Output_tests_T2R/HLM_testset.nc",
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
