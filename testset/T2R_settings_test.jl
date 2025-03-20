function ter2rad_settings(basefolder="empty"::String)

    par_in = Dict(

        # terrain settings
        "terrain_highres" => true, # include high high resolution local terrain (1-5 m)
        "highres_peri"    => 300, # radius for including high-res terrain 
        "terrain_lowres"  => true, # include low resolution regional terrain (> 25 m)
        "lowres_peri"     => 25000, # radius around points to include lowres terrain
        "buildings"       => false, # include buildings as opague features (requires a building height model, bhmf, like the chm)

        # calculation settings
        "calc_trans" => true,
        "calc_swr" => 1, # 0 = off; 1 = potential swr (atm_trans = 1); 2 = real swr (needs swrf in dat_in)

        # time step settings
        "t1"  => "01.10.2019 00:00:00", # start time stamp; format "dd.mm.yyyy HH:MM:SS"
        "t2"  => "21.06.2020 23:50:00", # end time stamp
        "tstep" => 2, # in minutes; note this is the time step of the saved data; the model is calculated at 2 minute intervals for a complete solar track

        # location settings
        "time_zone"   => 1, # relative to UTC (e.g CET = UTC+1)
        "epsg_code" => 21781, # epsg code of dataset (as integer) 

        # image settings
        "image_height" => 2.0, #FLOAT NUMBER  !!!!!!!!

        # run settings
        "batch"        => false, # running in parallel or single process
        "save_images"  => true, # save the calculated images in a netcdf file
        "make_pngs"    => true, # create .png files from the SHIs
        "save_horizon" => true, # save the calculated terrain horizon line 
        "progress"     => true  # report individual step progress in /ProgressLastPoint

    )

    # all spatial input data must be in the same coordinate system
    dat_in = Dict(

        # dtm
        "dtmf" => basefolder*"/input/test_dtm.tif",
        # dem
        "demf" => basefolder*"/input/dem_50m.tif",
        # swr
        "swrf" => basefolder*"/input/test_swr.txt"

    )

    return dat_in, par_in


end
