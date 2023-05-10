function ter2rad_settings(basefolder="empty"::String)

    par_in = Dict(

    # terrain settings
    "terrain_highres" => true, # include high high resolution local terrain (1-5 m)
    "terrain_lowres" => true, # include low resolution regional terrain (> 25 m)
    "terrain_tile" => false, # calculate terrain_lowres once for the centre point of the whole area

    "terrain_peri" => 10000, # radius around points to include lowres terrain

    "image_height" => 2.0, #FLOAT NUMBER  !!!!!!!!

    "calc_trans" => true, # calculate time-varying transmissivity
    "calc_swr" => 1, # 0 = off; 1 = potential swr (atm_trans = 1); 2 = real swr (needs swrf in dat_in)

    "t1"  => "01.10.2019 00:00:00", # start time stamp; format "dd.mm.yyyy HH:MM:SS"
    "t2"  => "21.06.2020 23:50:00", # end time stamp
    "tstep" => 2, # in minutes

    "time_zone" => 1, # relative to UTC (e.g CET = UTC+1)
    "coor_system" => "CH1903",
    # possible options "CH1903", "CH1903+" or e.g. "UTM 32 N" for UTM systems (spaces required)

    "save_hlm" => true, # save the calculated horizon line matrix 

    "save_images" => true, # to save the calculated terrain mask as an image

    "batch" => false,

    "buildings" => false, # include buildings as opague features (requires a building height model, bhmf, like the chm)

    "progress" => true  # report individual step progress in /ProgressLastPoint
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
