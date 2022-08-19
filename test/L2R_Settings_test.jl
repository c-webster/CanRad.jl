function LAS2Rad_Settings(basefolder::String)

    par_in = Dict(
    "trunks" => true, # true or false
    "branches" => true, # true or false
    "branch_spacing" => 0.1, # spacing interval for creating branches if left empty the default is 0.1

    # terrain settings
    "terrain_highres" => true, # include high high resolution local terrain (1-5 m)
    "terrain_lowres" => true, # include low resolution regional terrain (> 25 m)
    "terrain_tile" => true, # calculate terrain_lowres once for the whole tile

    "terrain_peri" => 10000, # radius around points to include lowres terrain

    "image_height" => 2.0, #FLOAT NUMBER  !!!!!!!!

    "calc_trans" => true,
    "calc_swr" => 1, # 0 = off; 1 = potential swr (atm_trans = 1); 2 = real swr (needs swrf in dat_in)

    "t1"  => "01.10.2019 00:00:00", # start time stamp; format "dd.mm.yyyy HH:MM:SS"
    "t2"  => "21.06.2020 23:50:00", # end time stamp
    "tstep" => 2, # in minutes
    # note this is the time step of the saved data; the model is calculated at 2 minute intervals for a complete solar track

    "time_zone" => 1,
    "coor_system" => "CH1903",
    # possible options "CH1903", "CH1903+" or e.g. "UTM 32 N" for UTM systems (spaces required)

    "save_images" => true, # save the calculated images in a netcdf file

    "surf_peri" => 100, # radius around points to include in image

    "tolerance" => 1.5, # how big to make the points in the output image (this needs to be calibrated for point cloud density)

    "batch" => false,

    "buildings" => false, # include buildings as opague features (requires a building height model, bhmf, like the chm)

    "progress" => true # report individual step progress in /ProgressLastPoint
    )

    # all spatial input data must be in the same coordinate system
    dat_in = Dict(
    # chm
    "dsmf" => basefolder*"/data/test.laz",
    # dtm
    "dtmf" => basefolder*"/data/test_dtm.tif",
    # dem
    "demf" => basefolder*"/data/dem_50m.tif",
    # dbh
    "dbhf" => basefolder*"/data/test_treeinfo.txt",
    # ltc
    "ltcf" => basefolder*"/data/test_lastreeclass.laz",
    # swr
    "swrf" => basefolder*"/data/test_swr.txt"
    )

    return dat_in, par_in


end
