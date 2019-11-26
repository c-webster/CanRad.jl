function Settings_makeDTM()

    set_in = Dict(

    "largelas" => "D:/Data/LidarStuff/Datasets/Fortress_Ridge/Fortress_Ridge_Lidar.las",

    "clip" => false, # true or false; clip the lidar to smaller limits (if true, set x and y limits below)

    "xlimits"       => (625450,625850),
    "ylimits"       => (5632700,5633250),

    "buffer"        => 300,

    "ground" => true, # existing lasground from tile (if true, clip = 0)

    "intm" => 1, # sample interval for DTM in same units as .las

    "wrkdir" => "C:/Users/webster/GDrive/Data/LidarStuff/temp",

    "lpath"  => "C:/Workspace/LASTools/bin",

    "exfile" => "C:/Users/webster"

    )

    return set_in




end
