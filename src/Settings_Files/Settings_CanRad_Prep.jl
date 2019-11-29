function Settings_CanRad_Prep()

    set_in = Dict(

    "lastoolpath"   => "C:/Workspace/LASTools/bin/",

    "rpath"         => "C:/'Program Files'/R/R-3.5.3/bin/",

    "preppath"      => "C:/Users/webster/GDrive/Data/L2HEval/Scripts_Associated/L2HE_Prep_R.R",
    # full path to the prep R script (to be integrated in future versions)

    "largelidar"    => "C:/Users/webster/GDrive/Data/L2HEval/Data_Surface/DSM/FlowerSite.las",

    "outfolder"     => "D:/LAS2Rad/Data_Areas/FlowerSite",

    "lidartilesize" => 1000,

    "makeCHM"       => false,

    "chmfmt"        => ".tif",

    "lasfmt"        => ".laz",

    "chmstep"       => 0.25,

    "xlimits"       => (484350,484490),
    "ylimits"       => (7472265,7472405),

    "buffer"        => 100,

    "trunks"        => true,
    "biome"         => 17,

    "branches"      => true,

    "species"       => 1,
    # 0: NA -> uses horizontal branches
    # 1: Norway Spruce
    # 2: Douglas Fir
    # 3: Scott Pint
    # 4: Silver Birch

    "gengrid"       => true,
    "spacing"       => 1,

    "epsgstr"       => "25835",
    # 25835 = UTM35N for Finland
    # 2056  = Swiss Grid

    "batch"         => true,

    "sdims"         => (2,2),
    "tilesize"      => 25

    )

    return set_in

end
