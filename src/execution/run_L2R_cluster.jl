
# bfdr = "M:/WORK_ALL/Data/LAS2Rad/"
# bfdr = "C:/Users/v1cwebst/GoogleDrive/Data/LAS2Rad/"
bfdr = "/home/v1cwebst/WORK_ALL/Data/LAS2Rad/"
# bfdr = "/exports/csce/eddie/geos/groups/geos_forsnow/LAS2Rad/"

wrkdir     = bfdr*"Scripts/"

basefolder = bfdr*"Data_Areas/LaretLarge"

setf       = bfdr*"Settings/LAS2Rad_Settings_LaretLarge.jl"


###############################################################################
### BEGIN

include(wrkdir*"LAS2Rad_PKG_check.jl")
include(wrkdir*"LAS2Rad_FunctionCall.jl")
include(wrkdir*"LAS2Rad.jl")

include(wrkdir*"extract_settings.jl")
include(setf)

inputs     = readdlm(ARGS[1])
segname    = inputs[2]

par_in, dat_in = LAS2Rad_Settings(basefolder,segname)

pts = readdlm(basefolder*"/Tiles/"*inputs[1])

exdir = basefolder*"/Output/"
if !ispath(exdir)
    mkdir(exdir)
end

LAS2Rad(pts,dat_in,par_in,exdir,ARGS[1])
