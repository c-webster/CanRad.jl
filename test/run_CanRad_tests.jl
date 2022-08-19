###############################################################################
### User Input

model      = "C2R" #L2R or C2R

###############################################################################
### Define file paths

using CanRad, DelimitedFiles

# get path to test dataset
basefolder = joinpath(join(split(pathof(CanRad),"/")[1:end-2],"/"),"test")

setf       = joinpath(basefolder,model*"_Settings_test.jl")

ptsf       = joinpath(basefolder,"test_gridpts.txt")

exdir      = joinpath(basefolder,"Output_tests_TEST")

###############################################################################
### BEGIN

include(setf)

pts = readdlm(ptsf)

exdir = exdir*"_"*model

if model == "L2R"

    dat_in, par_in = LAS2Rad_Settings(basefolder)
    dat_in, par_in = LAS2Rad(pts,dat_in,par_in,exdir)

elseif model == "C2R"

    dat_in, par_in = CHM2Rad_Settings(basefolder)
    dat_in, par_in = CHM2Rad(pts,dat_in,par_in,exdir)

end

write_metadata(exdir,dat_in,par_in)

if par_in["save_images"] == true

    make_SHIs(exdir)

end
