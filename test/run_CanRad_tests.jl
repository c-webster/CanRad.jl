
basefolder = "D:/CanRad/Domains/TEST/"

model      = "LAS" #LAS or CHM

setf       = basefolder*model*"2Rad_Settings_test.jl"

ptsf       = joinpath(basefolder,"test_gridpts.txt")

exdir      = joinpath(basefolder,"Output_tests_readltclas")

benchmark  = false # checks Vf values are the same as benchmark values

###############################################################################
### BEGIN

using CanRad, DelimitedFiles

include(setf)

pts = readdlm(ptsf)

if model == "LAS"

    dat_in, par_in = LAS2Rad_Settings(basefolder)
    exdir = exdir*"_L2R"
    LAS2Rad(pts[1:8,:],dat_in,par_in,exdir)

elseif model == "CHM"

    dat_in, par_in = CHM2Rad_Settings(basefolder)
    exdir = exdir*"_C2R"
    CHM2Rad(pts,dat_in,par_in,exdir)

end

if benchmark

    using NCDatasets

    exdir_bm = joinpath(basefolder,"Output_tests_L2R")
    bmk_dd = Dataset(joinpath(exdir_bm,"BENCHMARK","Output_TEST.nc"))
    tst_dd = Dataset(joinpath(exdir,"Output_TEST.nc"))

    if sum(sum(tst_dd["Vf_planar"][:].*0.01 .- bmk_dd["Vf_weighted"][:])) < 0.001
        println("Test passed Vf check")
    else
        println("Test failed Vf check")
    end

     close(bmk_dd)
     close(tst_dd)

end
