
model      = "CHM" #LAS or CHM

setf       = "D:/LAS2Rad/Data_Areas/TEST/"*model*"2Rad_Settings_test.jl"

ptsf       = "D:/LAS2Rad/Data_Areas/TEST/test_gridpts.txt"

exdir      = "D:/LAS2Rad/Data_Areas/TEST/Output_tests"

benchmark  = false # checks Vf values are the same as benchmark values

###############################################################################
### BEGIN

using CanRad, DelimitedFiles

include(setf)

pts = readdlm(ptsf)

if model == "LAS"

    par_in, dat_in = LAS2Rad_Settings()
    exdir = exdir*"_L2R"
    LAS2Rad(pts,dat_in,par_in,exdir)

elseif model == "CHM"

    par_in, dat_in = CHM2Rad_Settings()
    exdir = exdir*"_C2R"
    CHM2Rad(pts,dat_in,par_in,exdir)

end

if benchmark

    using NCDatasets

    bmk_dd = Dataset(joinpath(exdir,"BENCHMARK","Output_TEST.nc"))
    tst_dd = Dataset(joinpath(exdir,"Output_TEST.nc"))

    if sum(sum(tst_dd["Vf_weighted"][:].*0.01 .- bmk_dd["Vf_weighted"][:])) < 0.001
        println("Test passed Vf check")
    else
        println("Test failed Vf check")
    end

     close(bmk_dd)
     close(tst_dd)

end