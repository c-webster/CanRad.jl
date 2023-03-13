using Test, DelimitedFiles, NCDatasets, CanRad

println("Packages compiled")

cd("../testset")

datdir = pwd()

ptsf  = "test_pts.txt"

bmk_dd = readdlm(joinpath(datdir,"real_HPs","info.txt"),header=true)[1][:,4]

@testset "all_models" begin

    for model in ["T2R","C2R","L2R"]

        setf       = joinpath(datdir,model*"_settings_test.jl")

        include(setf)

        pts = readdlm(ptsf)

        exdir      = joinpath(datdir,"Output_tests_"*model)

        if model == "T2R"

            dat_in, par_in = ter2rad_settings(datdir)
            println("Running T2R ...")
            ter2rad!(pts,dat_in,par_in,exdir,"T2R test")
            
        else

            if model == "C2R"

                dat_in, par_in = chm2rad_settings(datdir)
                println("Testing C2R ...")
                chm2rad!(pts,dat_in,par_in,exdir,"C2R test")

            elseif model == "L2R"
    
                dat_in, par_in = las2rad_settings(datdir)
                println("Testing L2R ...")
                # las2rad!(pts,dat_in,par_in,exdir,"L2R test")

            end

            make_SHIs(exdir) # make the images to visually compare if test fails

            ncds = Dataset(joinpath(exdir,"Output_testset.nc"))
            tst_dd = ncds["Vf_planar"][:].*0.01
            close(ncds)

            # using Plots
            # p1 = plot(bmk_dd,tst_dd,seriestype=:scatter)
            # plot!(0:1,0:1,seriestype=:line)

            # test rmse is <0.05
            @test sum(sqrt.((tst_dd .- bmk_dd).^2))/size(pts,1) < 0.05

        end

    end

end

cd("..")