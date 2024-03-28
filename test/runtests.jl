using DelimitedFiles, NCDatasets, CanRad, Test

println("Packages compiled")

# cd("../testset")

cd("testset")

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

            ncds = Dataset(joinpath(exdir,"Output_testset.nc"))
            tst_dd = ncds["svf_planar"][:].*0.01
            close(ncds)

            @test 14.2 < sum(tst_dd) < 14.3

        else

            if model == "C2R"

                dat_in, par_in = chm2rad_settings(datdir)
                println("Testing C2R ...")
                chm2rad!(pts,dat_in,par_in,exdir,"C2R test")

                ncds = Dataset(joinpath(exdir,"Output_testset.nc"))
                tst_dd = ncds["svf_planar_e"][:].*0.01
                close(ncds)

                # test rmse is <0.05
                @test sum(sqrt.((tst_dd .- bmk_dd).^2))/size(pts,1) < 0.05

            elseif model == "L2R"
    
                dat_in, par_in = las2rad_settings(datdir)
                println("Testing L2R ...")
                las2rad!(pts,dat_in,par_in,exdir,"L2R test")

                ncds = Dataset(joinpath(exdir,"Output_testset.nc"))
                tst_dd = ncds["svf_planar"][:].*0.01
                close(ncds)

                # test rmse is <0.05
                @test sum(sqrt.((tst_dd .- bmk_dd).^2))/size(pts,1) < 0.05
                
            end


        end

    end

end

cd("..")

for model in ["T2R","C2R","L2R"]
    rm(joinpath(datdir,"Output_tests_"*model),force=true,recursive=true)
end