
println("new")

using Distributed

addprocs(Sys.CPU_THREADS)

@everywhere basefolder   = "D:/LAS2Rad/Data_Areas/Ukraine"

@everywhere ptsf         = basefolder*"/Ukraine_gridpts.txt"

@everywhere include("C:\\Users\\webster\\.julia\\dev\\CanRad\\src\\Settings_Files/LAS2Rad_Settings_Ukraine.jl")

###############################################################################
### BEGIN

@everywhere using CanRad, DelimitedFiles
            using Formatting

tsize = readdlm(basefolder*"/TileSize.txt")

segs = readdir(basefolder*"/Segments")

if !ispath(basefolder*"/Tiles")
    mkdir(basefolder*"/Tiles")
end

pts_all = Int.(readdlm(ptsf))

tdx = 0

@everywhere ptsfname     = String[]
@everywhere inputsegname = String[]

for sdx in eachindex(segs)

    limits = readdlm(basefolder*"/Segments"*"/"*segs[sdx]*"/"*segs[sdx]*"_analysisarea.txt")

    limx = Int.(collect(limits[1]:tsize[1]:limits[2]))
    limy = Int.(collect(limits[3]:tsize[1]:limits[4]))

    for x in eachindex(limx[1:end-1]), y in eachindex(limy[1:end-1])
        # create the tasks
        idx = (limx[x] .<= pts_all[:,1] .< limx[x+1]) .&
                (limy[y] .<= pts_all[:,2] .< limy[y+1])

        if sum(idx) > 0
            global tdx += 1
            push!(ptsfname,sprintf1.("%03.$(0)f", tdx)*".txt")
            # writedlm(basefolder*"/Tiles"*"/"*ptsfname[end],Int.(pts_all[idx,:]))

            push!(inputsegname,segs[sdx])
        end
    end

end

@everywhere exdir = basefolder*"/Output/"
if !ispath(exdir)
    mkdir(exdir)
end

println(basefolder)

@sync @distributed for x = 1:1:5# in eachindex(ptsfname)

        par_in, dat_in = LAS2Rad_Settings(basefolder,inputsegname[tdx])

        pts = readdlm(basefolder*"/Tiles"*"/"*ptsfname[tdx])

        LAS2Rad(pts,dat_in,par_in,exdir)

end
