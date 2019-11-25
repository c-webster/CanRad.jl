
wrkdir     = "M:/WORK_ALL/Data/LAS2Rad/Scripts/"

basefolder = "M:/WORK_ALL/Data/LAS2Rad/Data_Areas/LaretLarge"

using DelimitedFiles, Formatting

tsize = readdlm(basefolder*"/TileSize.txt")

segs = readdir(basefolder*"/Segments")

if !ispath(basefolder*"/Tiles")
    mkdir(basefolder*"/Tiles")
end

if !ispath(basefolder*"/Input")
    mkdir(basefolder*"/Input")
end

tdx = 0

for sdx in eachindex(segs)

    limits = readdlm(basefolder*"/Segments"*"/"*segs[sdx]*"/"*segs[sdx]*"_analysisarea.txt")

    limx = collect(limits[1]:tsize[1]:limits[2])
    limy = collect(limits[3]:tsize[1]:limits[4])

    pts_all = readdlm(basefolder*"/Segments"*"/"*segs[sdx]*"/"*segs[sdx]*"_gridpts.txt")

    for x in eachindex(limx[1:end-1]), y in eachindex(limy[1:end-1])
        # create the tasks
        idx = (limx[x] .<= pts_all[:,1] .< limx[x+1]) .&
                (limy[y] .<= pts_all[:,2] .< limy[y+1])

        pts = pts_all[idx,:]

        # create txtfile (1:sdx.txt)
        global tdx += 1
        exfname = sprintf1.("%03.$(0)f", tdx)*".txt"
        writedlm(basefolder*"/Tiles"*"/"*exfname,pts)

        exfnametask = "task."*sprintf1.("%03.$(0)f", tdx)
        writedlm(basefolder*"/Input/"*exfnametask,(exfname,segs[sdx]))

    end

end
