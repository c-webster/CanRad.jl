
bdr          = "M:/WORK_ALL/Data/LAS2Rad/"

bdr          = "C:/Users/v1cwebst/GoogleDrive/Data/LAS2Rad/"

wrkdir       = bdr*"Scripts/"

basefolder   = bdr*"Data_Areas/LaretLarge"

segname      = "SegmentA"

setf         = bdr*"Settings/LAS2Rad_Settings_LaretLarge.jl"

lastoolspath = "C:/Workspace/LAStools/bin/"


###############################################################################
### BEGIN

include(wrkdir*"LAS2Rad_PKG_check.jl")
include(wrkdir*"LAS2Rad_FunctionCall.jl")
include(wrkdir*"LAS2Rad.jl")

include(wrkdir*"extract_settings.jl")
include(setf)

tsize = readdlm(basefolder*"/TileSize.txt")

par_in, dat_in = LAS2Rad_Settings(basefolder,segname)

limits = readdlm(basefolder*"/"*segname*"/"*segname*"_analysisarea.txt")

limx = collect(limits[1]:tsize[1]:limits[2])
limy = collect(limits[3]:tsize[1]:limits[4])

pts_all = readdlm(basefolder*"/"*segname*"/"*segname*"_gridpts.txt")

exdir = basefolder*"/"*segname*"/"*"Output"
if !ispath(exdir)
    mkdir(exdir)
end

# for x in eachindex(limx[1:end-1]), y in eachindex(limy[1:end-1])
x=1; y=1
    # create the tasks
    idx = (limx[x] .<= pts_all[:,1] .< limx[x+1]) .&
            (limy[y] .<= pts_all[:,2] .< limy[y+1])

    pts = pts_all[idx,:]

    start = time()

    LAS2Rad(pts,dat_in,par_in,exdir)

    elapsed = time() - start
    println("Segment elapsed time: ",round(elapsed1,digits=2)," seconds")

# end
