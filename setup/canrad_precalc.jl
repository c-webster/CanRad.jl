"""
canrad_precalc.jl

by Clare Webster, 2024

uses pycrown (https://github.com/manaakiwhenua/pycrown) to retrieve individual tree crown information
   calculates specific datasets required by the CanRad models L2R or C2R

example uses the same domain as the test data in CanRad (see ../CanRad/testset/)

"""

#################################################################################

# requirements

# 1. pycrown installed within a conda environment (run this script within this environment)
# 2. add PolygonOps, StaticArrays, DelimitedFiles, Statistics, SpatialFileIO, PyCall
# 3. add https://github.com/c-webster/SpatialFileIO.jl
# 4. add GeoJSON@0.5.1 !!! cannot be higher or will fail
# 5. add Conda; using Conda; Conda.add("numpy")


#####################################################################################
# user settings

wrkdir = "../setup/example"
shp2json = "../setup/shp2json.py"
lastoolspath = "../software/LAStools/bin" # not provided

model = "L2R" # can be L2R (for trunks or branches) or C2R (for crown-specific leaf area volume density)

pycrownscript = "example.py"

forest_type = "TCG" # for dbh calculations using relationships from Jucker et al. 2017 (https://doi.org/10.1111/gcb.13388)
# TCG = temperate coniferous gymnosperm (Palearctic) / e.g. Spruce Switzerland
# BG  = boreal gymnosperm (Palearctic) / e.g. Finland Pine
# TMG = temperate mixed gymnosperm (Palearctic) / e.g. Beech

model = "C2R"
# turns on calculation of trunks or lavd etc.

tree_species   = "Spruce" # for calculation of leaf area for C2R
# "Spruce", "Beech", "Scot Pine"


####################################################################################
####################################################################################
####################################################################################


pycrowndir = joinpath(wrkdir,"pycrown")
inpycrown  = joinpath(pycrowndir,"data")
outpycrown = joinpath(pycrowndir,"result")
fname      = joinpath(outpycrown,"tree")
chmf       = joinpath(inpycrown,"CHM.tif")

runpycrown = joinpath(pycrowndir,pycrownscript)

# PART 1

run(`$"python" $runpycrown $pycrowndir`)

run(`$"python" $shp2json $pycrowndir`)


# PART 2

using PolygonOps, StaticArrays, GeoJSON, DelimitedFiles, Statistics, SpatialFileIO, PyCall, ArchGDAL
const GJ = GeoJSON
const AG = ArchGDAL
np = pyimport("numpy")

if model == "C2R"

    chm_x, chm_y, chm_z = read_griddata(chmf,true,false)

    chm_lavd = fill(0.0,size(vec_chm_z))
    # chm_base = fill(0.0,size(chm_z))

end

tops = GJ.read(read(fname*"_location_top_cor.geojson"))
tops_dict = GJ.geo2dict(tops)

crowns = GJ.read(read(fname*"_crown_poly_raster.geojson"))
crowns_dict = GJ.geo2dict(crowns)

tx = Vector{Float64}(); ty = Vector{Float64}(); tz = Vector{Float64}(); tn = Vector{Float64}();

for x = 1:size(tops_dict["features"],1)
    append!(tx,(tops_dict["features"][x]["geometry"]["coordinates"])[1]) # tree x coordinate
    append!(ty,(tops_dict["features"][x]["geometry"]["coordinates"])[2]) # tree y coordinate
    append!(tz,tops_dict["features"][x]["properties"]["TH"])             # tree height
    append!(tn,Base.parse(Int64,tops_dict["features"][x]["id"]))         # tree number from pycrown calcs
end

# x, y, height + number of tree crowns (new numbers determined numerically)
cx =  Vector{Float64}(); cy =  Vector{Float64}(); cz =  Vector{Float64}(); cn = Vector{Int64}();
# cd_m = Vector{Float64}(); # crown diameter in m
dbh =  Vector{Float64}();

for x = 1:size(crowns_dict["features"],1) # loop through the crowns

    # global cx,cy,cz,cd_m
    out_x = Array(hcat(crowns_dict["features"][x]["geometry"]["coordinates"][1][:]...)')[:,1]
    out_y = Array(hcat(crowns_dict["features"][x]["geometry"]["coordinates"][1][:]...)')[:,2]

    polygon = SVector.(push!(out_x,out_x[end]),push!(out_y,out_y[end]))
    # points  = vec(SVector.((floor(minimum(out_x)):0.05:ceil(maximum(out_x)))',(floor(minimum(out_y)):0.05:ceil(maximum(out_y)))))

    # find which tree top is within the crown polygon
    tops = vec(SVector.(tx,ty))
    inside_tops  = [PolygonOps.inpolygon(p, polygon; in=true, on=true, out=false) for p in tops]

    dx = findall(inside_tops.==1)

    if ~isempty(dx)
        
        append!(cx,tx[dx[1]]); append!(cy,ty[dx[1]]); append!(cz,tz[dx[1]])
        cz_t = tz[dx[1]]

        ca_m2 = 0.5*np.abs(np.dot(out_x,np.roll(out_y,1))-np.dot(out_y,np.roll(out_x,1)));
        cd_m  = 2*sqrt(ca_m2/pi)

        append!(cn,Base.parse(Int64,tops_dict["features"][x]["id"]))

        # equations taken from Jucker et al., 2017
        if forest_type == "TCG"  # temperate coniferous gymnosperm (Palearctic)
            dbh_t = 0.974 * (cz_t*cd_m)^0.748
        elseif forest_type == "BG"  # boreal gymnosperm (Palearctic)
            dbh_t = 1.43  * (cz_t*cd_m)^0.649
        elseif forest_type == "TMG" # temperate mixed gymnosperm (Palearctic)
            dbh_t = 1.001 * (cz_t*cd_m)^0.73
        else
            error("forest type unknown")
        end

        append!(dbh,dbh_t)

    end

    if model == "C2R"

        if tree_species == "Spruce"
            LA = exp((-8.31)+(2.61*(log(dbh_t*10)))+(-0.07*cz_t)) # norway spruce (Goude et al. 2019)
        elseif tree_species == "Scot Pine"
            LA = 0.0091 * (dbh_t^2.363) * (cz[x]^0.106) # scot pine
        elseif tree_species == "Beech"
            LA = 8.56 + (0.0286*(dbh_t^2.623)) # DOI 10.1007/s10342-009-0345-8
        else
            println("tree species not defined")
        end

        cv = 0.8 * pi * ((cd_m/2)^2) * (cz_t*0.75) # norway spruce - crown height 80% of total tree height; scot pine = 75%

        dx1 = findall((minimum(out_x) .<= chm_x .<= maximum(out_x)) .& (minimum(out_y) .<= chm_y .<= maximum(out_y)))

        if sum(dx1) > 0 
            chm = SVector.(chm_x[dx1],chm_y[dx1])
            dx2  = [PolygonOps.inpolygon(p, polygon; in=true, on=true, out=false) for p in chm]

            chm_lavd[dx1[dx2]] .= LA/cv
            # chm_base[dx1] .= cz_t * 0.25 ###
        end

    end
    
end

if model == "L2R"

    lasf   = joinpath(wrkdir,filter(x->endswith(x,".laz"),readdir(wrkdir))[1])

    open(lasf[1:end-4]*"_treeinfo.txt";write=true) do f
        write(f,"treeptsX\ttreeptsY\ttreeptsH_m\ttreeptsdbh_cm\ttree_num\n")
        writedlm(f,round.(hcat(cx,cy,cz,dbh,cn),digits=3))
    end

    lasfile    = fname*"s.las"
    if isfile(lasfile) && (@isdefined(lastoolspath) && ispath(lastoolspath))
        if Sys.iswindows(); wine = [];
        elseif Sys.islinux(); wine = "wine"; end
        lasgroundcmd  = joinpath(lastoolspath,"lasground_new.exe")
        lazfile    = lasf[1:end-4]*"_lastreeclass.laz"
        run(`$wine $lasgroundcmd $"-i" $lasfile $"-wilderness" $"-o" $lazfile`)
    else
        @warn "finished with prep but still need to create *_lastreeclass.laz. See documentation"
    end

elseif model == "C2R"

    gt, crs, nrows, ncols = AG.read(chmf) do egdat
        gt = AG.getgeotransform(egdat)
        crs = AG.getproj(egdat)
        nrows = Int64(AG.height(egdat))
        ncols = Int64(AG.width(egdat))
        return gt, crs, nrows, ncols
    end

    chm_lavd[isnan.(chm_lavd)] .= 0
    chm_lavd[(chm_lavd.==0) .& (chm_z .> 2)] .= mean(chm_lavd[chm_lavd .> 0])
    chm_lavd[chm_z .<= 2]      .= 0
    chm_lavd[isnan.(chm_z)]    .= 0

    lavdf = joinpath(pycrowndir,"lavd.tif")
    lavd_out = Float64.(transpose(reverse(reshape(chm_lavd,nrows,ncols),dims=1)))

    AG.create(
        lavdf,driver = AG.getdriver("GTiff"),
        width=Int(ncols),height=Int(nrows),nbands=1,dtype=Float64
    )do dataset
        AG.write!(dataset,lavd_out, 1)
        AG.setgeotransform!(dataset, gt)
        AG.setproj!(dataset, crs)
    end

end
