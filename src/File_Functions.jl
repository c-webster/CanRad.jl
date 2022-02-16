function loaddbh(fname::String,limits::Array{Float64,2},peri::Int64)
    dat, _ = readdlm(fname,'\t',header=true)

    dat_x, dat_y, dat_z, rmdx = clipdat(dat[:,1],dat[:,2],dat[:,3],limits,peri)

    dat_r = deleteat!(dat[:,4],rmdx) # diameter at breast height
    dat_r /= 2
    dat_r /= 100

    # dat_h = deleteat!(dat[:,5],rmdx) # height

    if size(dat,2) == 5
        dat_tc = deleteat!(dat[:,5],rmdx)
    else
        dat_tc = fill(NaN,size(dat_x))
    end

    return dat_x, dat_y, dat_z, dat_r, dat_tc
end

function loadltc_txt(fname::String,limits::Array{Float64,2},peri::Int64)
    ltc, _ = readdlm(fname,'\t',header=true)
    replace!(ltc, -9999=>NaN)
    _, _, _, rmdx = clipdat(ltc[:,1],ltc[:,2],ltc[:,3],limits,peri)
    ltc = ltc[setdiff(1:end, findall(rmdx.==1)), :]
    return ltc
end



function loadltc_laz(fname::String,limits::Array{Float64,2},peri::Int64,
            dbh_x::Array{Float64,1},dbh_y::Array{Float64,1},dbh_e::Array{Float64,1},
            lastc::Array{Float64,1})

        # calculate random branch angle for each tree
        ang = rand(Uniform(60,100),size(lastc,1)) # Norway Spruce

        header, ltcdat = LazIO.load(fname)
        ltc = fill(NaN,size(ltcdat,1),8)
        ltc_c = fill(NaN,size(ltcdat,1),1)
        for dx in eachindex(ltcdat)
            ltc_c[dx] = trunc(Int,float(classification(ltcdat[dx])))
            if ltc_c[dx] == 2
                continue
            else
                ltc[dx,1] = xcoord(ltcdat[dx],header) #lasx
                ltc[dx,2] = ycoord(ltcdat[dx],header) #lasy
                ltc[dx,3] = zcoord(ltcdat[dx],header) #lasx
                ltc[dx,7] = Int64.(LasIO.intensity(ltcdat[dx])) # treeclass
                tree_dx = (lastc .== ltc[dx,7])
                if sum(tree_dx) > 0
                    ltc[dx,4] = dbh_x[tree_dx][1] #treex
                    ltc[dx,5] = dbh_y[tree_dx][1] #treey
                    ltc[dx,6] = dist(ltc[dx,1],ltc[dx,2],ltc[dx,4],ltc[dx,5])
                    ltc[dx,8] = ang[tree_dx][1]
                end
            end
        end
        ltc = ltc[setdiff(1:end, findall(isnan.(ltc[:,4]).==1)), :]
        _, _, _, rmdx = clipdat(ltc[:,1],ltc[:,2],ltc[:,3],limits,peri)
        ltc = ltc[setdiff(1:end, findall(rmdx.==1)), :]
    return ltc
end

function load_hlm(hlmf::String,taskID::String)

    xllcorner = parse(Int,(split(taskID,"_")[2]))[1]
    yllcorner = parse(Int,(split(taskID,"_")[3]))[1]

    if xllcorner % 250 .!== 0.0
        xllcorner = 250 * round(xllcorner/250)
    end
    if yllcorner % 250 .!== 0.0
        yllcorner = 250 * round(yllcorner/250)
    end

    hlmds = NCDataset(hlmf,"r")
    xdx = findall(hlmds["easting"][:,:] .== xllcorner)
    ydx = findall(hlmds["northing"][:,:] .== yllcorner)

    # phi = vec(hlmds["phi"][:])
    phi = (-pi:(pi - -pi)/89:pi) .- pi/2
    tht = (vec(hlmds["tht"][ydx,xdx,:]))
    close(hlmds)

    rphi = collect(-pi:pi/1080:pi)  .- pi/2
    rtht = Array{Float64,1}(undef,size(rphi,1))

    rtht = LinearInterpolation(phi,vec(tht))(rphi)

    pol_phi, pol_tht = fillterrain(rphi,rtht,0.0)

    return pol2cart(pol_phi,pol_tht)

end


function createfiles(outdir::String,outstr::String,pts::Array{Float64,2},calc_trans::Bool,calc_swr::Int64,
            append_file::Bool,loc_time=nothing)

    outfile  = outdir*"/Output_"*outstr*".nc"

    # append_file option removed September 2021 (eventually should be included)

    ds = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(ds,"Coordinates",size(pts,1))
    defVar(ds,"easting",pts[:,1],("Coordinates",))
    defVar(ds,"northing",pts[:,2],("Coordinates",))

    Vf_planar = defVar(ds,"Vf_planar",Int32,("Coordinates",),deflatelevel=5,
                            attrib=["scale_factor"=>0.01, "comments" =>
                            "perspective of a horizontal flat uplooking surface;
                            zenith rings weighted by surface area projected onto a horizontal flat surface",])
    Vf_hemi     = defVar(ds,"Vf_hemi",Int32,("Coordinates",),deflatelevel=5,
                            attrib=["scale_factor"=>0.01, "comments" =>
                            "perspective of hemipherically shaped surface or plant;
                            zenith rings weighted by surface area on the hemisphere",])

    if calc_trans
        defDim(ds,"datetime",size(loc_time,1))
        defVar(ds,"datetime",loc_time,("datetime",))

        for_tau = defVar(ds,"Forest_Transmissivity",Int32,("datetime","Coordinates",),
                            deflatelevel=5,attrib=["scale_factor"=>0.01,])

        if calc_swr > 0
            swr_tot = defVar(ds,"SWR_total",Int32,("datetime","Coordinates",),
                                deflatelevel=5,fillvalue = Int32(-9999),
                                attrib=["units"=>"Watts per metre squared", "comments" =>
                                "total incoming shortwave radiation (diffuse + direct)"])
            swr_dir = defVar(ds,"SWR_direct",Int32,("datetime","Coordinates",),
                                deflatelevel=5,fillvalue = Int32(-9999),
                                attrib=["units"=>"Watts per metre squared",])
            return swr_tot, swr_dir, for_tau, Vf_planar, Vf_hemi, ds
        else
            return for_tau, Vf_planar, Vf_hemi, ds
        end
    else
        return Vf_planar, Vf_hemi, ds
    end

end


function create_exmat(outdir::String,outstr::String,pts::Array{Float64,2},g_img::Array{Int64,2},append_file::Bool)

    outfile  = outdir*"/SHIs_"*outstr*".nc"

    # append_file option removed September 2021 (eventually should be included)

        images   = NCDataset(outfile,"c",format=:netcdf4_classic)

        dims     = size(g_img)

        defDim(images,"img_x",size(g_img,1))
        defDim(images,"img_y",size(g_img,2))
        defDim(images,"Coordinates",size(pts,1))

        SHIs = defVar(images,"SHI",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
        defVar(images,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
        defVar(images,"northing",pts[:,2],("Coordinates",),deflatelevel=1)

    return SHIs, images

end

function write_metadata(exdir::String,dat_in::Dict,par_in::Dict)

    open(exdir*"_metadata.txt";write=true) do f
        write(f,"DateTime\n")
        write(f,string(now())*"\n")
        write(f,"\nInput Files\n")
        for key in sort(collect(keys(dat_in)))
            write(f,"$key => $(dat_in[key])"*"\n")
        end
        write(f,"\nInput Parameters\n")
        for k in sort(collect(keys(par_in)))
            write(f,"$k => $(par_in[k])"*"\n")
        end
        write(f,"\nModel Versions\n")
        # writedlm(f,Pkg.installed())
        write(f,string("CanRad"*" v"*string(get_pkg_version("CanRad"))*"\n"))
        write(f,string("SpatialFileIO"*" v"*string(get_pkg_version("SpatialFileIO"))*"\n"))
    end

end

get_pkg_version(name::AbstractString) =

    @chain Pkg.dependencies() begin
       values
       [x for x in _ if x.name == name]
       only
       _.version

end
