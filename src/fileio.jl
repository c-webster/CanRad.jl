function loaddbh(fname::String,limits::Vector{Float64},peri=0::Int64)
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

function loadltc_txt(fname::String,limits::Vector{Float64},peri::Int64)
    ltc, _ = readdlm(fname,'\t',header=true)
    replace!(ltc, -9999=>NaN)
    rmdx = clipdat(ltc[:,1],ltc[:,2],ltc[:,3],limits,peri)[4]
    ltc = ltc[setdiff(1:end, findall(rmdx.==1)), :]
    return ltc
end

function loadltc_laz(fname::String,limits::Vector{Float64},
    dbh_x::Vector{Float64},dbh_y::Vector{Float64},
    lastc::Vector{Float64},season::String=nothing)

    # calculate random branch angle for each tree
    ang = rand(Uniform(60,100),size(lastc,1)) # Norway Spruce

    header, ltcdat = LazIO.load(fname)

    dat = DataFrame(ltcdat)

    ignore = indexin(dat.intensity, lastc)

    dx = ((limits[1] .<= (dat.x .* header.x_scale .+ header.x_offset) .<= limits[2]) .&
        (limits[3] .<= (dat.y .* header.y_scale .+ header.y_offset) .<= limits[4])) .&
        (dat.raw_classification .!= 2) .& (ignore .!= nothing)

    ltc = DataFrame()
    
    ltc.x = dat.x[dx] .* header.x_scale .+ header.x_offset
    ltc.y = dat.y[dx] .* header.y_scale .+ header.y_offset
    ltc.z = dat.z[dx] .* header.z_scale .+ header.z_offset
    ltc.num = Int32.(dat.intensity[dx]) # tree number

    # get the tree class (0=evergreen, 1=deciduous)
    if season == "complete"
        ltc.cls = Int32.(dat.pt_src_id[dx]) # tree class
    else
        ltc.cls = Int32.(fill(0,(size(ltc.x)))) # leave empty
    end

    ltc.tx = fill(NaN,size(ltc.x,1))
    ltc.ty = fill(NaN,size(ltc.x,1))
    ltc.hd = fill(NaN,size(ltc.x,1))
    ltc.ang = fill(NaN,size(ltc.x,1))
    for tx in eachindex(lastc)
        t_dx = (lastc[tx] .== ltc.num)
        ltc.tx[t_dx] .= dbh_x[tx] # tree x
        ltc.ty[t_dx] .= dbh_y[tx] # tree y
        ltc.hd[t_dx] .= dist(ltc.x[t_dx],ltc.y[t_dx],dbh_x[tx],dbh_y[tx]) # horizontal distance between point and trunk
        ltc.ang[t_dx] .= ang[tx] # branch angle
    end

    # return ltc_x, ltc_y, ltc_z, tree_num, ltc_tx, ltc_ty, ltc_hd, ltc_ang
    return ltc

end

function createfiles(outdir::String,outstr::String,pts::Matrix{Float64},calc_trans::Bool,calc_swr::Int64,
    loc_time=nothing,time_zone=nothing)

    outfile  = joinpath(outdir,"Output_"*outstr*".nc")

    ds = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(ds,"Coordinates",size(pts,1))
    defVar(ds,"easting",pts[:,1],("Coordinates",))
    defVar(ds,"northing",pts[:,2],("Coordinates",))

    defVar(ds,"svf_planar",Int32,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
        attrib=[
            "long_name"=> "sky-view fraction planar",
            "comments" =>
            "perspective of a flat planar",])

    defVar(ds,"svf_hemi",Int32,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
        attrib=[
            "long_name"=> "sky-view fraction hemispherical",
            "comments" =>
            "perspective of hemipherically shaped surface or plant",])

    if calc_trans

        defDim(ds,"datetime",length(loc_time))

        if time_zone >= 0
            dt_comment = "time zone UTC+"*string(time_zone)
        elseif time_zone < 0
            dt_comment = "time zone UTC-"*string(time_zone)
        end

        defVar(ds,"datetime",loc_time,("datetime",),attrib=["comments" => dt_comment])

        defVar(ds,"for_trans",Int32,("datetime","Coordinates",),fillvalue = Int8(-99),
                    deflatelevel=5,                            
                    attrib=["long_name"=> "forest transmissivity",])

        if calc_swr > 0

            defVar(ds,"swr_total",Int32,("datetime","Coordinates",),
                        deflatelevel=5,fillvalue = Int32(-9999),
                        attrib=[
                            "units"=>"Watts per metre squared",
                            "comments" => " sum of calculated diffuse + direct radiation components"])
            defVar(ds,"swr_direct",Int32,("datetime","Coordinates",),
                        deflatelevel=5,fillvalue = Int32(-9999),
                        attrib=["units"=>"Watts per metre squared",])

        end

    end

    return ds

end

function createfiles(outdir::String,outstr::String,pts::Matrix{Float64},calc_trans::Bool,calc_swr::Int64,
    forest_type::String,season::String,calc_terrain::Bool,loc_time=nothing,time_zone=nothing)

    outfile  = joinpath(outdir,"Output_"*outstr*".nc")

    ds = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(ds,"Coordinates",size(pts,1))
    defVar(ds,"easting",pts[:,1],("Coordinates",))
    defVar(ds,"northing",pts[:,2],("Coordinates",))

    if forest_type == "evergreen"

        defVar(ds,"svf_planar_e",Int32,("Coordinates",),deflatelevel=5,fillvalue = Int32(-99),
            attrib=[
                "long_name"=> "sky-view fraction planar evergreen forest",
                "comments" =>
                "perspective of a flat planar surface
                calculated for 100% evergreen forest"
                    ,])
        defVar(ds,"svf_hemi_e",Int32,("Coordinates",),deflatelevel=5,fillvalue = Int32(-99),
            attrib=[
                "long_name"=> "sky-view fraction hemispherical evergreen forest",
                "comments" =>
                "perspective of hemipherically shaped surface or plant
                calculated for 100% evergreen forest",])

    elseif (forest_type == "deciduous") || (forest_type == "mixed")

        if (season == "summer") || (season == "both")

            defVar(ds,"svf_planar_s",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                attrib=[
                    "long_name"=> "sky-view fraction planar summer",
                    "comments" =>
                    "perspective of a flat planar surface
                    calculated for deciduous or mixed forests in summer canopy conditions"
                    ,])
            defVar(ds,"svf_hemi_s",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                attrib=[
                    "long_name"=> "sky-view fraction hemi summer",
                    "comments" =>
                    "perspective of hemipherically shaped surface or plant
                    calculated for deciduous or mixed forests in summer canopy conditions"
                    ,])

        end

        if (season == "winter") || (season == "both")

            defVar(ds,"svf_planar_w",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                attrib=[
                    "long_name"=> "sky-view fraction planar winter",
                    "comments" =>
                    "perspective of a flat planar surface
                    calculated for deciduous or mixed forests in winter canopy conditions"
                    ,])

            defVar(ds,"svf_hemi_w",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                attrib=[
                    "long_name"=> "sky-view fraction hemi winter",
                    "comments" =>
                    "perspective of hemipherically shaped surface or plant
                    calculated for deciduous or mixed forests in winter canopy conditions"
                    ,])

        end

    end

    if calc_terrain

        defVar(ds,"svf_planar_t",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
            attrib=[
                "long_name" => "sky-view fraction planar terrain",
                "comments" =>
                "perspective of a flat planar surface
                calcualated for terrain only"
                ,])
        defVar(ds,"svf_hemi_t",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
            attrib=[
                "long_name" => "sky-view fraction hemi terrain",
                "comments" =>
                "perspective of hemipherically shaped surface or plant
                calcualated for terrain only"
                ,])

    end


    if calc_trans

        defDim(ds,"datetime",length(loc_time))

        if time_zone >= 0
            dt_comment = "time zone UTC+"*string(time_zone)
        elseif time_zone < 0
            dt_comment = "time zone UTC-"*string(time_zone)
        end

        defVar(ds,"datetime",loc_time,("datetime",),attrib=["comments" => dt_comment])

        if forest_type == "evergreen"

            defVar(ds,"for_trans_e",Int32,("datetime","Coordinates",),fillvalue = Int8(-99),
                deflatelevel=5,
                attrib=[
                    "long_name"=> "forest transmissivity evergreen forest",
                    "comments" =>
                    "calculated for 100% evergreen forest"
                    ,])

        elseif (forest_type == "deciduous") || (forest_type == "mixed")

            if (season == "summer") || (season == "both")

                defVar(ds,"for_trans_s",Int32,("datetime","Coordinates",),fillvalue = Int8(-99),
                    deflatelevel=5,
                    attrib=[
                        "long_name"=> "forest transmissivity summer",
                        "comments" =>
                        "calculated for deciduous or mixed forests in summer canopy conditions"
                        ,])

            end

            if (season == "winter") || (season == "both")

                        defVar(ds,"for_trans_w",Int32,("datetime","Coordinates",),fillvalue = Int8(-99),
                            deflatelevel=5,
                            attrib=[
                                "long_name"=> "forest transmissivity winter",
                                "comments" =>
                                "calculated for deciduous or mixed forests in winter canopy conditions"
                                ,])

            end

        end

        if calc_terrain

            defVar(ds,"trans_t",Int8,("datetime","Coordinates",),fillvalue = Int8(-99),
                deflatelevel=5,
                attrib=[
                    "long_name" => "transmissivity terrain",
                    "comments"  => "calcualated for terrain only"
                    ,])
        
        end

        if calc_swr > 0

            if forest_type == "evergreen"

                defVar(ds,"swr_total_e",Int32,("datetime","Coordinates",),
                    deflatelevel=5,fillvalue = Int32(-9999),
                    attrib=[
                        "long_name" => "total incoming shortwave radiation evergreen forest",    
                        "units"=>"Watts per metre squared", 
                        "comments" =>
                        "calculated for 100% evergreen forest
                        sum of calculated diffuse + direct"
                        ,])

                defVar(ds,"swr_direct_e",Int32,("datetime","Coordinates",),
                    deflatelevel=5,fillvalue = Int32(-9999),
                    attrib=[
                        "long_name" => "direct incoming shortwave radiation evergreen forest",    
                        "units"=>"Watts per metre squared",
                        "comments" =>
                        "calculated for 100% evergreen forest"
                        ,])


            elseif (forest_type == "deciduous") || (forest_type == "mixed")

                if (season == "summer") || (season == "both")

                    defVar(ds,"swr_total_s",Int32,("datetime","Coordinates",),
                        deflatelevel=5,fillvalue = Int32(-9999),
                        attrib=[
                            "long_name" => "total incoming shortwave radiation summer",  
                            "units"=>"Watts per metre squared", 
                            "comments" =>
                            "calculated for deciduous or mixed forests in summer canopy conditions
                            sum of calculated diffuse + direct radiation components"
                            ,])

                    defVar(ds,"swr_direct_s",Int32,("datetime","Coordinates",),
                        deflatelevel=5,fillvalue = Int32(-9999),
                        attrib=[
                            "long_name" => "direct incoming shortwave radiation summer",  
                            "units"=>"Watts per metre squared",
                            "comments" =>
                            "calculated for deciduous or mixed forests in summer canopy conditions"
                            ,])

                end

                if (season == "winter") || (season == "both")

                    defVar(ds,"swr_total_w",Int32,("datetime","Coordinates",),
                        deflatelevel=5,fillvalue = Int32(-9999),
                        attrib=[
                            "long_name" => "total incoming shortwave radiation winter",
                            "units"=>"Watts per metre squared", 
                            "comments" =>
                            "calculated for deciduous or mixed forests in winter canopy conditions
                            sum of calculated diffuse + direct radiation components"
                            ,])

                    defVar(ds,"swr_direct_w",Int32,("datetime","Coordinates",),
                        deflatelevel=5,fillvalue = Int32(-9999),
                        attrib=[
                            "long_name" => "direct incoming shortwave radiation winter",    
                            "units"=>"Watts per metre squared",
                            "comments" =>
                            "calculated for deciduous or mixed forests in winter canopy conditions"
                            ,])

                end
            end

            if calc_terrain

            defVar(ds,"swr_total_t",Int32,("datetime","Coordinates",),
                deflatelevel=5,fillvalue = Int32(-9999),
                attrib=[
                    "long_name" => "total incoming shortwave radiation terrain",    
                    "units"=>"Watts per metre squared", 
                    "comments" =>
                    "calculated for terrain only
                    sum of calculated diffuse + direct radiation components"
                    ,])

            defVar(ds,"swr_direct_t",Int32,("datetime","Coordinates",),
                deflatelevel=5,fillvalue = Int32(-9999),
                attrib=[
                    "long_name" => "direct incoming shortwave radiation terrain",    
                    "units"=>"Watts per metre squared",
                    "comments" =>
                    "calculated for terrain only"
                    ,])

            end

        end

    end

	return ds

end

function create_exmat(outdir::String,outstr::String,pts::Matrix{Float64},g_img::Matrix{Int64})

    outfile  = joinpath(outdir,"SHIs_"*outstr*".nc")

    images = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(images,"img_x",size(g_img,1))
    defDim(images,"img_y",size(g_img,2))
    defDim(images,"Coordinates",size(pts,1))

    defVar(images,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
    defVar(images,"northing",pts[:,2],("Coordinates",),deflatelevel=1)
    defVar(images,"SHI",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)

    return images

end

function create_exmat(outdir::String,outstr::String,pts::Matrix{Float64},g_img::Matrix{Int64},forest_type::String,
                    season::String,calc_terrain::Bool)

    outfile  = joinpath(outdir,"SHIs_"*outstr*".nc")

    images = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(images,"img_x",size(g_img,1))
    defDim(images,"img_y",size(g_img,2))
    defDim(images,"Coordinates",size(pts,1))

    defVar(images,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
    defVar(images,"northing",pts[:,2],("Coordinates",),deflatelevel=1)

    if forest_type == "evergreen"
        defVar(images,"SHI_evergreen",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)

    elseif (forest_type == "deciduous") || (forest_type == "mixed")
    
        if (season == "summer") || (season == "both")
            defVar(images,"SHI_summer",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
        end
    
        if (season == "winter") || (season == "both")
            defVar(images,"SHI_winter",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
        end
    end
    
    if calc_terrain
        defVar(images,"SHI_terrain",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
    end

    return images

end

function create_exhlm(outdir::String,outstr::String,pts::Matrix{Float64},ter2rad::TER2RAD)

    outfile  = joinpath(outdir,"HLM_"*outstr*".nc")

    @unpack phi_bins, dx1, dx2 = ter2rad

    phi_bins = phi_bins[dx1:dx2-1]

    phi_bins .= ((pi/2) .- phi_bins) .* (180/pi);
    phi_bins[phi_bins.>=360] .= phi_bins[phi_bins.>=360] .- 360;              
    phi_bins[phi_bins.<0]    .= phi_bins[phi_bins.<0] .+ 360;         # phi now is the azimuth angle (phi) in degree, so that 0° is North, 90° is East, 180° is South, and 270° is West

    # p_srt = sortperm(phi_bins)

    hlm = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(hlm,"Coordinates",size(pts,1))
    defDim(hlm,"phi",size(phi_bins,1))

    defVar(hlm,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
    defVar(hlm,"northing",pts[:,2],("Coordinates",),deflatelevel=1)
    defVar(hlm,"phi",(round.(phi_bins)),("phi",),deflatelevel=1)
    # defVar(hlm,"phi_oshd",(round.(phi_bins))[p_srt],("phi",),deflatelevel=1)

    defVar(hlm,"tht",Int32,("phi","Coordinates",),fillvalue = Int32(-9999),
            deflatelevel=5,attrib=["scale_factor"=>0.01,"comments"=>"zenith angle of horizon line"],)
    # defVar(hlm,"tht_oshd",Int32,("phi","Coordinates",),fillvalue = Int32(-9999),
    #         deflatelevel=5,attrib=["scale_factor"=>0.01,"comments"=>"zenith angle of horizon line"],)

    return hlm 

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

function organise_outf(taskID::String,exdir::String,batch::Bool)

	if batch
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3]
        global outdir = exdir*"/"*outstr
    else
        outstr = splitpath(exdir)[end-1]
        global outdir = exdir
    end

    !ispath(outdir) && mkpath(outdir)

	# set start point within tile
	fx = readdir(outdir)[findall(startswith.(readdir(outdir),"Processing"))]
	!isempty(fx) && rm(joinpath(outdir,fx[1]))

	return outdir, outstr

end

function make_SHIs(datdir::String)

    files = readdir(datdir)
    fx = findall(startswith.(files,"SHIs") .& endswith.(files,".nc"))

    images = NCDataset(joinpath(datdir,files[fx][1]),"r")

    coords_x = images["easting"][:]
    coords_y = images["northing"][:]

    odir = joinpath(datdir,files[fx][1][1:end-3])

    !ispath(odir) && (mkpath(odir))

    fstr = "%07.$(2)f"

    for ix in eachindex(coords_x)

        outf = joinpath(odir,"SHI_"*cfmt.(fstr,coords_x[ix])*"_"*cfmt.(fstr,coords_y[ix])*".png")
        save(outf,colorview(Gray,float.(images["SHI"][:,:,ix])))

    end

end

function make_SHIs(datdir::String,forest_type::String,season::String,calc_terrain::Bool)

    files = readdir(datdir)
    fx = findall(startswith.(files,"SHIs") .& endswith.(files,".nc"))

    images = NCDataset(joinpath(datdir,files[fx][1]),"r")

    coords_x = images["easting"][:]
    coords_y = images["northing"][:]

    odir = joinpath(datdir,files[fx][1][1:end-3])

    if !ispath(odir); mkpath(odir); end

    if forest_type == "evergreen"
        odir_e = joinpath(odir,"terrain")
        !ispath(odir_e) && mkpath(odir_e)
    elseif (forest_type == "deciduous") || (forest_type == "mixed")
        if (season == "summer") || (season == "both")
            odir_s = joinpath(odir,"summer")
            !ispath(odir_s) && mkpath(odir_s)
        end
        if (season == "winter") || (season == "both")
            odir_w = joinpath(odir,"winter")
            !ispath(odir_w) && mkpath(odir_w)
        end
    end
    
    if calc_terrain
        odir_t = joinpath(odir,"terrain")
        !ispath(odir_t) && mkpath(odir_t)
    end

    fstr = "%07.$(2)f"

    for ix in eachindex(coords_x)

        if forest_type == "evergreen"
            outf = joinpath(odir,"SHI_"*cfmt.(fstr,coords_x[ix])*"_"*cfmt.(fstr,coords_y[ix])*"_evergreen.png")
            save(outf,colorview(Gray,float.(images["SHI_evergreen"][:,:,ix])))

        elseif (forest_type == "deciduous") || (forest_type == "mixed")
        
            if (season == "summer") || (season == "both")
                outf_s = joinpath(odir_s,"SHI_"*cfmt.(fstr,coords_x[ix])*"_"*cfmt.(fstr,coords_y[ix])*"_summer.png")
                save(outf_s,colorview(Gray,float.(images["SHI_summer"][:,:,ix])))
            end
        
            if (season == "winter") || (season == "both")
                outf_w = joinpath(odir_w,"SHI_"*cfmt.(fstr,coords_x[ix])*"_"*cfmt.(fstr,coords_y[ix])*"_winter.png")
                save(outf_w,colorview(Gray,float.(images["SHI_winter"][:,:,ix])))
            end
        end
        
        if calc_terrain
            outf_t = joinpath(odir_t,"SHI_"*cfmt.(fstr,coords_x[ix])*"_"*cfmt.(fstr,coords_y[ix])*"_terrain.png")
            save(outf_t,colorview(Gray,float.(images["SHI_terrain"][:,:,ix])))
        end

    end

    close(images)

end

function getterrainmask(canrad::CANRAD,terf::String,pts_x::Vector{Float64},pts_y::Vector{Float64})

    terrain_mask = Array{Int8}(undef,(size(canrad.mat2ev,1),size(canrad.mat2ev,2),size(pts_x,1)))

    tm_ds = NCDataset(terf)
    easting = tm_ds["easting"][:]
    northing = tm_ds["northing"][:]

    SHIs = tm_ds["SHI"][:,:,:]

    if (size(pts_x,1) .== size(easting,1)) && (sum(abs.(easting .- pts_x)) .== 0.0) && (sum(abs.(northing .- pts_y)) .== 0.0)
        terrain_mask .= tm_ds["SHI"][:,:,:]
    else
        for dx in eachindex(pts_x)
            tmx = (abs.(easting .- pts_x[dx]) .== 0.0) .& (abs.(northing .- pts_y[dx]) .== 0.0)
            terrain_mask[:,:,dx] = SHIs[:,:,tmx]
        end
    end

    close(tm_ds)

    return terrain_mask

end


function load_hlm(ter2rad::TER2RAD,hlmf::String,pts_x::Vector{Float64},pts_y::Vector{Float64})

    hlmds = NCDataset(hlmf)
    easting = hlmds["easting"][:]
    northing = hlmds["northing"][:]
    tht = identity.(hlmds["tht"][:]) ./ 100
    close(hlmds)

    # correct if coordinates aren't in the correct order
    if !(size(pts_x,1) .== size(easting,1)) || ((sum(abs.(easting .- pts_x)) .== 0.0) && (sum(abs.(northing .- pts_y))) .== 0.0)
        tht_new = Matrix{Float64}(undef,(size(tht,1),size(pts_x,1)))
        for dx in eachindex(pts_x)
            tmx = (abs.(easting .- pts_x[dx]) .== 0.0) .& (abs.(northing .- pts_y[dx]) .== 0.0)
            tht_new[:,dx] = tht[:,tmx]
            ###!!!! reorder tht to match coordinates of pts (not easting/northing of input file)
        end
    else
        tht_new = tht
    end

    @unpack rphi, phi_bins, dx1, dx2 = ter2rad
    phi = phi_bins[dx1:dx2]
    rtht = Matrix{Float64}(undef,(size(rphi,1),size(easting,1)))

    for dx in eachindex(pts_x)

        rtht[:,dx] = LinearInterpolation(phi,push!(tht_new[:,dx],tht_new[end,dx]))(rphi)

    end

    return rtht

end


function getterrainmask(canrad::CANRAD,terf::String,pts_x::Float64,pts_y::Float64)

    terrain_mask = Array{Int8}(undef,(size(canrad.mat2ev,1),size(canrad.mat2ev,2),1))

    tm_ds = NCDataset(terf)
    easting = tm_ds["easting"][:]
    northing = tm_ds["northing"][:]

    SHIs = tm_ds["SHI"][:,:,:]

    tmx = (abs.(easting .- pts_x) .== 0.0) .& (abs.(northing .- pts_y) .== 0.0)
    terrain_mask[:,:,dx] = SHIs[:,:,tmx]
    
    close(tm_ds)

    return terrain_mask

end

function load_hlm(ter2rad::TER2RAD,hlmf::String,pts_x::Float64,pts_y::Float64)

    hlmds = NCDataset(hlmf)
    easting = hlmds["easting"][:]
    northing = hlmds["northing"][:]
    tht = identity.(hlmds["tht"][:]) ./ 100
    close(hlmds)

    # correct if coordinates aren't in the correct order
    tmx = (abs.(easting .- pts_x) .== 0.0) .& (abs.(northing .- pts_y) .== 0.0)
    tht_new = tht[:,tmx]

    @unpack rphi, phi_bins, dx1, dx2 = ter2rad
    phi = phi_bins[dx1:dx2]
    rtht = LinearInterpolation(phi,push!(tht_new,tht_new[end]))(rphi)

    return rtht

end

function load_hlm_oshd(hlmf::String,taskID::String)

    error("load_hlm_oshd function not updated yet")

    # xllcorner = parse(Int,(split(taskID,"_")[2]))[1]
    # yllcorner = parse(Int,(split(taskID,"_")[3]))[1]

    # if xllcorner % 250 .!== 0.0
    #     xllcorner = 250 * round(xllcorner/250)
    # end
    # if yllcorner % 250 .!== 0.0
    #     yllcorner = 250 * round(yllcorner/250)
    # end

    # hlmds = NCDataset(hlmf,"r")
    # xdx = findall(hlmds["easting"][:,:] .== xllcorner)
    # ydx = findall(hlmds["northing"][:,:] .== yllcorner)

    # # phi = vec(hlmds["phi"][:])
    # phi = (-pi:(pi - -pi)/89:pi) .- pi/2
    # tht = (vec(hlmds["tht"][ydx,xdx,:]))
    # close(hlmds)

    # rphi = collect(-pi:pi/1080:pi)  .- pi/2
    # rtht = Array{Float64,1}(undef,size(rphi,1))

    # rtht = LinearInterpolation(phi,vec(tht))(rphi)

    # pol_phi, pol_tht = fillterrain(rphi,rtht,0.0)

    # return pol2cart(pol_phi,pol_tht)

end



function collate2tilefile(outdir::String,limits::Matrix{Int64},input::String,ptsx::Vector{Float64},
                            ptsy::Vector{Float64},exdir::String,hlmexdir::String,par_in::Dict{String, Any},
                            ptdx::Vector{Float64})

            if (!par_in["special_implementation"] == "swissrad") || (!par_in["special_implementation"] == "oshd")
                error("collate2tilefile() only works for swissrad implementation of CanRad")
            end

            !ispath(hlmexdir) && mkpath(hlmexdir)

            # combine tile output to one file
            tiles = readdir(outdir)

            xlim = collect(limits[1]:10:limits[2]-10)
            ylim = collect(limits[3]:10:limits[4]-10)

            # create the Radiation output dataset
            ds = NCDataset(joinpath(exdir,"Output_"*input*".nc"),"c",format=:netcdf4_classic)
            defDim(ds,"locX",size(xlim,1)); defDim(ds,"locY",size(ylim,1))
            defVar(ds,"Easting",Float32.(unique(ptsx)),("locX",)); defVar(ds,"Northing",Float32.(reverse(unique(ptsy))),("locY",))
            loc_time = NCDataset(joinpath(outdir,tiles[1],"Output_"*tiles[1]*".nc"))["datetime"][:]
            defDim(ds,"DateTime",size(loc_time,1))
            time_zone = par_in["time_zone"]
            if  time_zone >= 0
				dt_comment = "time zone UTC+"*string(time_zone)
			elseif time_zone < 0
				dt_comment = "time zone UTC-"*string(time_zone)*"; timestamp is for beginning of averaged period"
			end
			defVar(ds,"datetime",loc_time,("datetime",),attrib=["comments" => dt_comment])

            # create the horizon line output dataset
            hlm = NCDataset(joinpath(hlmexdir,"HLM_"*input*".nc"),"c",format=:netcdf4_classic)
            phi_bins = NCDataset(joinpath(outdir,tiles[1],"HLM_"*tiles[1]*".nc"))["phi"][:]
            defDim(hlm,"locX",size(xlim,1)); defDim(hlm,"locY",size(ylim,1))
            defVar(hlm,"Easting",Float32.(unique(ptsx)),("locX",)); defVar(hlm,"Northing",Float32.(reverse(unique(ptsy))),("locY",))
            defDim(hlm,"phi",size(phi_bins,1))
            defVar(hlm,"Horizon_Line_phi",Int16.(phi_bins),("phi",),deflatelevel=1,attrib=["units" =>
                "degrees clockwise from north"])
            defVar(hlm,"phi_sortdx",Int16.(NCDataset(joinpath(outdir,tiles[1],"HLM_"*tiles[1]*".nc"))["phi_sortdx"][:]),("phi",),deflatelevel=1,
                attrib=["comments" => "sort index for phi/tht variables for phi values to be sorted 0-360degrees",])

            for tx in tiles

                tds = NCDataset(joinpath(outdir,tx,"Output_"*tx*".nc"))
                pds = NCDataset(joinpath(outdir,tx,"HLM_"*tx*".nc"))

                if tx == tiles[1]

                    global east  = tds["easting"][:]
                    global north = tds["northing"][:]

                    global svf_p_w = tds["svf_planar_winter"][:]
                    global svf_p_s = tds["svf_planar_summer"][:]
                    global svf_p_t = tds["svf_planar_terrain"][:]

                    global svf_h_w = tds["svf_hemi_winter"][:]
                    global svf_h_s = tds["svf_hemi_summer"][:]
                    global svf_h_t = tds["svf_hemi_terrain"][:]

                    global for_tau_w = tds["Forest_Transmissivity_winter"][:]
                    global for_tau_s = tds["Forest_Transmissivity_summer"][:]
                    global for_tau_t = tds["Forest_Transmissivity_terrain"][:]

                    global tht   = pds["tht"][:]

                else

                    append!(east,tds["easting"][:])
                    append!(north,tds["northing"][:])

                    append!(svf_p_w,tds["svf_planar_winter"][:])
                    append!(svf_p_s,tds["svf_planar_summer"][:])
                    append!(svf_p_t,tds["svf_planar_terrain"][:])

                    append!(svf_h_w,tds["svf_hemi_winter"][:])
                    append!(svf_h_s,tds["svf_hemi_summer"][:])
                    append!(svf_h_t,tds["svf_hemi_terrain"][:])

                    for_tau_w = hcat(for_tau_w,tds["Forest_Transmissivity_winter"][:])
                    for_tau_s = hcat(for_tau_s,tds["Forest_Transmissivity_summer"][:])
                    for_tau_t = hcat(for_tau_t,tds["Forest_Transmissivity_terrain"][:])

                    tht  = hcat(tht,pds["tht"][:])

                end

                close(tds)
                close(pds)

            end

            B = hcat(east,north)
            p = sortperm(view.(Ref(B), 1:size(B,1), :));

            locs = (ptdx .> 0)

            tht_all = fill(-99,360,10000)
            tht_all[:,locs] .= tht[:,p]

            tht_newdim = fill(-99,360,100,100)
            for x in 1:1:size(tht_all,1)
                tht_newdim[x,:,:] = reverse(reshape(tht_all[x,:],100,100),dims=1)
            end

            defVar(hlm,"Horizon_Line_tht",Int16.(tht_newdim),("phi","locY","locX",),deflatelevel=1,attrib=["units" =>
                        "zenith angle of horizon line"])
            close(hlm)

            svf_p_all_w, svf_h_all_w, ft_newdim_w = getnewdims(svf_p_w,svf_h_w,for_tau_w,p,locs)
            svf_p_all_s, svf_h_all_s, ft_newdim_s = getnewdims(svf_p_s,svf_h_s,for_tau_s,p,locs)
            svf_p_all_t, svf_h_all_t, ft_newdim_t = getnewdims(svf_p_t,svf_h_t,for_tau_t,p,locs)

	    defVar(ds,"svf_planar_winter",Int8.(reverse(reshape(svf_p_all_w,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for winter (leaf-off) canopy conditions",])
	    defVar(ds,"svf_hemi_winter",Int8.(reverse(reshape(svf_h_all_w,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for winter (leaf-off) canopy conditions",])

        defVar(ds,"svf_planar_summer",Int8.(reverse(reshape(svf_p_all_s,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for summer (leaf-on) canopy conditions",])
	    defVar(ds,"svf_hemi_summer",Int8.(reverse(reshape(svf_h_all_s,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for summer (leaf-on) canopy conditions",])

        defVar(ds,"svf_planar_terrain",Int8.(reverse(reshape(svf_p_all_t,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for terrain only",])
	    defVar(ds,"svf_hemi_terrain",Int8.(reverse(reshape(svf_h_all_t,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for terrain only",])

        defVar(ds,"Forest_Transmissivity_winter",Int8.(ft_newdim_w),("datetime","locY","locX",),fillvalue = Int8(-99),
                    deflatelevel=5,attrib=["comments" => 
                    "calcualated for winter (leaf-off) canopy conditions",])
        defVar(ds,"Forest_Transmissivity_summer",Int8.(ft_newdim_s),("datetime","locY","locX",),fillvalue = Int8(-99),
                    deflatelevel=5,attrib=["comments" => 
                    "calcualated for summer (leaf-on) canopy conditions",])
        defVar(ds,"Forest_Transmissivity_terrain",Int8.(ft_newdim_t),("datetime","locY","locX",),fillvalue = Int8(-99),
                    deflatelevel=5,attrib=["comments" => 
                    "calcualated for terrain only",])

        close(ds)     




end


function getnewdims(svf_p,svf_h,for_tau,p,locs)

    svf_p[ismissing.(svf_p)] .= -99
    svf_p_all = fill(-99,10000,1)
    svf_p_all[locs] .= (svf_p[p])

    svf_h[ismissing.(svf_h)] .= -99
    svf_h_all = fill(-99,10000,1)
    svf_h_all[locs] .= (svf_h[p])

    for_tau[ismissing.(for_tau)] .= -99
    ft_all = fill(-99,8784,10000)
    ft_all[:,locs] .= for_tau[:,p]
    ft_newdim = fill(-99,8784,100,100)
    for x in 1:1:size(ft_all,1)
        ft_newdim[x,:,:] = reverse(reshape(ft_all[x,:],100,100),dims=1)
    end
    
    return svf_p_all, svf_h_all, ft_newdim

end