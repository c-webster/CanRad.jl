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
    _, _, _, rmdx = clipdat(ltc[:,1],ltc[:,2],ltc[:,3],limits,peri)
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
            append_file::Bool,dimensions::Int64=2,loc_time=nothing,time_zone=nothing)

    outfile  = outdir*"/Output_"*outstr*".nc"

    # append_file option removed September 2021 (eventually should be included)

	if append_file

		ds = NCDataset(outfile,"a")

	else

	    ds = NCDataset(outfile,"c",format=:netcdf4_classic)

	    defDim(ds,"Coordinates",size(pts,1))
	    defVar(ds,"easting",pts[:,1],("Coordinates",))
	    defVar(ds,"northing",pts[:,2],("Coordinates",))

        if dimensions == 3
            defVar(ds,"height",pts[:,3],("Coordinates",))
        end

	    defVar(ds,"Vf_planar",Int32,("Coordinates",),deflatelevel=5,fillvalue = Int32(-9999),
                    attrib=["scale_factor"=>0.01, "comments" =>
                    "perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface",])
	    defVar(ds,"Vf_hemi",Int32,("Coordinates",),deflatelevel=5,fillvalue = Int32(-9999),
                    attrib=["scale_factor"=>0.01, "comments" =>
                    "perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere",])

	    if calc_trans
	        defDim(ds,"datetime",length(loc_time))

			if time_zone >= 0
				dt_comment = "time zone UTC+"*string(time_zone)
			elseif time_zone < 0
				dt_comment = "time zone UTC-"*string(time_zone)
			end
			defVar(ds,"datetime",loc_time,("datetime",),attrib=["comments" => dt_comment])

	        defVar(ds,"Forest_Transmissivity",Int32,("datetime","Coordinates",),fillvalue = Int32(-9999),
                        deflatelevel=5,attrib=["scale_factor"=>0.01,])

	        if calc_swr > 0
	            defVar(ds,"SWR_total",Int32,("datetime","Coordinates",),
                            deflatelevel=5,fillvalue = Int32(-9999),
                            attrib=["units"=>"Watts per metre squared", "comments" =>
                            "total incoming shortwave radiation (diffuse + direct)"])
	            defVar(ds,"SWR_direct",Int32,("datetime","Coordinates",),
	                        deflatelevel=5,fillvalue = Int32(-9999),
	                        attrib=["units"=>"Watts per metre squared",])
			end
	    end

	end

	return ds

end


function create_exmat(outdir::String,outstr::String,pts::Matrix{Float64},
    g_img::Matrix{Int64},append_file::Bool,dimensions::Int64=2)

    outfile  = outdir*"/SHIs_"*outstr*".nc"

    if append_file

		images = NCDataset(outfile,"a")

	else

        images = NCDataset(outfile,"c",format=:netcdf4_classic)

        defDim(images,"img_x",size(g_img,1))
        defDim(images,"img_y",size(g_img,2))
        defDim(images,"Coordinates",size(pts,1))

        defVar(images,"SHI",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
        defVar(images,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
        defVar(images,"northing",pts[:,2],("Coordinates",),deflatelevel=1)

        if dimensions == 3
            defVar(images,"height",pts[:,3],("Coordinates",))
        end

	end

    return images

end

function create_exhlm(outdir::String,outstr::String,pts::Matrix{Float64},ter2rad::TER2RAD)

    outfile  = outdir*"/HLM_"*outstr*".nc"

    @unpack phi_bins, dx1, dx2 = ter2rad

    phi_bins = phi_bins[dx1:dx2-1]

    phi_bins .= ((pi/2) .- phi_bins) .* (180/pi);
    phi_bins[phi_bins.>=360] .= phi_bins[phi_bins.>=360] .- 360;              
    phi_bins[phi_bins.<0]    .= phi_bins[phi_bins.<0] .+ 360;         # phi now is the azimuth angle (phi) in degree, so that 0° is North, 90° is East, 180° is South, and 270° is West

    p_srt = sortperm(phi_bins)

    hlm = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(hlm,"Coordinates",size(pts,1))
    defDim(hlm,"phi",size(phi_bins,1))

    defVar(hlm,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
    defVar(hlm,"northing",pts[:,2],("Coordinates",),deflatelevel=1)
    defVar(hlm,"phi",(round.(phi_bins)),("phi",),deflatelevel=1)
    defVar(hlm,"phi_oshd",(round.(phi_bins))[p_srt],("phi",),deflatelevel=1)

    defVar(hlm,"tht",Int32,("phi","Coordinates",),fillvalue = Int32(-9999),
            deflatelevel=5,attrib=["scale_factor"=>0.01,"comments"=>"zenith angle of horizon line"],)
    defVar(hlm,"tht_oshd",Int32,("phi","Coordinates",),fillvalue = Int32(-9999),
            deflatelevel=5,attrib=["scale_factor"=>0.01,"comments"=>"zenith angle of horizon line"],)

    return hlm, p_srt

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

function organise_outf(taskID::String,exdir::String,batch::Bool,numpts::Int)

	if batch
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3]
        global outdir = exdir*"/"*outstr
    else
        outstr = splitpath(exdir)[end-1]
        global outdir = exdir
    end

    if !ispath(outdir)
        mkpath(outdir)
    end

	# set start point within tile
	fx = readdir(outdir)[findall(startswith.(readdir(outdir),"Processing"))]
	# if tile never started, tile failed on intialisation last time, or tile completed but is being re-run
	if isempty(fx) || parse(Int,split(fx[1])[4]) == 0 || parse(Int,split(fx[1])[4]) == numpts
		crxstart = 1
		append_file = false
		percentdone = 0
	else # restart tile at next point
		crxstart = parse(Int,split(fx[1])[4])+1
		append_file = true
		percentdone = Int(floor((crxstart / numpts) * 100))
	end

	if !isempty(fx)
		rm(joinpath(outdir,fx[1]))
	end

	return outdir, outstr, crxstart, append_file, percentdone

end

function make_SHIs(datdir::String)

    files = readdir(datdir)
    fx = findall(startswith.(files,"SHIs") .& endswith.(files,".nc"))

    images = NCDataset(joinpath(datdir,files[fx][1]),"r")

    coords_x = images["easting"][:]
    coords_y = images["northing"][:]
    coords_z = images["height"][:]

    odir = joinpath(datdir,files[fx][1][1:end-3])

    if !ispath(odir)
        mkpath(odir)
    end

    fstr = "%07.$(2)f"

    for ix in eachindex(coords_x)
        
        if haskey(images,"height")
            outf = joinpath(odir,"SHI_"*sprintf1.(fstr,coords_x[ix])*"_"*sprintf1.(fstr,coords_y[ix])*"_"*sprintf1.(fstr,coords_z[ix])*".png")

        else
            outf = joinpath(odir,"SHI_"*sprintf1.(fstr,coords_x[ix])*"_"*sprintf1.(fstr,coords_y[ix])*".png")
        end
        save(outf,colorview(Gray,float.(images["SHI"][:,:,ix])))

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

    error("you haven't updated the load_hlm_oshd function yet")

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