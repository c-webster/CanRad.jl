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
            append_file::Bool,loc_time=nothing,time_zone=nothing)

    outfile  = outdir*"/Output_"*outstr*".nc"

    # append_file option removed September 2021 (eventually should be included)

	if append_file

		ds = NCDataset(outfile,"a")

	else

	    ds = NCDataset(outfile,"c",format=:netcdf4_classic)

	    defDim(ds,"Coordinates",size(pts,1))
	    defVar(ds,"easting",Float32.(pts[:,1]),("Coordinates",))
	    defVar(ds,"northing",Float32.(pts[:,2]),("Coordinates",))

	    defVar(ds,"Vf_planar_winter",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for winter canopy conditions",])
	    defVar(ds,"Vf_hemi_winter",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for winter canopy conditions",])

        defVar(ds,"Vf_planar_summer",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for summer canopy conditions",])
	    defVar(ds,"Vf_hemi_summer",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for summer canopy conditions",])

        defVar(ds,"Vf_planar_terrain",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for terrain only",])
	    defVar(ds,"Vf_hemi_terrain",Int8,("Coordinates",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for terrain only",])

	    if calc_trans
	        defDim(ds,"datetime",length(loc_time))

			if time_zone >= 0
				dt_comment = "time zone UTC+"*string(time_zone)
			elseif time_zone < 0
				dt_comment = "time zone UTC-"*string(time_zone)*"; timestamp is for beginning of averaged period"
			end
			defVar(ds,"datetime",loc_time,("datetime",),attrib=["comments" => dt_comment])

	        defVar(ds,"Forest_Transmissivity_winter",Int8,("datetime","Coordinates",),fillvalue = Int8(-99),
                        deflatelevel=5,attrib=["comments" => 
                        "calcualated for winter canopy conditions",])
            defVar(ds,"Forest_Transmissivity_summer",Int8,("datetime","Coordinates",),fillvalue = Int8(-99),
                        deflatelevel=5,attrib=["comments" => 
                        "calcualated for summer canopy conditions",])
            defVar(ds,"Forest_Transmissivity_terrain",Int8,("datetime","Coordinates",),fillvalue = Int8(-99),
                        deflatelevel=5,attrib=["comments" => 
                        "calcualated for terrain only",])

	        if calc_swr > 0
	            defVar(ds,"SWR_total_winter",Int16,("datetime","Coordinates",),
                            deflatelevel=5,fillvalue = Int16(-999),
                            attrib=["units"=>"Watts per metre squared", "comments" =>
                            "total incoming shortwave radiation (diffuse + direct);
                            calcualated for winter canopy conditions",])
	            defVar(ds,"SWR_direct_winter",Int16,("datetime","Coordinates",),
	                        deflatelevel=5,fillvalue = Int16(-999),
	                        attrib=["units"=>"Watts per metre squared;
                            calcualated for winter canopy conditions",])

                defVar(ds,"SWR_total_summer",Int16,("datetime","Coordinates",),
                            deflatelevel=5,fillvalue = Int16(-999),
                            attrib=["units"=>"Watts per metre squared", "comments" =>
                            "total incoming shortwave radiation (diffuse + direct);
                            calcualated for summer canopy conditions",])
	            defVar(ds,"SWR_direct_summer",Int16,("datetime","Coordinates",),
	                        deflatelevel=5,fillvalue = Int16(-999),
	                        attrib=["units"=>"Watts per metre squared;
                            calcualated for summer canopy conditions",])

                defVar(ds,"SWR_total_terrain",Int16,("datetime","Coordinates",),
                            deflatelevel=5,fillvalue = Int16(-999),
                            attrib=["units"=>"Watts per metre squared", "comments" =>
                            "total incoming shortwave radiation (diffuse + direct);
                            calcualated for terrain only",])
	            defVar(ds,"SWR_direct_terrain",Int16,("datetime","Coordinates",),
	                        deflatelevel=5,fillvalue = Int16(-999),
	                        attrib=["units"=>"Watts per metre squared", "comments" =>"
                            calcualated for terrain only",])

			end
	    end

	end

	return ds

end




function create_exmat(outdir::String,outstr::String,pts::Matrix{Float64},g_img::Matrix{Int64},append_file::Bool)

    outfile  = outdir*"/SHIs_"*outstr*".nc"

    if append_file

		images = NCDataset(outfile,"a")

	else

        images = NCDataset(outfile,"c",format=:netcdf4_classic)

        defDim(images,"img_x",size(g_img,1))
        defDim(images,"img_y",size(g_img,2))
        defDim(images,"Coordinates",size(pts,1))

        defVar(images,"SHI_winter",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
        defVar(images,"SHI_summer",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
        defVar(images,"SHI_terrain",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)

        defVar(images,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
        defVar(images,"northing",pts[:,2],("Coordinates",),deflatelevel=1)

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

    defVar(hlm,"easting",Float32.(pts[:,1]),("Coordinates",),deflatelevel=1)
    defVar(hlm,"northing",Float32.(pts[:,2]),("Coordinates",),deflatelevel=1)
    defVar(hlm,"phi",Int16.((round.(phi_bins))),("phi",),deflatelevel=1)
    # defVar(hlm,"phi_oshd",Int16.((round.(phi_bins))[p_srt]),("phi",),deflatelevel=1)

    defVar(hlm,"tht",Int8,("phi","Coordinates",),fillvalue = Int8(-99),
            deflatelevel=5,attrib=["comments"=>"zenith angle of horizon line"],)
    # defVar(hlm,"tht_oshd",Int8,("phi","Coordinates",),fillvalue = Int8(-99),
    #         deflatelevel=5,attrib=["scale_factor"=>0.01,"comments"=>"zenith angle of horizon line"],)

    defVar(hlm,"phi_sortdx",Int16.(p_srt),("phi",),deflatelevel=1,
            attrib=["comments" => "sort index for phi/tht variables for phi values to be sorted 0-360degrees",])

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

    odir = joinpath(datdir,files[fx][1][1:end-3])
    odir_w = joinpath(odir,"winter")
    odir_s = joinpath(odir,"summer")
    odir_t = joinpath(odir,"terrain")

    if !ispath(odir); mkpath(odir); end
    if !ispath(odir_w); mkpath(odir_w); end
    if !ispath(odir_s); mkpath(odir_s); end
    if !ispath(odir_t); mkpath(odir_t); end

    fstr = "%07.$(2)f"

    for ix in eachindex(coords_x)

        outf_w = joinpath(odir_w,"SHI_"*sprintf1.(fstr,coords_x[ix])*"_"*sprintf1.(fstr,coords_y[ix])*"_winter.png")
        save(outf_w,colorview(Gray,float.(images["SHI_winter"][:,:,ix])))

        outf_s = joinpath(odir_s,"SHI_"*sprintf1.(fstr,coords_x[ix])*"_"*sprintf1.(fstr,coords_y[ix])*"_summer.png")
        save(outf_s,colorview(Gray,float.(images["SHI_summer"][:,:,ix])))

        outf_t = joinpath(odir_t,"SHI_"*sprintf1.(fstr,coords_x[ix])*"_"*sprintf1.(fstr,coords_y[ix])*"_terrain.png")
        save(outf_t,colorview(Gray,float.(images["SHI_terrain"][:,:,ix])))

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


function collate2tilefile(outdir::String,limits::Matrix{Int64},input::String,ptsx::Vector{Float64},
                            ptsy::Vector{Float64},exdir::String,hlmexdir::String,par_in::Dict{String, Any},
                            ptdx::Vector{Float64})


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

                    global Vf_p_w = tds["Vf_planar_winter"][:]
                    global Vf_p_s = tds["Vf_planar_summer"][:]
                    global Vf_p_t = tds["Vf_planar_terrain"][:]

                    global Vf_h_w = tds["Vf_hemi_winter"][:]
                    global Vf_h_s = tds["Vf_hemi_summer"][:]
                    global Vf_h_t = tds["Vf_hemi_terrain"][:]

                    global for_tau_w = tds["Forest_Transmissivity_winter"][:]
                    global for_tau_s = tds["Forest_Transmissivity_summer"][:]
                    global for_tau_t = tds["Forest_Transmissivity_terrain"][:]

                    global tht   = pds["tht"][:]

                else

                    append!(east,tds["easting"][:])
                    append!(north,tds["northing"][:])

                    append!(Vf_p_w,tds["Vf_planar_winter"][:])
                    append!(Vf_p_s,tds["Vf_planar_summer"][:])
                    append!(Vf_p_t,tds["Vf_planar_terrain"][:])

                    append!(Vf_h_w,tds["Vf_hemi_winter"][:])
                    append!(Vf_h_s,tds["Vf_hemi_summer"][:])
                    append!(Vf_h_t,tds["Vf_hemi_terrain"][:])

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

            Vf_p_all_w, Vf_h_all_w, ft_newdim_w = getnewdims(Vf_p_w,Vf_h_w,for_tau_w,p,locs)
            Vf_p_all_s, Vf_h_all_s, ft_newdim_s = getnewdims(Vf_p_s,Vf_h_s,for_tau_s,p,locs)
            Vf_p_all_t, Vf_h_all_t, ft_newdim_t = getnewdims(Vf_p_t,Vf_h_t,for_tau_t,p,locs)

	    defVar(ds,"Vf_planar_winter",Int8.(reverse(reshape(Vf_p_all_w,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for winter (leaf-off) canopy conditions",])
	    defVar(ds,"Vf_hemi_winter",Int8.(reverse(reshape(Vf_h_all_w,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for winter (leaf-off) canopy conditions",])

        defVar(ds,"Vf_planar_summer",Int8.(reverse(reshape(Vf_p_all_s,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for summer (leaf-on) canopy conditions",])
	    defVar(ds,"Vf_hemi_summer",Int8.(reverse(reshape(Vf_h_all_s,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calcualated for summer (leaf-on) canopy conditions",])

        defVar(ds,"Vf_planar_terrain",Int8.(reverse(reshape(Vf_p_all_t,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calcualated for terrain only",])
	    defVar(ds,"Vf_hemi_terrain",Int8.(reverse(reshape(Vf_h_all_t,100,100),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-99),
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


function getnewdims(Vf_p,Vf_h,for_tau,p,locs)

    Vf_p[ismissing.(Vf_p)] .= -99
    Vf_p_all = fill(-99,10000,1)
    Vf_p_all[locs] .= (Vf_p[p])

    Vf_h[ismissing.(Vf_h)] .= -99
    Vf_h_all = fill(-99,10000,1)
    Vf_h_all[locs] .= (Vf_h[p])

    for_tau[ismissing.(for_tau)] .= -99
    ft_all = fill(-99,8784,10000)
    ft_all[:,locs] .= for_tau[:,p]
    ft_newdim = fill(-99,8784,100,100)
    for x in 1:1:size(ft_all,1)
        ft_newdim[x,:,:] = reverse(reshape(ft_all[x,:],100,100),dims=1)
    end
    
    return Vf_p_all, Vf_h_all, ft_newdim

end