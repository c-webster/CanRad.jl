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
    lastc::Vector{Float64},phenology::String)

    # Random branch angle for each tree
    ang = rand(Uniform(60, 100), length(lastc))

    # Load point cloud
    if phenology != "complete"
        pc = getfield(PointCloud(fname; attributes = (classification, :tree_num => intensity)), :data)
        pc.fortype .= Int32.(0)
    else
        pc = getfield(PointCloud(fname; attributes = (classification, :tree_num => intensity, :fortype => pt_src_id)), :data)
        # !!! needs validating but option currently disabled in main anyway
        # pc.fortype should be Int32
    end
    # in this case the intensity field might have been replaced with a forest type classification for
    # each tree crown (i.e. deciduous or evergreen) - allowing distinction for leaf on and off conditions
    # in deciduous forests.  

    pc.tree_num = Int16.(round.(pc.tree_num.*65535)) # scale back intensity to the raw value
        
    mask = (pc.x .>= limits[1]) .& (pc.x .<= limits[2]) .&
        (pc.y .>= limits[3]) .& (pc.y .<= limits[4]) .&
        (pc.classification .!= 2) .& (in.(pc.tree_num, Ref(Set(lastc))))

    keepat!(pc, mask)

    # Add new columns directly to pc
    n = nrow(pc)
    pc.tx = fill(NaN, n)
    pc.ty = fill(NaN, n)
    pc.hd = fill(NaN, n)
    pc.ang = fill(NaN, n)

    # Assign trunk coordinates and branch angles
    for tx in eachindex(lastc)
        t_dx = (lastc[tx] .== pc.tree_num)
        pc.tx[t_dx] .= dbh_x[tx]
        pc.ty[t_dx] .= dbh_y[tx]
        pc.hd[t_dx] .= sqrt.((pc.x[t_dx] .- dbh_x[tx]).^2 .+ (pc.y[t_dx] .- dbh_y[tx]).^2)
        pc.ang[t_dx] .= ang[tx]
    end

    return pc

end

function createfiles(outdir::String,outstr::String,pts::Matrix{Float64},st::SETTINGS,fileio::FILEIO,loc_time=nothing)

    envstrings = environments_flag(st)

    outfile  = joinpath(outdir,"Output_"*outstr*".nc")

    ds = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(ds,"Coordinates",size(pts,1))
    defVar(ds,"easting",pts[:,1],("Coordinates",))
    defVar(ds,"northing",pts[:,2],("Coordinates",))

    if st.calc_trans
        defDim(ds,"datetime",length(loc_time))
        @unpack dt_desc = fileio
        defVar(ds,"datetime",loc_time,("datetime",),attrib=["time_zone" => dt_desc])
    end

    defl_val = 5

    @unpack percent_desc, 
        svf_planar_desc, 
        svf_planar_calc,
        svf_hemi_desc,
        svf_hemi_calc = fileio

    for envname in envstrings

        outsuf = "_$(envname)"

        (envname == "terrain") && (env_desc = fileio.terrain_desc)
        (envname == "leafon") && (env_desc = fileio.leafon_desc)
        (envname == "leafoff") && (env_desc = fileio.leafoff_desc)
        (envname == "evergreen") && (env_desc = fileio.evergreen_desc)
 
        defVar(ds,"svf_planar$(outsuf)",Int32,("Coordinates",),deflatelevel=defl_val,fillvalue = Int8(-1),
            attrib=[
                "long_name"   => "sky-view factor planar",
                "scaling"     => percent_desc,
                "description" => svf_planar_desc,
                "calculation" => svf_planar_calc,
                "environment" => env_desc,])

        defVar(ds,"svf_hemi$(outsuf)",Int32,("Coordinates",),deflatelevel=defl_val,fillvalue = Int8(-1),
            attrib=[
                "long_name"   => "sky-view factor hemispherical",
                "scaling"     => percent_desc,
                "description" => svf_hemi_desc,
                "calculation" => svf_hemi_calc,
                "environment" => env_desc,])

    
        if st.calc_trans

            defVar(ds,"tvt$(outsuf)",Int32,("datetime","Coordinates",),fillvalue = Int8(-1),
                        deflatelevel=defl_val,                            
                        attrib=[
                            "long_name"   => "time-varying direct-beam transmissivity",
                            "scaling"     => percent_desc,
                            "environment" => env_desc,])

            if st.calc_swr > 0

                @unpack swr_tot_desc, 
                    swr_tot_calc,
                    swr_dir_calc,
                    swr_calc = fileio

                defVar(ds,"swr_total$(outsuf)",Int32,("datetime","Coordinates",),
                            deflatelevel=defl_val,fillvalue = Int32(-1),
                            attrib=[
                                "long_name"   => "total incoming shortwave radiation",
                                "units"       => "Watts per metre squared",
                                "description" => swr_tot_desc,
                                "calculation" => "$swr_tot_calc and $swr_calc",
                                "environment" => env_desc,])
                defVar(ds,"swr_direct$(outsuf)",Int32,("datetime","Coordinates",),
                            deflatelevel=defl_val,fillvalue = Int32(-1),
                            attrib=[
                                "long_name"   => "direct incoming shortwave radiation",
                                "units"       => "Watts per metre squared",
                                "calculation" => "$swr_dir_calc and $swr_calc",
                                "environment" => env_desc,])

            end
        end

    end

    return ds

end

function environments_flag(st::SETTINGS)

    envstrings = String[]
    st.calc_terrain && push!(envstrings, "terrain")
    (st.phenology == "leafon" || st.phenology == "both") && push!(envstrings, "leafon")
    (st.phenology == "leafoff" || st.phenology == "both") && push!(envstrings, "leafoff")
    st.forest_type == "evergreen" && push!(envstrings, "evergreen")

    return envstrings

end


function create_exmat(outdir::String,outstr::String,pts::Matrix{Float64},g_img::Matrix{Int64},st::SETTINGS)

    envstrings = environments_flag(st)

    outfile  = joinpath(outdir,"SHIs_"*outstr*".nc")

    images = NCDataset(outfile,"c",format=:netcdf4_classic)

    defDim(images,"img_x",size(g_img,1))
    defDim(images,"img_y",size(g_img,2))
    defDim(images,"Coordinates",size(pts,1))

    defVar(images,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
    defVar(images,"northing",pts[:,2],("Coordinates",),deflatelevel=1)

    for envname in envstrings

        outsuf = "_$(envname)"

        defVar(images,"SHI$(outsuf)",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)

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


function write_metadata(outdir::String,outstr::String,st::SETTINGS,fp::FILEPATHS)

    settings_file = joinpath(outdir, "Settings_$(outstr).txt")
    open(settings_file, "w") do io
        println(io, "=== SETTINGS (st) ===")
        for field in fieldnames(typeof(st))
            println(io, "$field: $(getfield(st, field))")
        end
        println(io, "\n=== DATA (fp) ===")
        for field in fieldnames(typeof(fp))
            println(io, "$field: $(getfield(fp, field))")
        end
        println(io, "\n=== VERSION ===")
        println(io,string("CanRad"*" v"*string(get_pkg_version("CanRad"))))
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

    files = filter(x -> startswith(x, "SHIs") && endswith(x, ".nc"), readdir(datdir))
    images = NCDataset(joinpath(datdir, files[1]), "r")

    coords_x = images["easting"][:]
    coords_y = images["northing"][:]

    varnames = filter(name -> startswith(String(name), "SHI_"), keys(images))

    odir = joinpath(datdir,files[1][1:end-3])

    !ispath(odir) && mkpath(odir)

    fstr = "%07.$(2)f"
    for varname in varnames

        varodir = size(varnames)[1] == 1 ? odir : joinpath(odir, varname[5:end])
        outsuf  = varname[4:end]
        for ix in eachindex(coords_x)
            outf = joinpath(varodir,"SHI_$(cfmt.(fstr,coords_x[ix]))_$(cfmt.(fstr,coords_y[ix]))$(outsuf).png")
            save(outf,colorview(Gray,float.(images["SHI$(outsuf)"][:,:,ix])))
        end

    end

    close(images)

end

function getterrainmask(canrad::CANRAD,terf::String,pts_x::Vector{Float64},pts_y::Vector{Float64})

    terrain_mask = Array{Int8}(undef,(size(canrad.mat2ev,1),size(canrad.mat2ev,2),size(pts_x,1)))

    tm_ds = NCDataset(terf)
    easting = tm_ds["easting"][:]
    northing = tm_ds["northing"][:]

    SHIs = tm_ds["SHI_terrain"][:,:,:]

    if (size(pts_x,1) .== size(easting,1)) && (sum(abs.(easting .- pts_x)) .== 0.0) && (sum(abs.(northing .- pts_y)) .== 0.0)
        terrain_mask .= tm_ds["SHI_terrain"][:,:,:]
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

    SHIs = tm_ds["SHI_terrain"][:,:,:]

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

    error("load_hlm_oshd function deprecated for now.
        use output from save_hlm instead")

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
                            ptsy::Vector{Float64},exdir::String,par_in::Dict{String, Any},
                            ptdx::Vector{Float64},pt_spacing::Int64,fulltilesize::Int64,sub_tilesize::Int64,hlmexdir::String="",
                            shiexdir::String="")

    if (par_in["special_implementation"] !== "swissrad") && (par_in["special_implementation"] !== "oshd") && (par_in["special_implementation"] !== "oshd-alps")
        @warn "collate2tilefile() only written and tested for swissrad and oshd implementation of CanRad"
    end

    if haskey(par_in,"save_horizon") && par_in["save_horizon"]
        save_hlm = true
        !ispath(hlmexdir) && mkpath(hlmexdir)
    else
        save_hlm = false
    end

    if haskey(par_in,"save_images") && par_in["save_images"]
        move_shis = true
        !ispath(joinpath(shiexdir,input)) && mkpath(joinpath(shiexdir,input))
    else
        move_shis = false
    end

    numptsgrid = div(fulltilesize,pt_spacing)
    numptsvec = numptsgrid^2

    # combine tile output to one file
    subtiles = filter(!startswith("."), readdir(outdir)) # avoids pesky .DS_Store files when using MacOS

    xlim = collect(limits[1]:pt_spacing:limits[2]-pt_spacing)
    ylim = collect(limits[3]:pt_spacing:limits[4]-pt_spacing)

    # create the Radiation output dataset
    loc_time = NCDataset(joinpath(outdir,subtiles[1],"Output_"*subtiles[1]*".nc"))["datetime"][:]

    ds = NCDataset(joinpath(exdir,"Output_"*input*".nc"),"c",format=:netcdf4_classic)
    defDim(ds,"locX",size(xlim,1)); defDim(ds,"locY",size(ylim,1))
    defVar(ds,"Easting",Float64.(unique(ptsx)),("locX",)); defVar(ds,"Northing",Float64.(reverse(unique(ptsy))),("locY",))
    defDim(ds,"DateTime",size(loc_time,1))
    time_zone = par_in["time_zone"]
    if  time_zone >= 0
        dt_comment = "time zone UTC+"*string(time_zone)
    elseif time_zone < 0
        dt_comment = "time zone UTC-"*string(time_zone)*"; timestamp is for beginning of averaged period"
    end
    defVar(ds,"datetime",loc_time,("datetime",),attrib=["comments" => dt_comment])

    numptstime = size(loc_time)[1]

    if save_hlm
        # create the horizon line output dataset
        hlm = NCDataset(joinpath(hlmexdir,"HLM_"*input*".nc"),"c",format=:netcdf4_classic)
        phi_bins = NCDataset(joinpath(outdir,subtiles[1],"HLM_"*subtiles[1]*".nc"))["phi"][:]
        defDim(hlm,"locX",size(xlim,1)); defDim(hlm,"locY",size(ylim,1))
        defVar(hlm,"Easting",Float64.(unique(ptsx)),("locX",)); defVar(hlm,"Northing",Float64.(reverse(unique(ptsy))),("locY",))
        defDim(hlm,"phi",size(phi_bins,1))
        defVar(hlm,"Horizon_Line_phi",Int16.(phi_bins),("phi",),deflatelevel=1,attrib=["units" =>
            "degrees clockwise from north"])
    end

    # create collate flags for terrain, leafon and leafoff
    if (haskey(par_in,"calc_terrain") && par_in["calc_terrain"]) || (haskey(par_in,"SHI_terrain") && par_in["SHI_terrain"])
       terrain = true
       svf_p_t = Float64[]
       svf_h_t = Float64[]
       for_tau_t = Float64[]
    else
        terrain = false
    end

    if (haskey(par_in,"phenology") && ((par_in["phenology"] == "both") || (par_in["phenology"] == "leafon"))) ||
         (haskey(par_in,"SHI_leafon") && par_in["SHI_leafon"])
       leafon = true
       svf_p_s = Float64[]
       svf_h_s = Float64[]
       for_tau_s = Float64[]
    else
        leafon = false
    end

    if (haskey(par_in,"phenology") && ((par_in["phenology"] == "both") || (par_in["phenology"] == "leafoff"))) ||
        (haskey(par_in,"SHI_leafoff") && par_in["SHI_leafoff"])
       leafoff = true
       svf_p_leafoff = Float64[]
       svf_h_leafoff = Float64[]
       for_tau_leafoff = Float64[]
    else
        leafoff = false
    end

    pts_all = hcat(ptsx[ptdx.>0],ptsy[ptdx.>0],ptdx[ptdx.>0])

    limx = Int.((floor(limits[1]/sub_tilesize))*sub_tilesize):sub_tilesize:(Int.(floor((limits[2]/sub_tilesize))*sub_tilesize))
    limy = (Int.((floor(limits[3]/sub_tilesize))*sub_tilesize):sub_tilesize:(Int.(floor(((limits[4])/sub_tilesize))*sub_tilesize)))

    # initilise other variables
    east = Float64[]
    north = Float64[]
    save_hlm && (tht = Float64[])

    for tx in subtiles

        tds = NCDataset(joinpath(outdir,tx,"Output_"*tx*".nc"))
        save_hlm && (pds = NCDataset(joinpath(outdir,tx,"HLM_"*tx*".nc")))

        xdx = findall(parse(Int,(split(tx,"_")[1]))[1] .== limx)[1]
        ydx = findall(parse(Int,(split(tx,"_")[2]))[1] .== limy)[1]

        idx = (limx[xdx] .<= pts_all[:,1] .< limx[xdx+1]) .&
                (limy[ydx] .<= pts_all[:,2] .< limy[ydx+1])

        pts = pts_all[idx,:]

        append!(east, tds["easting"][:])
        append!(north, tds["northing"][:])

        if terrain
            (append!(svf_p_t,tds["svf_planar_terrain"][:]))
            (append!(svf_h_t,tds["svf_hemi_terrain"][:]))
            if isempty(for_tau_t)
                for_tau_t = tds["trans_terrain"][:,:]
            else
                for_tau_t = hcat(for_tau_t, tds["trans_terrain"][:,:])
            end
        end

        if !(sum(pts[:,3]) == (size(pts,1) * 2)) # checks whether this sub-tile was calculated with ter2rad

            leafoff && (append!(svf_p_leafoff,tds["svf_planar_leafoff"][:]))
            leafon && (append!(svf_p_leafon,tds["svf_planar_leafon"][:]))

            leafoff && (append!(svf_h_leafoff,tds["svf_hemi_leafoff"][:]))
            leafon && (append!(svf_h_leafon,tds["svf_hemi_leafon"][:]))

            if leafoff && isempty(for_tau_leafoff)
                for_tau_leafoff = tds["for_trans_leafoff"][:,:]
            elseif leafoff
                for_tau_leafoff = hcat(for_tau_leafoff, tds["for_trans_leafoff"][:,:])
            end

            if leafon && isempty(for_tau_leafon)
                for_tau_leafon = tds["for_trans_leafon"][:,:]
            elseif leafon
                for_tau_leafon = hcat(for_tau_leafon, tds["for_trans_leafon"][:,:])
            end

        else # if subtile is terrain only (calculated with ter2rad), fill forest variables with nodata

            dummvf = fill(-1,size(tds["svf_planar_t"][:]))
            leafoff && (append!(svf_p_leafoff,dummvf))
            leafon && (append!(svf_p_leafon,dummvf))

            leafoff && (append!(svf_h_leafoff,dummvf))
            leafon && (append!(svf_h_leafon,dummvf))

            dummytau = fill(-1,size(tds["trans_t"][:,:]))
            if leafoff && isempty(for_tau_leafoff)
                for_tau_leafoff = tds["trans_leafoff"][:,:]
            elseif leafoff
                for_tau_leafoff = hcat(for_tau_leafoff, dummytau)
            end

            if leafon && isempty(for_tau_leafon)
                for_tau_leafon = tds["trans_leafon"][:,:]
            elseif leafon
                for_tau_leafon = hcat(for_tau_leafon, dummytau)
            end

        end

        if save_hlm && isempty(tht)
            tht = pds["tht"][:,:]
        elseif save_hlm
            tht = hcat(tht,pds["tht"][:,:])
        end

        close(tds)
        save_hlm && close(pds)

        if move_shis
            cp(joinpath(outdir,tx,"SHIs_"*tx*".nc"),joinpath(shiexdir,input,"SHIs_"*tx*".nc"),force=true)
        end

    end

    B = hcat(east,north)
    p = sortperm(view.(Ref(B), 1:size(B,1), :));

    locs = (ptdx .> 0)

    if save_hlm
        tht_all = fill(-1,360,numptsvec)
        tht_all[:,locs] .= tht[:,p]

        tht_newdim = fill(-1,360,numptsgrid,numptsgrid)
        for x in 1:1:size(tht_all,1)
            tht_newdim[x,:,:] = reverse(reshape(tht_all[x,:],numptsgrid,numptsgrid),dims=1)
        end

        defVar(hlm,"Horizon_Line_tht",Int16.(tht_newdim),("phi","locY","locX",),deflatelevel=1,attrib=["units" =>
                    "zenith angle of horizon line"])
        close(hlm)
    end

    if leafoff 
        svf_p_all_leafoff, svf_h_all_leafoff, ft_newdim_leafoff = getnewdims(svf_p_leafoff,svf_h_leafoff,for_tau_leafoff,p,locs,numptsvec,numptsgrid,numptstime); end

    if leafon
        svf_p_all_leafon, svf_h_all_leafon, ft_newdim_leafon = getnewdims(svf_p_leafon,svf_h_leafon,for_tau_leafon,p,locs,numptsvec,numptsgrid,numptstime); end

    if terrain
        svf_p_all_t, svf_h_all_t, ft_newdim_t = getnewdims(svf_p_t,svf_h_t,for_tau_t,p,locs,numptsvec,numptsgrid,numptstime); end

    if leafoff
        defVar(ds,"svf_planar_leafoff",Int8.(reverse(reshape(svf_p_all_leafoff,numptsgrid,numptsgrid),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-1),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calculated for leafoff canopy conditions",])
        defVar(ds,"svf_hemi_leafoff",Int8.(reverse(reshape(svf_h_all_leafoff,numptsgrid,numptsgrid),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-1),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calculated for leafoff canopy conditions",])
        defVar(ds,"transmissivity_leafoff",Int8.(ft_newdim_leafoff),("datetime","locY","locX",),fillvalue = Int8(-1),
                    deflatelevel=5,attrib=["comments" => 
                    "values are represented as percentage;
                    calculated for leafoff canopy conditions",])
    end

    if leafon
        defVar(ds,"svf_planar_leafon",Int8.(reverse(reshape(svf_p_all_leafon,numptsgrid,numptsgrid),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-1),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calculated for leafon canopy conditions",])
        defVar(ds,"svf_hemi_leafon",Int8.(reverse(reshape(svf_h_all_leafon,numptsgrid,numptsgrid),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-1),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calculated for leafon canopy conditions",])
        defVar(ds,"transmissivity_leafon",Int8.(ft_newdim_leafon),("datetime","locY","locX",),fillvalue = Int8(-1),
                    deflatelevel=5,attrib=["comments" => 
                    "values are represented as percentage;
                    calculated for leafon canopy conditions",])
    end

    if terrain
        defVar(ds,"svf_planar_terrain",Int8.(reverse(reshape(svf_p_all_t,numptsgrid,numptsgrid),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-1),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of a horizontal flat uplooking surface;
                    zenith rings weighted by surface area projected onto a horizontal flat surface;
                    calculated for terrain only (no canopy)",])
        defVar(ds,"svf_hemi_terrain",Int8.(reverse(reshape(svf_h_all_t,numptsgrid,numptsgrid),dims=1)),("locY","locX",),deflatelevel=5,fillvalue = Int8(-1),
                    attrib=["comments" =>
                    "values are represented as percentage;
                    perspective of hemipherically shaped surface or plant;
                    zenith rings weighted by surface area on the hemisphere;
                    calculated for terrain only (no canopy)",])
        defVar(ds,"transmissivity_terrain",Int8.(ft_newdim_t),("datetime","locY","locX",),fillvalue = Int8(-1),
                    deflatelevel=5,attrib=["comments" => 
                    "values are represented as percentage;
                    calculated for terrain only (no canopy)",])
    end

    close(ds)     

end


function getnewdims(svf_p,svf_h,for_tau,p,locs,numptsvec,numptsgrid,numptstime)

    svf_p[ismissing.(svf_p)] .= -1
    svf_p_all = fill(-1,numptsvec,1)
    svf_p_all[locs] .= (svf_p[p])

    svf_h[ismissing.(svf_h)] .= -1
    svf_h_all = fill(-1,numptsvec,1)
    svf_h_all[locs] .= (svf_h[p])

    for_tau[ismissing.(for_tau)] .= -1
    ft_all = fill(-1,numptstime,numptsvec)
    ft_all[:,locs] .= for_tau[:,p]
    ft_newdim = fill(-1,numptstime,numptsgrid,numptsgrid)
    for x in 1:1:size(ft_all,1)
        ft_newdim[x,:,:] = reverse(reshape(ft_all[x,:],numptsgrid,numptsgrid),dims=1)
    end
    
    return svf_p_all, svf_h_all, ft_newdim

end