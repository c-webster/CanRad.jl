function loaddbh(fname::String,limits::Matrix{Float64},peri=0::Int64)
    dat, _ = readdlm(fname,'\t',header=true)

    dat_x, dat_y, dat_z, kpdx = clipdat(dat[:,1],dat[:,2],dat[:,3],limits,peri)

    dat_r = (dat[kpdx,4] ./2) ./100

    if size(dat,2) == 5
        dat_tc = dat[kpdx,5]
    else
        dat_tc = fill(NaN,size(dat_x))
    end

    return dat_x, dat_y, dat_z, dat_r, dat_tc
end

function loadltc_laz(fname::String,limits::Matrix{Float64},
    dbh_x::Vector{Float64},dbh_y::Vector{Float64},
    lastc::Vector{Float64},season::String)

    # calculate random branch angle for each tree
    ang = rand(Uniform(60,100),size(lastc,1)) # Norway Spruce

    header, ltcdat = LazIO.load(fname)

    dat = DataFrame(ltcdat)

    dx = ((limits[1] .<= (dat.x .* header.x_scale .+ header.x_offset) .<= limits[2]) .&
        (limits[3] .<= (dat.y .* header.y_scale .+ header.y_offset) .<= limits[4])) .&
        (dat.raw_classification .!= 2) #.& (ignore .!= nothing)

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


function createfiles(outdir::String,outstr::String,pts::Matrix{Float64},calc_trans::Bool,calc_swr::Int64,
            append_file::Bool,loc_time::Vector{DateTime}=nothing,time_zone::Int64=nothing,season::String="none")

    outfile  = outdir*"/Output_"*outstr*".nc"

    # append_file option removed September 2021 (eventually should be included)

	if append_file

		ds = NCDataset(outfile,"a")

	else

	    ds = NCDataset(outfile,"c",format=:netcdf4_classic)

	    defDim(ds,"Coordinates",size(pts,1))
	    defVar(ds,"easting",pts[:,1],("Coordinates",))
	    defVar(ds,"northing",pts[:,2],("Coordinates",))

        if season != "complete"
            defVar(ds,"Vf_planar",Int32,("Coordinates",),deflatelevel=5,
                        attrib=["scale_factor"=>0.01, "comments" =>
                        "perspective of a horizontal flat uplooking surface;
                        zenith rings weighted by surface area projected onto a horizontal flat surface",])
            defVar(ds,"Vf_hemi",Int32,("Coordinates",),deflatelevel=5,
                        attrib=["scale_factor"=>0.01, "comments" =>
                        "perspective of hemipherically shaped surface or plant;
                        zenith rings weighted by surface area on the hemisphere",])
        end

	    if calc_trans
	        defDim(ds,"datetime",length(loc_time))

			if time_zone >= 0
				dt_comment = "time zone UTC+"*string(time_zone)
			elseif time_zone < 0
				dt_comment = "time zone UTC-"*string(time_zone)
			end
			defVar(ds,"datetime",loc_time,("datetime",),attrib=["comments" => dt_comment])

            if season == "complete"
                defVar(ds,"Vf_planar",Int32,("datetime","Coordinates",),deflatelevel=5,
                            attrib=["scale_factor"=>0.01, "comments" =>
                            "perspective of a horizontal flat uplooking surface;
                            zenith rings weighted by surface area projected onto a horizontal flat surface",])
                defVar(ds,"Vf_hemi",Int32,("datetime","Coordinates",),deflatelevel=5,
                            attrib=["scale_factor"=>0.01, "comments" =>
                            "perspective of hemipherically shaped surface or plant;
                            zenith rings weighted by surface area on the hemisphere",])
            end

	        defVar(ds,"Forest_Transmissivity",Int32,("datetime","Coordinates",),
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


function create_exmat(outdir::String,outstr::String,pts::Matrix{Float64},g_img::Matrix{Int64},
                append_file::Bool,season::String="none")

    outfile  = outdir*"/SHIs_"*outstr*".nc"

    if append_file

		images = NCDataset(outfile,"a")

	else

        images = NCDataset(outfile,"c",format=:netcdf4_classic)

        defDim(images,"img_x",size(g_img,1))
        defDim(images,"img_y",size(g_img,2))
        defDim(images,"Coordinates",size(pts,1))

        defVar(images,"easting",pts[:,1],("Coordinates",),deflatelevel=1)
        defVar(images,"northing",pts[:,2],("Coordinates",),deflatelevel=1)

        if season != "complete"
            defVar(images,"SHI",Int8,("img_y","img_x","Coordinates",),deflatelevel=1)
        else
            defDim(images,"Season",2)
            defVar(images,"SHI",Int8,("img_y","img_x","Coordinates","Season",),deflatelevel=1,
                    attrib=["seasonal dimension"=>"1. winter, 2. summer",])
        end

	end

    return images

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
        outstr = String(split(exdir,"/")[end-1])
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
		percentdone = percentdone = Int(floor((crxstart / numpts) * 100))
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

    if !ispath(odir)
        mkpath(odir)
    end

    fstr = "%07.$(2)f"

    for ix in eachindex(coords_x)

        outf = joinpath(odir,"SHI_"*sprintf1.(fstr,coords_x[ix])*"_"*sprintf1.(fstr,coords_y[ix])*".png")
        save(outf,colorview(Gray,float.(images["SHI"][:,:,ix])))

    end

    close(images)
end