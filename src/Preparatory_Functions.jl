function compatability_check(dat_in,par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    if !@isdefined(buildings) # include buildings as solid objects (not floating carpets)
        par_in["buildings"] = false
    end

    if !@isdefined(horizon_line) # add in pre-calculated horizon line for OSHD compatability
        par_in["horizon_line"] = 0
    end

    if !@isdefined(pt_corr) # correct for probability of transmissivity
        par_in["pt_corr"] = true
    end

    # renamed terrain options
    if @isdefined(terrain)
        if terrain == 1
            par_in["terrain_tile"] =  false
        elseif terrain == 2
            par_in["terrain_tile"] = true
        end
    end

    if @isdefined(tershad)
        if tershad == 1
            par_in["terrain_highres"] = true
            par_in["terrain_lowres"] = true
        elseif tershad == 2
            par_in["terrain_highres"] = true
            par_in["terrain_lowres"] = false
        elseif tershad == 3
            par_in["terrain_highres"] = false
            par_in["terrain_lowres"] = false
        end
    end

    if !@isdefined(horizon_line)
        par_in["horizon_line"] = false
    end

    # disable tilt option since it is not implemented properly
    if !@isdefined(tilt) || tilt 
        par_in["tilt"] = false
    end

    if !@isdefined(OSHD)
        par_in["OSHD"] = false
    end

	# renamed variables with version 0.7.0
	if !@isdefined(image_height)
		par_in["image_height"] = ch
	end

	if @isdefined(utm_zone) && coor_system == "UTM"
		par_in["coor_system"] = "UTM"*" "*utm_zone
	end

	if !@isdefined(b_space) && !@isdefined(branch_spacing) # spacing between branch points in cm
        par_in["branch_spacing"] = 0.1
	elseif @isdefined(b_space)
		par_in["branch_spacing"] = b_space
    end

    if !@isdefined(season)
        par_in["season"] = "none"
    end

    return dat_in, par_in

end

function extract(d::Dict)
    expr = quote end
    for (k, v) in d
       push!(expr.args, :($(Symbol(k)) = $v))
    end
    # eval(expr)
    return expr
end

function clipdat(pc_x::Vector{Float64},pc_y::Vector{Float64},pc_z::Vector{Float64},limits,peri=0::Int64)
    
    kpdx = (limits[1]-peri .<= pc_x .<= limits[2]+peri) .&& 
                (limits[3]-peri .<= pc_y .<= limits[4]+peri) 

    return pc_x[kpdx], pc_y[kpdx], pc_z[kpdx], kpdx

end

function clipdat(dat::DataFrame,limits,peri=0::Int64)
    
    kpdx = (limits[1]-peri .<= dat.x .<= limits[2]+peri) .&& 
                (limits[3]-peri .<= dat.y .<= limits[4]+peri) 

    return dat[kpdx,:] 

end


function create_tiles(basefolder::String,ptsf::String,settings_fun::Function)

    tsize = readdlm(basefolder*"/TileSize.txt")

    segs = readdir(basefolder*"/Segments")

    try
        mkdir(basefolder*"/Tiles"); catch
    end

    pts_full = (readdlm(ptsf))

    analysis_limits = readdlm(basefolder*"/AnalysisAreaLimits.txt")

    pts_idx = (analysis_limits[1] .<= pts_full[:,1] .< analysis_limits[2]) .&
            (analysis_limits[3] .<= pts_full[:,2] .< analysis_limits[4])

    pts_all = pts_full[pts_idx,:]

    tdx = 0

    ptsfname     = String[]
    inputsegname = String[]

    for segname in segs

        limits = readdlm(basefolder*"/Segments"*"/"*segname*"/"*segname*"_analysisarea.txt")

        limx = (collect(limits[1]:tsize[1]:limits[2]))
        limy = (collect(limits[3]:tsize[1]:limits[4]))

        if limy[end] < limits[4]; push!(limy,limy[end].+tsize[1]); end
        if limx[end] < limits[2]; push!(limx,limx[end].+tsize[1]); end

        for x in eachindex(limx[1:end-1]), y in eachindex(limy[1:end-1])
            # create the tasks
            idx = (limx[x] .<= pts_all[:,1] .< limx[x+1]) .&
                    (limy[y] .<= pts_all[:,2] .< limy[y+1])

            if sum(idx) > 0
                tdx += 1
                # push!(ptsfname,sprintf1.("%03.$(0)f",tdx)*".txt")
                push!(ptsfname,sprintf1.("%03.$(0)f",tdx)*"_"*sprintf1.("%03.$(0)f", limx[x])*"_"*sprintf1.("%03.$(0)f", limy[y])*".txt")
                f = open(basefolder*"/Tiles"*"/"*ptsfname[end],"w")
                    writedlm(f,sprintf1.("%7.$(2)f",pts_all[idx,:]))
                close(f)
                push!(inputsegname,segname)
            end
        end

    end

    return ptsfname, inputsegname, basefolder*"/Output/", collect(1:1:size(ptsfname,1))

end

function check_output(exdir::String,pts::Matrix{Float64},batch::Bool,taskID::String)

    if batch
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3]
        outdir = exdir*"/"*outstr
    else
        outstr = String(split(exdir,"/")[end])
        outdir = exdir
    end

    if ispath(outdir)
        try
            crxstart = parse(Int,(split(reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])))[4])
            if crxstart == size(pts,1)
                out = true
            else
                out = false
            end
        catch
            out = false
        end
    else
        out = false
    end
    return out
end


function get_constants(g_img::Matrix{Int64},loc_time::Vector{DateTime})

    drad      = [0.533 1.066 2.132]./2
    im_centre = size(g_img,1)./2

    lens_profile_tht  = (0:10:90)
    lens_profile_rpix = (0:1/9:1)

    trans_for = fill(0.0,size(loc_time,1))

    return drad, im_centre, lens_profile_tht, lens_profile_rpix, trans_for

end
