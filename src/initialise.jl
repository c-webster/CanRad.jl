
function compatability_check!(par_in::Dict)

    if !haskey(par_in,"buildings") # include buildings as solid objects (not floating carpets)
        par_in["buildings"] = false
    end

    if !haskey(par_in,"pt_corr") # correct for probability of transmissivity
        par_in["pt_corr"] = true
    end

    # renamed terrain options
    if haskey(par_in,"terrain")
        if par_in["terrain"] == 1
            par_in["terrain_tile"] = false
        elseif par_in["terrain"] == 2
            par_in["terrain_tile"] = true
        end
    end

    if haskey(par_in,"tershad")
        if par_in["tershad"] == 1
            par_in["terrain_highres"] = true
            par_in["terrain_lowres"] = true
        elseif par_in["tershad"] == 2
            par_in["terrain_highres"] = true
            par_in["terrain_lowres"] = false
        elseif par_in["tershad"] == 3
            par_in["terrain_highres"] = false
            par_in["terrain_lowres"] = false
        end
    end

    if !haskey(par_in,"horizon_line")
        par_in["horizon_line"] = false
    end

    # disable tilt option since it is not implemented properly
    if !haskey(par_in,"tilt") || par_in["tilt"]
        par_in["tilt"] = false
    end

    if !haskey(par_in,"OSHD")
        par_in["OSHD"] = false
    end

	# renamed variables with version 0.7.0
	if !haskey(par_in,"image_height")
		par_in["image_height"] = par_in["ch"]
	end

	if haskey(par_in,"utm_zone") && par_in["coor_system"] == "UTM"
		par_in["coor_system"] = "UTM"*" "*utm_zone
	end

	if !haskey(par_in,"b_space") && !haskey(par_in,"branch_spacing") # spacing between branch points in cm
        par_in["branch_spacing"] = 0.1
	elseif @isdefined(b_space)
		par_in["branch_spacing"] = b_space
    end

    if !haskey(par_in,"season")
        par_in["season"] = "evergreen"
    end

    if haskey(par_in,"tolerance")
        if length(par_in["tolerance"]) < 2
            error("update tolerance value") # accomodate tolerance changing to pixels
        end
    end

end

function extract(d::Dict)
    expr = quote end
    for (k, v) in d
       push!(expr.args, :($(Symbol(k)) = $v))
    end
    # eval(expr)
    return expr
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


