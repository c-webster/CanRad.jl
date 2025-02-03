
function compatability_check!(par_in::Dict)

    # basic warnings/errors
    # write check for presence of required user settings and datasets

    # check point size has two values
    (haskey(par_in,"point_size") && length(par_in["point_size"]) < 2) && error("point size requires two values") 


    ### fill default settings
    !haskey(par_in,"terrainmask_precalc") && (par_in["terrainmask_precalc"] = false)
    !haskey(par_in,"calc_terrain") && (par_in["calc_terrain"] = false)
    !haskey(par_in,"buildings") && (par_in["buildings"] = false)
    !haskey(par_in,"hlm_precalc") && (par_in["hlm_precalc"] = false)

    !haskey(par_in,"save_images") && (par_in["save_images"] = false)
    !haskey(par_in,"save_horizon") && (par_in["save_horizon"] = false)
    !haskey(par_in,"progress") && (par_in["progress"] = false)

    !haskey(par_in,"branch_spacing") && (par_in["branch_spacing"] = 0.1)

    !haskey(par_in,"cbh") && (par_in["cbh"] = 2.0)
    par_in["cbh"] = Float64(par_in["cbh"])

    ### disable tilt option since it is not implemented properly
    (!haskey(par_in,"tilt") || par_in["tilt"]) && (par_in["tilt"] = false)


    ### back compatibility with perimeter variable names in <v0.10
    haskey(par_in,"terrain_peri") && (par_in["lowres_peri"] = par_in["terrain_peri"])

    haskey(par_in,"dtm_peri") && (par_in["highres_peri"] = par_in["dtm_peri"])
    !haskey(par_in,"highres_peri") && (par_in["highres_peri"] = 300)

    haskey(par_in,"surf_peri") && (par_in["forest_peri"] = par_in["surf_peri"])

    haskey(par_in,"save_hlm") && (par_in["save_horizon"] = par_in["save_hlm"])

    !haskey(par_in,"oshd_flag") && (par_in["oshd_flag"] = false)

    (!haskey(par_in,"make_pngs") && par_in["save_images"]) && (par_in["make_pngs"] = true)

    !haskey(par_in,"keep_points") && (par_in["keep_points"] = "canopy")


end

function extract(d::Dict)
    expr = quote end
    for (k, v) in d
       push!(expr.args, :($(Symbol(k)) = $v))
    end
    # eval(expr)
    return expr
end


function check_output(exdir::String,pts::Matrix{Float64},batch::Bool,taskID::String)

    if batch
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3]
        outdir = joinpath(exdir,outstr)
    else
        outstr = String(split(exdir,"/")[end])
        outdir = exdir
    end

    if ispath(outdir)
        try
            crxstart = parse(Int,(split(reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])))[4])
            (crxstart == size(pts,1)) ? (out = true) : (out = false)
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


