
function update_deprecated_settings!(par_in::Dict{String, Any})
    update_deprecated_settings!(par_in,Dict{String, String}())
end

function update_deprecated_settings!(par_in::Dict{String, Any},dat_in::Dict{String, String})

    # basic warnings/errors/backwards compatibilities with older settings files

    ### disable tilt option for now since it is not implemented properly
    haskey(par_in,"tilt") && delete!(par_in,"tilt")

    ### back compatibilities with old versions (just in case)
    haskey(par_in,"terrain_peri") && (par_in["lowres_peri"] = par_in["terrain_peri"]; delete!(par_in,"terrain_peri"))
    haskey(par_in,"dtm_peri") && (par_in["highres_peri"] = par_in["dtm_peri"]; delete!(par_in,"dtm_peri"))
    haskey(par_in,"surf_peri") && (par_in["forest_peri"] = par_in["surf_peri"]; delete!(par_in, "surf_peri"))
    haskey(par_in,"save_hlm") && (par_in["save_horizon"] = par_in["save_hlm"]; delete!(par_in,"save_hlm"))
    haskey(par_in,"horizon_line") && (par_in["hlm_precalc"] = par_in["horizon_line"]; delete!(par_in,"horizon_line"))

    # added epsg_code option for better compatibility with different coordinate systems (April 2025)
    if haskey(par_in,"coor_system")
        (par_in["coor_system"] == "CH1903+") && (par_in["epsg_code"] = 2056)
        (par_in["coor_system"] == "CH1903") && (par_in["epsg_code"] = 21781)

        if startswith(par_in["coor_system"], "UTM")
            @error "UTM code in coor_system setting. Update input to use epsg_code setting."
        end

        delete!(par_in, "coor_system")

    end
    !haskey(par_in,"epsg_code") && error("epsg_code must be specified in par_in")

    ### v0.12->v0.13 (November 2025)
    # replacement of "season" with "phenology"
    if haskey(par_in,"season") && par_in["season"] != "none"
        par_in["phenology"] = par_in["season"] == "summer" ? "leafon" : 
                            par_in["season"] == "winter" ? "leafoff" : 
                            par_in["season"] == "both" ? "both" : "none"
    else
        par_in["phenology"] = "none"
    end

    # shi2rad settings - replace SHI* information that is compatible with settings and therefore output files
    if haskey(par_in,"SHI_summer") # only do this if it's a shi2rad model run
        
        shi_summer = haskey(par_in,"SHI_summer") && par_in["SHI_summer"]
        shi_winter = haskey(par_in,"SHI_winter") && par_in["SHI_winter"]
        par_in["phenology"] = if shi_summer && shi_winter
            "both"
        elseif shi_summer
            "leafon"
        elseif shi_winter
            "leafoff"
        else
            "none"
        end
        (par_in["phenology"] == "none") ? (par_in["forest_type"] = "evergreen") : (par_in["forest_type"] = "deciduous")
        
        par_in["calc_terrain"] = par_in["SHI_terrain"]
        
        # Clean up old SHI settings
        for key in ["SHI_summer", "SHI_winter", "SHI_evergreen", "SHI_terrain"]
            haskey(par_in, key) && delete!(par_in, key)
        end

    end

    # replacement of "tree_species" with "leaf_type"
    haskey(par_in,"tree_species") && (par_in["leaf_type"] = par_in["tree_species"])

    # deal with input filename changes
    haskey(dat_in,"dtmf") && (dat_in["hrdtmf"] = dat_in["dtmf"]; delete!(dat_in,"dtmf")) # dtmf->hrdtmf
    haskey(dat_in,"demf") && (dat_in["lrdtmf"] = dat_in["demf"]; delete!(dat_in,"demf")) # demf->lrdtmf
    haskey(dat_in,"terf") && (dat_in["termf"] = dat_in["terf"]; delete!(dat_in,"terf")) # terf->termf
    haskey(dat_in,"dsmf") && (dat_in["lasf"] = dat_in["dsmf"]; delete!(dat_in,"dsmf"))  # dsmf->lasf

    # other old settings
    haskey(par_in,"progress") && (par_in["step_progress"] = par_in["progress"]; delete!(par_in,"progress"))  # dsmf->lasf

end

function terrain_defaults!(par_in)

    par_in["phenology"] = "none"
    par_in["leaf_type"] = "none"
    par_in["forest_type"] = "none"
    par_in["calc_terrain"] = true

end

# """
#     check_conflicts(st::SETTINGS)

# Check for conflicting settings in the user input SETTINGS struct.
# Doesn't check for everything, but checks for basic issues that will cause
#   the model to eventually fail 

# # Arguments
# - `st::SETTINGS`: The SETTINGS struct containing user-defined settings.

# """
function check_conflicts(st::SETTINGS,fp::FILEPATHS,model::String="none")

    ### create some early errors to save time

    if model == "l2r"
        # check point size length is only 2
        length(st.point_size) != 2 && 
                error("point size requires two values") 
        # and that the first value is higher than second
        (st.point_size[1] < st.point_size[2]) && 
                error("first value in point_size should be larger or equal to the second value")

        # you can't calculate terrain within l2r because I don't want to overcomplicate things
        st.calc_terrain && error("cannot save calculated terrain in L2R . . . 
                set calc_terrain=false, run ter2rad first, then rerun L2R with terrainmask_precalc=true")
    end

    if model == "c2r"
        # can't calculate both leaf-on and leaf-off when lavdf or cbh input files are used:
        if ((fp.lavdf != "") && (st.phenology == "both")) || ((fp.cbhf != "") && (st.phenology == "both"))
            error("user request for leaf on and/or leaf off phenologies is not possible with
                    only one lavdf or cbhf . . . suggest creating separate model runs or special_implementation")
        end
        # if forest type is evergreen then phenology has to be none
        (st.forest_type == "evergreen" && st.phenology != "none") &&
            error("forest_type = \"evergreen\" is only possible when phenology = \"none\"")
        # check if phenology is leafon, leafoff or both that either deciduous or mixed are the options for forest type
        (st.phenology in ("leafon","leafoff","both")) && !(st.forest_type in ("deciduous","mixed")) &&
            error("phenology = \"leafon\", \"leafoff\" or \"both\" is only possible when forest_type = \"deciduous\" or \"mixed\"") 
        # warning for special_implementation oshd
        ((st.special_implementation == "oshd") && (st.hlm_precalc)) && @error("special_implementation oshd and hlm_precalc is deprecated. Use \"swissrad\" methods instead.")
    end

    # make_pngs only works if save_images is enabled
    (st.make_pngs && !st.save_images) &&
        @warn "\"make_pngs\" requested but \"save_images\" is not enabled"

    # time_zone needed if calc_trans and calc_swr are requested
    (((st.calc_swr > 0) || !st.calc_trans) && (st.time_zone == -99)) &&
        error("\"time_zone\" must be specified if calc_trans or calc_swr enabled")

    # time_zone needed if calc_trans and calc_swr are requested
    (((st.calc_swr > 0) || !st.calc_trans) && (st.epsg_code == -99)) &&
        error("\"epsg_code\" of datasets must be specified if calc_trans or calc_swr enabled")

    # calc_trans needs to be enabled if calc_swr is enabled
    ((st.calc_swr > 0) && !st.calc_trans) &&
        (@warn "\"calc_trans\"=false in settings file; shortwave radiation will not be calculated")

end

# function extract(d::Dict)
#     expr = quote end
#     for (k, v) in d
#        push!(expr.args, :($(Symbol(k)) = $v))
#     end
#     # eval(expr)
#     return expr
# end


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


