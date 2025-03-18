function calcringratios(canrad::CANRAD,mat2ev::Matrix{Int64})

    @unpack w2all, relevant_pix, ring_radius = canrad

    fill!(w2all,NaN)
    for rix = 1:1:size(ring_radius,1)-1
        w2all[rix]   = sum(mat2ev[relevant_pix[:,rix]].==1) / sum(relevant_pix[:,rix])
    end

    return w2all

end

function calc_svf(canrad::CANRAD,mat2ev::Matrix{Int64})

    w2all = calcringratios(canrad,mat2ev)

    @unpack surf_area_p, surf_area_h = canrad

    return sum(w2all.*surf_area_p), sum(w2all.*surf_area_h)
end

function getsundxs!(window::Matrix{Float64},trans_for::Vector{Float64},
    mini::Float64,maxi::Float64,nolp::Int64,noap::Int64,ix::Int64)

    window .= min.(window,maxi)
    trans_for[ix] = (sum(mini .<= window .<= maxi)-nolp)/(length(window)-noap)
    return trans_for

end

function calc_transmissivity!(canrad::CANRAD,solar::SOLAR,trans_for::Vector{Float64},
    mat2ev::Matrix{Float64},sol_phi::Vector{Float64},sol_tht::Vector{Float64},
    imcX=0.0::Float64,imcY=0.0::Float64) # imcX and imcY are offsets from centre if 'tilt' is enabled

    @unpack lens_profile_tht, lens_profile_rpix, im_centre, radius, imsize = canrad
    @unpack prad, dpix  = solar

    fill!(prad,0)
    fill!(trans_for,0)

    keepix       = findall(sol_tht .<= 90)
    prad[keepix] = linear_interpolation(lens_profile_tht,lens_profile_rpix*radius)(sol_tht[keepix])

    x = ((im_centre .+ sin.(deg2rad.(sol_phi[keepix])) .* prad[keepix]) .- (imcX / (90 / radius)))
    y = ((im_centre .+ cos.(deg2rad.(sol_phi[keepix])) .* prad[keepix]) .- (imcY / (90 / radius)))

    x_lower = Int.(max.(1,round.(x .- dpix)))
    y_lower = Int.(max.(1,round.(y .- dpix)))
    x_upper = Int.(min.(imsize[2],round.(x .+ dpix)))
    y_upper = Int.(min.(imsize[2],round.(y .+ dpix)))

    # probably doesn't need to be looped but can be mapped instead...
    # might not make a difference in speed though
    @inbounds @simd for six in eachindex(keepix)
        getsundxs!(mat2ev[y_lower[six]:y_upper[six],x_lower[six]:x_upper[six]],
            trans_for,0.6,0.6,0,0,keepix[six])
    end

    trans_for[isnan.(trans_for)] .== 0

    return trans_for # note there is no weighting because only one solar disc size is used

end

function aggregate_data(solar::SOLAR,dat::Vector{Float64})

    @unpack agg_dat, tdx_step, t_agg_int = solar

    fill!(agg_dat,0.0)

    for tdx = 1:size(agg_dat,1)-1
        agg_dat[tdx] = mean(dat[tdx_step[tdx]:tdx_step[tdx]+t_agg_int])
    end

    agg_dat[isnan.(agg_dat)] .= -9999

    return agg_dat

end


function get_latlon(loc_x::Float64,loc_y::Float64,epsg_code::Number)

    source_epsg = epsg_code
    target_epsg = 4326

    trans = Proj.Transformation("EPSG:$source_epsg", "EPSG:$target_epsg"; always_xy=true)

    # Transform the coordinates using broadcasting
    lon, lat = trans.(loc_x, loc_y)

    return lat, lon

end

function calc_solar_track(solar::SOLAR,lat::Float64,lon::Float64,time_zone::Int64)

    @unpack true_solar_time_min, tod, eq_of_time_minutes, sun_declin_deg, idx  = solar
    @unpack hour_angle_deg, solar_zenith_angle_deg, solar_elev_angle_deg = solar
    @unpack approx_atm_refrac_deg, solar_elev_corr_atm_ref_deg, solar_azimuth_angle = solar

    # calculate solar elevation angle
    true_solar_time_min .= mod.(tod .* 1440 .+ eq_of_time_minutes .+ 4 .* lon .- 60 .* time_zone,1440)

    idx                 .= (true_solar_time_min ./4) .< 0
    hour_angle_deg[idx] .= true_solar_time_min[idx] ./4 .+ 180
    idx                 .= (true_solar_time_min ./4).>= 0
    hour_angle_deg[idx] .= true_solar_time_min[idx] ./4 .- 180

    solar_zenith_angle_deg .= rad2deg.(acos.(sin.(deg2rad.(lat)) .* sin.(deg2rad.(sun_declin_deg)) .+
                                    cos.(deg2rad.(lat)) .* cos.(deg2rad.(sun_declin_deg)) .* cos.(deg2rad.(hour_angle_deg))))
    solar_elev_angle_deg   .= 90 .- solar_zenith_angle_deg

    # calculate atmospheric diffraction
    approx_atm_refrac_deg[solar_elev_angle_deg .> 85] .= 0

    idx                        .= (5 .< solar_elev_angle_deg .<= 85)
    approx_atm_refrac_deg[idx] .= (58.1 ./ tan.(deg2rad.(solar_elev_angle_deg[idx])) .- 0.07 ./
                                    (tan.(deg2rad.(solar_elev_angle_deg[idx]))).^3 .+ 0.000086 ./
                                    (tan.(deg2rad.(solar_elev_angle_deg[idx]))).^5) ./ 3600

    idx                        .= (-0.757 .< solar_elev_angle_deg .<= 5)
    approx_atm_refrac_deg[idx] .= (1735 .+ solar_elev_angle_deg[idx] .* (-518.2 .+ solar_elev_angle_deg[idx]
                                    .* (103.4 .+ solar_elev_angle_deg[idx] .*
                                    (-12.79 .+ solar_elev_angle_deg[idx] .* 0.711)))) ./ 3600

    idx                        .= (solar_elev_angle_deg .<= -0.757)
    approx_atm_refrac_deg[idx] .= (-20.772 ./ tan.(deg2rad.(solar_elev_angle_deg[idx]))) ./ 3600

    solar_elev_corr_atm_ref_deg .= solar_elev_angle_deg .+ approx_atm_refrac_deg

    # calculate solar azimuth angle depending on hour angle
    idx                      .= hour_angle_deg .> 0
    solar_azimuth_angle[idx] .= mod.(round.((rad2deg.(acos.(round.(((sin.(deg2rad.(lat)) .*
                                    cos.(deg2rad.(solar_zenith_angle_deg[idx]))) .-
                                    sin.(deg2rad.(sun_declin_deg[idx]))) ./ (cos.(deg2rad.(lat)) .*
                                    sin.(deg2rad.(solar_zenith_angle_deg[idx]))),digits=10))) .+ 180).*100000)./100000,360)

    idx                      .= hour_angle_deg .<= 0
    solar_azimuth_angle[idx] .= mod.(round.((540 .- rad2deg.(acos.(round.(((sin.(deg2rad.(lat)) .*
                                    cos.(deg2rad.(solar_zenith_angle_deg[idx]))) .-
                                    sin.(deg2rad.(sun_declin_deg[idx]))) ./ (cos.(deg2rad.(lat)) .*
                                    sin.(deg2rad.(solar_zenith_angle_deg[idx]))),digits=10)))).*100000)./100000,360)

    # convert solar position to phi/tht
    solar_azimuth_angle[solar_azimuth_angle .>=360] .-= 360
    solar_azimuth_angle[solar_azimuth_angle .< 0] .+ 360

    # return sol_tht, sol_phi, sol_sin_elev
    return 90 .- solar_elev_corr_atm_ref_deg, solar_azimuth_angle, sin.(solar_elev_corr_atm_ref_deg./360 .*2 .* pi)

    # # convert solar position to phi/tht
    # sol_tht                     = 90 .- solar_elev_corr_atm_ref_deg          # convert to zenith angle in degree, i.e. sun @ zenith is tht = 0; sun @ horizon is tht = 90
    # sol_phi                     = solar_azimuth_angle
    # sol_phi[sol_phi.>=360]      = sol_phi[sol_phi .>= 360] .- 360
    # sol_phi[sol_phi.<0]         = sol_phi[sol_phi .< 0] .+ 360                            #convert to azimuth angle in degree, so that 0ï¿½ is North, 90ï¿½ is East, 180ï¿½ is South, and 270ï¿½ is East
    # sol_sin_elev                = sin.(solar_elev_corr_atm_ref_deg./360 .*2 .* pi) # sin of elevation angle

end

function calculateSWR(radiation::RADIATION,trans_for::Vector{Float64},sol_sinelev::Vector{Float64},
        svf::Float64,calc_swr::Int64)

    @unpack solar_const, dif_frac = radiation
    @unpack swr_open, trans_atm = radiation
    @unpack fix, dif_frac, dir_frac, trans_dirmax = radiation

    # check if swr_open is maximum potential and set trans_atm to 1 for all time steps
    if calc_swr == 1
        swr_open  .= max.(1367*sol_sinelev,0)
        # trans_atm = ones
    elseif calc_swr == 2
        # swr_open = measured values
        trans_atm .= swr_open./max.(1367*sol_sinelev,0)
    end

    fill!(dif_frac,1)
    fill!(dir_frac,1)

    fix           .= sol_sinelev .> 0
    dif_frac[fix] .= 0.165
    fix           .= (sol_sinelev .> 0) .& (trans_atm .< 0.80)
    dif_frac[fix] .= 0.9511 .- 0.1604.*trans_atm[fix] .+ 4.388.*trans_atm[fix].^2 .-
                       16.638.*trans_atm[fix].^3 .+ 12.336.*trans_atm[fix].^4
    fix           .= (sol_sinelev .> 0) .& (trans_atm .< 0.22)
    dif_frac[fix] .= 1 .- 0.09.*trans_atm[fix]

    dir_frac      .= 1 .- dif_frac

    max_dir        = min.(solar_const .* exp.(log(trans_dirmax) ./ max.(sol_sinelev,0)),
                        solar_const .* trans_atm .* dir_frac)

    # swr_for_dif      = dif_frac .* swr_open .* svf
    swr_for_dir_flat = min.(dir_frac .* swr_open,max_dir) .* trans_for

    # calculated incoming to a tilted surface (currently disabled)
    # loc_slp = NaN
    # loc_asp = NaN
    # cang_1           = max.(sin.((90 .- sol_tht)/360*2*pi),0.001)
    # temp = sin.((90 .- sol_tht)/360*2*pi).*cos.(loc_slp/360*2*pi) .+
    #                         cos.((90 .- sol_tht)/360*2*pi).*sin.(loc_slp/360*2*pi) .*
    #                         cos.((loc_asp .- sol_phi)/360*2*pi)
    # temp[isnan.(temp)] .= 0
    # cang_2           = max.(temp,0)
    # swr_for_dir_incl = min.(dir_frac .* swr_open ./ cang_1,max_dir) .* trans_for_wgt .* cang_2
    # swr_for_all_incl = swr_for_dif + swr_for_dir_incl

    return (dif_frac .* swr_open .* svf) + swr_for_dir_flat, swr_for_dir_flat
    # return swr_tot, swr_dir
end
