@with_kw struct CANRAD

    # image parameters
    radius::Int64 = 500
    diameter::Int64 = radius*2
    mat2ev::Matrix{Int64}   = ones(radius*2,radius*2)
    imsize::Tuple{Int64, Int64} = size(mat2ev)

    tgrid1::Matrix{Float64} = ((1:radius*2)' .* ones(radius*2))
    tgrid2::Matrix{Float64} = (ones(radius*2)' .* (1:radius*2))

    g_rad::Matrix{Float64}  = sqrt.((tgrid1 .- radius .- 1).^2 .+ (tgrid2 .- radius .- 1).^2)

    g_tht::Matrix{Float64}  = (g_rad./radius).*90

    g_phi::Matrix{Float64}  = reverse(map(atan,(tgrid2.-radius),(tgrid1.-radius)),dims=1)

    # g_coorpol::Matrix{Float64} = [vec(g_phi) vec(g_tht)]
    g_coorcrt::Matrix{Float64} = [vec(tgrid1) vec(tgrid2)]

    lia::BitVector = Bool.(zeros(size(g_coorcrt,1)))
    kdtreedims::Int64  = size(g_coorcrt,1)

    # parameters for Vf calculations
    im_centre::Float64 = radius

    lens_profile_tht::Vector{Int64}    = collect(0:10:90)
    lens_profile_rpix::Vector{Float64} = collect(0:1/9:1)

    ring_tht::Vector{Float64}    = collect(0:90/9:90)
    ring_radius::Vector{Float64} = pyinterp.interp1d(lens_profile_tht,lens_profile_rpix*radius)(ring_tht)

    w2all::Vector{Float64}       = fill(NaN,(size(ring_radius,1)-1))
    surf_area_p::Vector{Float64} = fill(NaN,(size(ring_radius,1)-1))
    surf_area_h::Vector{Float64} = fill(NaN,(size(ring_radius,1)-1))

    relevant_pix::BitMatrix  = Bool.(zeros(length(g_rad),(size(ring_radius,1)-1)))

    limits_canopy::Vector{Float64} = zeros(4)


end


@with_kw struct SOLAR

    # parameters for calc solar track

    ### parameters and constants for solar track and transmissivity calculations
    loc_time::Vector{DateTime}
    trans_for::Vector{Float64} = zeros(size(loc_time,1))
    time_zone::Int64

    loc_time_agg::Vector{DateTime}
    tstep::Int64
    agg_dat::Vector{Float64} = fill(0.0,(size(loc_time_agg,1)))
    tdx_step::StepRange{Int64, Int64} = (1:Int.((tstep/2)):size(loc_time,1))
    t_agg_int::Int64 = div(tstep,2)-1

    prad::Vector{Float64} = zeros(size(loc_time,1))

    radius::Int64
    # drad::Matrix{Float64} = [0.533,1.066,2.132]./2
    drad::Float64 = 0.2665 # = apparent radius of solar disc in degrees
    dpix::Float64 = drad / 90*radius*sqrt(pi)/2

    # time conversion
    doy::Vector{Int64}   = Dates.dayofyear.(loc_time)
    tod::Vector{Float64} = (Dates.hour.(loc_time) + (Dates.minute.(loc_time)./60) + Dates.second.(loc_time)) ./ 24
    jd::Vector{Float64}  = Dates.datetime2julian.(loc_time) .- time_zone/24
    jc::Vector{Float64}  = ((Dates.datetime2julian.(loc_time) .- time_zone/24) .- 2451545) ./ 36525

    # pre-calc of solar elevation angle
    geom_mean_long_sun_deg::Vector{Float64}  = mod.(280.46646 .+ jc .* (36000.76983 .+ jc .* 0.0003032),360)
    geom_mean_anom_sun_deg::Vector{Float64}  = 357.52911 .+ jc .* (35999.05029 .- 0.0001537 .* jc)
    eccent_earth_orbit::Vector{Float64}      = 0.016708634 .- jc .* (0.000042037 .+ 0.0000001267 .* jc)
    sun_eq_of_ctr::Vector{Float64}           = sin.(deg2rad.(geom_mean_anom_sun_deg)) .* (1.914602 .- jc .* (0.004817 .+ 0.000014 .* jc)) .+
                                                    sin.(deg2rad.(2 .* geom_mean_anom_sun_deg)) .* ( 0.019993 .- 0.000101 .* jc) .+
                                                    sin.(deg2rad.(3 .* geom_mean_anom_sun_deg)) .* 0.000289
    sun_true_long_deg::Vector{Float64}       = sun_eq_of_ctr .+ geom_mean_long_sun_deg
    sun_app_long_deg::Vector{Float64}        = sun_true_long_deg .- 0.00569 .- 0.00478 .* sin.(deg2rad.(125.04 .- 1934.136 .* jc))
    mean_obliq_ecliptic_deg::Vector{Float64} = 23 .+ (26 .+ ((21.448 .- jc .* (46.815 .+ jc .* (0.00059 .- jc .* 0.001813)))) ./ 60) ./ 60
    obliq_corr_deg::Vector{Float64}          = mean_obliq_ecliptic_deg .+ 0.00256 .* cos.(deg2rad.(125.04 .- 1934.136 .* jc))
    sun_declin_deg::Vector{Float64}          = rad2deg.(asin.(sin.(deg2rad.(obliq_corr_deg)) .* sin.(deg2rad.(sun_app_long_deg))))
    var_y::Vector{Float64}                   = tan.(deg2rad.(obliq_corr_deg ./ 2)) .* tan.(deg2rad.(obliq_corr_deg ./ 2))
    eq_of_time_minutes::Vector{Float64}      = 4 .* rad2deg.(var_y .* sin.(2 .* deg2rad.(geom_mean_long_sun_deg)) .- 2 .* eccent_earth_orbit .*
                                                  sin.(deg2rad.(geom_mean_anom_sun_deg)) .+ 4 .* eccent_earth_orbit .* var_y .*
                                                  sin.(deg2rad.(geom_mean_anom_sun_deg)) .* cos.(2 .* deg2rad.(geom_mean_long_sun_deg)) .-
                                                  0.5 .* var_y .* var_y .* sin.(4 .* deg2rad.(geom_mean_long_sun_deg)) .- 1.25 .*
                                                  eccent_earth_orbit .* eccent_earth_orbit .* sin.(2 .* deg2rad.(geom_mean_anom_sun_deg)))

    idx::BitVector = Bool.(zeros(size(loc_time)))
    true_solar_time_min::Vector{Float64}   = zeros(size(loc_time,1))

    hour_angle_deg::Vector{Float64}         = zeros(size(loc_time,1)) .* NaN
    solar_zenith_angle_deg::Vector{Float64} = zeros(size(loc_time,1)) .* NaN
    solar_elev_angle_deg::Vector{Float64}   = zeros(size(loc_time,1)) .* NaN

    approx_atm_refrac_deg::Vector{Float64}       = zeros(size(loc_time,1)) .* NaN
    solar_elev_corr_atm_ref_deg::Vector{Float64} = zeros(size(loc_time,1)) .* NaN
    solar_azimuth_angle::Vector{Float64}         = zeros(size(loc_time,1)) .* NaN

end

@with_kw struct RADIATION

    loc_time::Vector{DateTime} 
    solar_const::Int64 = 1367

    swr_open::Vector{Float64}  = zeros(size(loc_time))
    trans_atm::Vector{Float64} = ones(size(loc_time))

    trans_dirmax::Float64 = 1 - 0.165

    dif_frac::Vector{Float64}  = ones(size(loc_time))
    dir_frac::Vector{Float64}  = ones(size(loc_time))

    fix::BitVector = Bool.(zeros(size(loc_time)))

end

@with_kw struct TER2RAD

    pts_sz::Int64

    limits_highres::Vector{Float64} = zeros(4)
    pts_e::Vector{Float64} = zeros(pts_sz)

    limits_lowres::Vector{Float64} = zeros(4)
    pts_e_dem::Vector{Float64} = zeros(pts_sz)

    # terrain horizon line parameters
    phi_bins::Vector{Float64} = collect(-2*pi:pi/180:2*pi)

    fix1::Vector{Int64} = zeros(10000)
    tdx::Vector{Bool} = zeros(size(phi_bins,1))
    mintht::Vector{Float64} = fill(90.0,size(phi_bins,1))
    tempmintht::Vector{Float64} = fill(90.0,size(phi_bins,1))

    rphi::Vector{Float64} = collect(-pi:pi/1080:pi)
    rtht::Vector{Float64} = zeros(size(rphi,1))

    dx1::Int64 = findall(abs.(phi_bins .+ pi) .== 0)[1]
    dx2::Int64 = findall(abs.(phi_bins .- pi) .== 0)[1]

end

@with_kw struct CHM2RAD

    ### for the canopy horizon line + ptrans calcs
    phi_bins_long::Vector{Float64}  = collect(-2*pi:pi/360:2*pi)
    phi_bins_short::Vector{Float64} = collect(-pi:pi/360:pi)

    dx1_pbl::Int64 = findall(abs.(phi_bins_long .+ pi) .== 0)[1]
    dx2_pbl::Int64 = findall(abs.(phi_bins_long .- pi) .== 0)[1]

    tdx_chm::Vector{Bool} = zeros(size(phi_bins_long,1))
    tdx_bse::Vector{Bool} = zeros(size(phi_bins_short,1))

    chm_temptht_long::Vector{Float64}  = fill(90.0,size(phi_bins_long,1))
    chm_temptht_short::Vector{Int64} = fill(90.0,size(phi_bins_short,1))
    bse_temptht::Vector{Float64}       = fill(90.0,size(phi_bins_short,1))

    temp_lavd::Vector{Float64}   = zeros(size(phi_bins_long,1))
    temp_thick::Matrix{Float64}  = zeros(size(chm_temptht_short,1),90);

    sum_thick::Matrix{Float64}      = zeros((size(phi_bins_short,1),90))
    sum_lavd_thick::Matrix{Float64} = zeros(size(phi_bins_short,1),90)

    Ptrans::Matrix{Float64} = zeros(size(sum_lavd_thick))

    canopy_thick::Matrix{Float64} = zeros(size(sum_lavd_thick))

end

@with_kw struct LAS2RAD

    forest_peri::Float64
    
    rbins::Vector{Float64} = collect(0:(forest_peri-0)/5:forest_peri)

    point_size::Vector{Float64}
    knum::Vector{Float64} = collect(point_size[1]:(diff(point_size)/5)[1]:point_size[2])

    knum_t::Vector{Float64} = collect(15:(diff([10,1])/5)[1]:1)

    # eventually default parameters for the trunks instead of pre-allo trunks
end

