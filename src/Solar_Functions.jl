function calc_ringratios(ring_tht,ring_radius,mat2ev,g_rad)
    w2all     = zeros(size(ring_radius,1)-1) .* NaN
    surf_area = zeros(size(ring_radius,1)-1) .* NaN
    cos_corr  = zeros(size(ring_radius,1)-1) .* NaN

    for rix = 1:1:size(ring_radius,1)-1
        inner_radius   = ring_radius[rix]
        outer_radius   = ring_radius[rix+1]
        relevant_pix   = findall(inner_radius .< vec(g_rad) .< outer_radius)
        w2all[rix]     = sum(mat2ev[relevant_pix].==1) / size(relevant_pix,1)
        surf_area[rix] = 2*pi*(cos(ring_tht[rix]/360*2*pi) - cos(ring_tht[rix+1]/360*2*pi))/2/pi
        cos_corr[rix]  = (cos(ring_tht[rix]/360*2*pi) + cos(ring_tht[rix+1]/360*2*pi))/2
    end

    return w2all, surf_area, cos_corr
end

function calculateVf(mat2ev,g_rad)

    lens_profile_tht  = collect(0:10:90)
    lens_profile_rpix = collect(0:1/9:1)
    ring_tht          = collect(0:90/9:90)
    ring_radius       = pyinterp.interp1d(lens_profile_tht,lens_profile_rpix*radius)(ring_tht)

    w2all, surf_area, cos_corr = calc_ringratios(ring_tht,ring_radius,mat2ev,g_rad)

    Vf = zeros(1,2) .* NaN
    Vfweighted = sum(w2all.*surf_area.*cos_corr) / sum(surf_area.*cos_corr)
    Vfflat = sum(w2all.*surf_area)

    return Vfweighted, Vfflat
end

function calcmm(mat2ev,radius,drad,six,dix,x,y)
        dpix = drad[dix] / 90*radius*sqrt(pi)/2
        xmm  = Int.(collect(max(1,round(x[six]-dpix)):1:min(size(mat2ev,2),round(x[six]+dpix))))
        ymm  = Int.(collect(max(1,round(y[six]-dpix)):1:min(size(mat2ev,1),round(y[six]+dpix))))
    return xmm,ymm
end

function getsundxs(mat2ev,trans_for,xmm,ymm,mini,maxi,nolp,noap,dix,keepix,six)
        mat2ev[ymm,xmm] = min.(mat2ev[ymm,xmm],maxi)
        nolp1 = sum(mini .<= mat2ev[ymm,xmm] .<= maxi)
        noap1 = length(mat2ev[ymm,xmm])
        trans_for[keepix[six],dix] = (nolp1-nolp)/(noap1-noap)
    return mat2ev,nolp1,noap1,trans_for
end

function calculateSWR(mat2ev,loc_time,radius,sol_phi,sol_tht,Vf,g_coorpol,sol_sinelev)
    im_centre         = size(mat2ev)./2
    drad              = [0.533 1.066 2.132]./2
    trans_for         = zeros(size(loc_time,1),size(drad,2))

    lens_profile_tht  = collect(0:10:90)
    lens_profile_rpix = collect(0:1/9:1)
    dx                = findall(sol_tht .<= maximum(filter(!isnan,g_coorpol[:,2])))
    prad              = ones(size(sol_tht,1)) .* NaN
    prad[dx]          = pyinterp.interp1d(lens_profile_tht,lens_profile_rpix*radius)(sol_tht[dx])

    keepix            = findall(sol_tht .<= 90)
    x                 = im_centre[1] .+ sin.(deg2rad.(sol_phi[keepix])) .* prad[keepix]
    y                 = im_centre[2] .+ cos.(deg2rad.(sol_phi[keepix])) .* prad[keepix]

    for six = 1:1:size(x,1)
        # global mat2ev, trans_for
        xmm, ymm = calcmm(mat2ev,radius,drad,six,1,x,y)
        mat2ev,nolp1,noap1,trans_for = getsundxs(mat2ev,trans_for,xmm,ymm,0.6,0.6,0,0,1,keepix,six)

        xmm, ymm = calcmm(mat2ev,radius,drad,six,2,x,y)
        mat2ev,nolp2,noap2,trans_for = getsundxs(mat2ev,trans_for,xmm,ymm,0.6,0.61,nolp1,noap1,2,keepix,six)

        xmm, ymm = calcmm(mat2ev,radius,drad,six,3,x,y)
        mat2ev,_,_,trans_for = getsundxs(mat2ev,trans_for,xmm,ymm,0.61,0.62,nolp2,noap2,3,keepix,six)
    end

    trans_for[isnan.(trans_for)] .== 0
    trans_for_wgt = (trans_for[:,1].*60 + trans_for[:,2]*30 + trans_for[:,3].*10)./100

    swr_open =  max.(1367*sol_sinelev,0)

    trans_atm      = swr_open./max.(1367*sol_sinelev,0)
    dif_frac       = ones(size(loc_time))

    fix            = findall(sol_sinelev .> 0)
    dif_frac[fix]  .= 0.165
    fix            = (sol_sinelev .> 0) .& (trans_atm .< 0.80)
    dif_frac[fix]  .= 0.9511 .- 0.1604.*trans_atm[fix] .+ 4.388.*trans_atm[fix].^2 .-
                        16.638.*trans_atm[fix].^3 .+ 12.336.*trans_atm[fix].^4
    fix            = (sol_sinelev .> 0) .& (trans_atm .< 0.22)
    dif_frac[fix]  = 1 .- 0.09.*trans_atm[fix]
    dir_frac       = 1 .- dif_frac

    trans_dirmax = 1 - 0.165
    max_dir_1    = 1367 .* exp.(log(trans_dirmax) ./ max.(sol_sinelev,0))
    max_dir_2    = 1367 .* trans_atm .* dir_frac
    max_dir      = min.(max_dir_1,max_dir_2)

    swr_for_dif    = dif_frac .* swr_open .* Vf

    loc_slp = NaN
    loc_asp = NaN

    cang_1           = max.(sin.((90 .- sol_tht)/360*2*pi),0.001)
    temp = sin.((90 .- sol_tht)/360*2*pi).*cos.(loc_slp/360*2*pi) .+
                            cos.((90 .- sol_tht)/360*2*pi).*sin.(loc_slp/360*2*pi) .*
                            cos.((loc_asp .- sol_phi)/360*2*pi)
    temp[isnan.(temp)] .= 0
    cang_2           = max.(temp,0)
    swr_for_dir_flat = min.(dir_frac .* swr_open,max_dir) .* trans_for_wgt
    swr_for_dir_incl = min.(dir_frac .* swr_open ./ cang_1,max_dir) .* trans_for_wgt .* cang_2

    swr_for_all_flat = swr_for_dif + swr_for_dir_flat
    swr_for_all_incl = swr_for_dif + swr_for_dir_incl
    return swr_for_all_flat, swr_for_dir_flat, swr_for_dif, trans_for_wgt
end

function utm2deg(loc_x,loc_y,zone,hemi)
    # Credit: based on utm2lonlat found on Matlab file exchange; accessed on 2018/08/25
    # Copyright (c) 2013, Erwin N. All rights reserved.
    sa    = 6378137.000000
    sb    = 6356752.314245
    e2sq = ((((sa .^ 2 ) - (sb .^ 2 )).^ 0.5) ./ sb).^2
    c     = ( sa .^ 2 ) ./ sb;

    if @match(hemi,"N")
        X    = loc_x - 500000
        Y    = loc_y
    elseif @match(hemi,"S")
        X    = loc_x - 500000
        Y    = loc_y - 10000000
    end

    S      = ( ( zone .* 6 ) - 183 )
    lat    =  Y ./ ( 6366197.724 .* 0.9996 )
    v      = ( c ./ ( ( 1 + ( e2sq .* ( cos(lat) ) .^ 2 ) ) ) .^ 0.5 ) .* 0.9996
    a      = X ./ v
    a1     = sin( 2 .* lat )
    a2     = a1 .* ( cos(lat) ) .^ 2
    j2     = lat + ( a1 ./ 2 )
    j4     = ( ( 3 .* j2 ) + a2 ) ./ 4
    j6     = ( ( 5 .* j4 ) + ( a2 .* ( cos(lat) ) .^ 2) ) ./ 3
    alpha  = ( 3 ./ 4 ) .* e2sq
    beta   = ( 5 ./ 3 ) .* alpha .^ 2
    gamma  = ( 35 ./ 27 ) .* alpha .^ 3
    Bm     = 0.9996 .* c .* ( lat - alpha .* j2 + beta .* j4 - gamma .* j6 )
    b      = ( Y - Bm ) ./ v
    Epsi   = ( ( e2sq .* a.^2 ) ./ 2 ) .* ( cos(lat) ).^ 2
    Eps    = a .* ( 1 - ( Epsi ./ 3 ) )
    nab    = ( b .* ( 1 - Epsi ) ) + lat
    senoheps = ( exp(Eps) - exp(-Eps) ) ./ 2
    Delt   = atan(senoheps ./ (cos(nab) ) )
    TaO    = atan(cos(Delt) .* tan(nab))
    longitude = (Delt .* (180/pi) ) + S

    latitude = ( lat + ( 1 + e2sq .* (cos(lat).^2) - ( 3/2 )
      .* e2sq .* sin(lat) .* cos(lat) .* ( TaO - lat ) )
      .* ( TaO - lat ) ) .* (180/pi)

    return latitude, longitude
end

function calc_solar_track(pts,loc_time,time_zone,coor_system)

    loc_x = mean([maximum(filter(!isnan,pts[:,1])),minimum(filter(!isnan,pts[:,1]))])
    loc_y = mean([maximum(filter(!isnan,pts[:,2])),minimum(filter(!isnan,pts[:,2]))])

    if @match(coor_system,"CH1903")
        xd  = (loc_x - 600000)/1000000
        yd  = (loc_y - 200000)/1000000
        lon = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36
        lat = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36
    elseif @match(coor_system,"CH1903+")
        xd  = (loc_x - 2600000)/1000000
        yd  = (loc_y - 1200000)/1000000
        lon = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36
        lat = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36
    elseif @match(coor_system,"UTM")
        utmz     = map(join, groupby(identity, utm_zone))
        lat,lon  = utm2deg(loc_x,loc_y,parse(Float64,utmz[1]),utmz[2])
    end

    # time conversion
    doy = Dates.dayofyear.(loc_time)
    tod = (Dates.hour.(loc_time) + (Dates.minute.(loc_time)./60) + Dates.second.(loc_time)) ./ 24
    jd = Dates.datetime2julian.(loc_time) .- time_zone/24
    jc = ((Dates.datetime2julian.(loc_time) .- time_zone/24) .- 2451545) ./ 36525

    # calculate solar elevation angle
    geom_mean_long_sun_deg    = mod.(280.46646 .+ jc .* (36000.76983 .+ jc .* 0.0003032),360)
    geom_mean_anom_sun_deg    = 357.52911 .+ jc .* (35999.05029 .- 0.0001537 .* jc)
    eccent_earth_orbit        = 0.016708634 .- jc .* (0.000042037 .+ 0.0000001267 .* jc)
    sun_eq_of_ctr             = sin.(deg2rad.(geom_mean_anom_sun_deg)) .* (1.914602 .- jc .* (0.004817 .+ 0.000014 .* jc)) .+
                                    sin.(deg2rad.(2 .* geom_mean_anom_sun_deg)) .* ( 0.019993 .- 0.000101 .* jc) .+
                                     sin.(deg2rad.(3 .* geom_mean_anom_sun_deg)) .* 0.000289
    sun_true_long_deg         = sun_eq_of_ctr .+ geom_mean_long_sun_deg
    sun_app_long_deg          = sun_true_long_deg .- 0.00569 .- 0.00478 .* sin.(deg2rad.(125.04 .- 1934.136 .* jc))
    mean_obliq_ecliptic_deg   = 23 .+ (26 .+ ((21.448 .- jc .* (46.815 .+ jc .* (0.00059 .- jc .* 0.001813)))) ./ 60) ./ 60
    obliq_corr_deg            = mean_obliq_ecliptic_deg .+ 0.00256 .* cos.(deg2rad.(125.04 .- 1934.136 .* jc))
    sun_declin_deg            = rad2deg.(asin.(sin.(deg2rad.(obliq_corr_deg)) .* sin.(deg2rad.(sun_app_long_deg))))
    var_y                     = tan.(deg2rad.(obliq_corr_deg ./ 2)) .* tan.(deg2rad.(obliq_corr_deg ./ 2))
    eq_of_time_minutes        = 4 .* rad2deg.(var_y .* sin.(2 .* deg2rad.(geom_mean_long_sun_deg)) .- 2 .* eccent_earth_orbit .*
                                    sin.(deg2rad.(geom_mean_anom_sun_deg)) .+ 4 .* eccent_earth_orbit .* var_y .*
                                     sin.(deg2rad.(geom_mean_anom_sun_deg)) .* cos.(2 .* deg2rad.(geom_mean_long_sun_deg)) .-
                                      0.5 .* var_y .* var_y .* sin.(4 .* deg2rad.(geom_mean_long_sun_deg)) .- 1.25 .*
                                       eccent_earth_orbit .* eccent_earth_orbit .* sin.(2 .* deg2rad.(geom_mean_anom_sun_deg)))
    true_solar_time_min       = mod.(tod .* 1440 .+ eq_of_time_minutes .+ 4 .* lon .- 60 .* time_zone,1440)

    hour_angle_deg            = zeros(size(true_solar_time_min)) .* NaN
    idx1                      = (Int.(true_solar_time_min ./4 .< 0))
    hour_angle_deg[idx1.==1] .= true_solar_time_min[idx1.==1] ./4 .+ 180
    idx2                      = (Int.(true_solar_time_min ./4 .>= 0))
    hour_angle_deg[idx2.==1] .= true_solar_time_min[idx2.==1] ./4 .- 180

    solar_zenith_angle_deg    = rad2deg.(acos.(sin.(deg2rad.(lat)) .* sin.(deg2rad.(sun_declin_deg)) .+
                                    cos.(deg2rad.(lat)) .* cos.(deg2rad.(sun_declin_deg)) .* cos.(deg2rad.(hour_angle_deg))))
    solar_elev_angle_deg      = 90 .- solar_zenith_angle_deg

    # calculate atmospheric diffraction
    approx_atm_refrac_deg       = zeros(size(solar_elev_angle_deg)).* NaN
    approx_atm_refrac_deg[solar_elev_angle_deg .> 85] .= 0

    idx1                             = Int.(5 .< solar_elev_angle_deg .<= 85)
    approx_atm_refrac_deg[idx1.==1] .= (58.1 ./ tan.(deg2rad.(solar_elev_angle_deg[idx1.==1])) .- 0.07 ./
                                            (tan.(deg2rad.(solar_elev_angle_deg[idx1.==1]))).^3 .+ 0.000086 ./
                                             (tan.(deg2rad.(solar_elev_angle_deg[idx1.==1]))).^5) ./ 3600
    idx2                             = Int.(-0.757 .< solar_elev_angle_deg .<= 5)
    approx_atm_refrac_deg[idx2.==1] .= (1735 .+ solar_elev_angle_deg[idx2.==1] .* (-518.2 .+ solar_elev_angle_deg[idx2.==1]
                                            .* (103.4 .+ solar_elev_angle_deg[idx2.==1] .*
                                             (-12.79 .+ solar_elev_angle_deg[idx2.==1] .* 0.711)))) ./ 3600
    idx3                             = Int.(solar_elev_angle_deg .<= -0.757)
    approx_atm_refrac_deg[idx3.==1]  = (-20.772 ./ tan.(deg2rad.(solar_elev_angle_deg[idx3.==1]))) ./ 3600
    solar_elev_corr_atm_ref_deg      = solar_elev_angle_deg .+ approx_atm_refrac_deg

    # calculate solar azimuth angle depending on hour angle
    solar_azimuth_angle            = zeros(size(hour_angle_deg)) .* NaN
    idx1                           = hour_angle_deg .> 0
    solar_azimuth_angle[idx1.==1] .= mod.(round.((rad2deg.(acos.(((sin.(deg2rad.(lat)) .*
                                        cos.(deg2rad.(solar_zenith_angle_deg[idx1.==1]))) .-
                                         sin.(deg2rad.(sun_declin_deg[idx1.==1]))) ./ (cos.(deg2rad.(lat)) .*
                                          sin.(deg2rad.(solar_zenith_angle_deg[idx1.==1]))))) .+ 180).*100000)./100000,360)
    idx2                           = hour_angle_deg .<= 0
    solar_azimuth_angle[idx2.==1] .= mod.(round.((540 .- rad2deg.(acos.(((sin.(deg2rad.(lat)) .*
                                        cos.(deg2rad.(solar_zenith_angle_deg[idx2.==1]))) .-
                                         sin.(deg2rad.(sun_declin_deg[idx2.==1]))) ./ (cos.(deg2rad.(lat)) .*
                                          sin.(deg2rad.(solar_zenith_angle_deg[idx2.==1])))))).*100000)./100000,360)

    # convert solar position to phi/tht
    sol_tht                     = 90 .- solar_elev_corr_atm_ref_deg          # convert to zenith angle in degree, i.e. sun @ zenith is tht = 0; sun @ horizon is tht = 90
    sol_phi                     = solar_azimuth_angle
    sol_phi[sol_phi.>=360]      = sol_phi[sol_phi .>= 360] .- 360
    sol_phi[sol_phi.<0]         = sol_phi[sol_phi .< 0] .+ 360                            #convert to azimuth angle in degree, so that 0ï¿½ is North, 90ï¿½ is East, 180ï¿½ is South, and 270ï¿½ is East
    sol_sin_elev                = sin.(solar_elev_corr_atm_ref_deg./360 .*2 .* pi) # sin of elevation angle

    return sol_tht, sol_phi, sol_sin_elev
end
