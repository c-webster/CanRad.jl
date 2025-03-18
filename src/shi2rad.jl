"""

! calc_trans => "true" is assumed


"""

function shi2rad!(shif::String,par_in_shi::Dict{String, Any},
    exdir_shi::String,taskID="task",pts_m=nothing)

    eval(extract(par_in_shi))

    canrad  = CANRAD()

    # get model number (1 = canopy, 2 = terrain only)
    pts_m == nothing && (pts_m = ones(size(pts_x)))

    # load the shis and get pts from the file
    shi_ds = NCDataset(shif,"r")

    pts_x = shi_ds["easting"][:]
    pts_y = shi_ds["northing"][:]

    outdir, outstr  = organise_outf(taskID,exdir_shi,batch)
    global outtext = "Processing "*cfmt.("%.$(0)f", 0)*"% ... "*string(0)*" of "*string(size(pts_x,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    # get time/date info
    loc_time     = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
    loc_time_agg = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))

    dataset = createfiles_fromSHI(outdir,outstr,hcat(pts_x,pts_y),calc_swr,
                                        SHI_summer,SHI_winter,SHI_terrain,SHI_evergreen,loc_time_agg,time_zone)
    pts_lat, pts_lon = get_latlon(pts_x,pts_y,epsg_code)
    solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = tstep, radius = canrad.radius, time_zone = time_zone)

    @unpack radius, g_rad = canrad
    g_rad[g_rad .> radius] .= NaN

    @unpack ring_radius, ring_tht, surf_area_p, surf_area_h, relevant_pix = canrad
    for rix = 1:1:size(ring_radius,1)-1
        relevant_pix[:,rix] = (ring_radius[rix] .< vec(g_rad) .< ring_radius[rix+1])
        surf_area_p[rix] = sin(ring_tht[rix+1]/360*2*pi)^2 - sin(ring_tht[rix]/360*2*pi)^2
        surf_area_h[rix] = 2*pi*(cos(ring_tht[rix]/360*2*pi) - cos(ring_tht[rix+1]/360*2*pi))/2/pi
    end

    @simd for crx = 1:size(pts_x,1)

        sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],time_zone)
        @unpack trans_for = solar

        if SHI_summer && (pts_m[crx] .== 1)

            mat2ev = Int64.(shi_ds["SHI_summer"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_s"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_s"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["for_trans_s"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,calc_swr)
                dataset["swr_total_s"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_s"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if SHI_winter && (pts_m[crx] .== 1)

            mat2ev = Int64.(shi_ds["SHI_winter"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_w"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_w"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["for_trans_w"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,calc_swr)
                dataset["swr_total_w"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_w"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if SHI_terrain

            mat2ev = Int64.(shi_ds["SHI_terrain"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_t"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_t"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["trans_t"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,calc_swr)
                dataset["swr_total_t"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_t"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if SHI_evergreen && (pts_m[crx] .== 1)

            mat2ev = Int64.(shi_ds["SHI_evergreen"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_e"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_e"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["for_trans_e"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,calc_swr)
                dataset["swr_total_e"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_e"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        # save the progress
        percentdone = Int(floor((crx / size(pts_x,1)) * 100))
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*cfmt.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts_x,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end

    close(dataset)
    close(shi_ds)

end