"""

! calc_trans => "true" is assumed


"""

function shi2rad!(shif::String,par_in::Dict{String, Any},
    exdir::String,taskID="task",pts_m=nothing)

    # eval(extract(par_in_shi))
    update_deprecated_settings!(par_in)

    st = SETTINGS(; (Symbol.(keys(par_in)) .=> values(par_in))... )

    canrad = CANRAD()
    fileio = FILEIO(time_zone=st.time_zone,calc_swr=st.calc_swr)

    # load the shis and get pts from the file
    shi_ds = NCDataset(shif,"r")

    pts_x = shi_ds["easting"][:]
    pts_y = shi_ds["northing"][:]

    # get model number (1 = canopy, 2 = terrain only)
    pts_m == nothing && (pts_m = ones(size(pts_x)))

    outdir, outstr  = organise_outf(taskID,exdir,st.batch)
    global outtext = "Processing "*cfmt.("%.$(0)f", 0)*"% ... "*string(0)*" of "*string(size(pts_x,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    # get time/date info
    loc_time     = collect(Dates.DateTime(st.t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(st.t2,"dd.mm.yyyy HH:MM:SS"))
    loc_time_agg = collect(Dates.DateTime(st.t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(st.tstep):Dates.DateTime(st.t2,"dd.mm.yyyy HH:MM:SS"))

    dataset = createfiles(outdir,outstr,hcat(pts_x,pts_y),st,fileio,loc_time_agg)
    pts_lat, pts_lon = get_latlon(pts_x,pts_y,st.epsg_code)
    solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = st.tstep, radius = canrad.radius, time_zone = st.time_zone)

    @unpack radius, g_rad = canrad
    g_rad[g_rad .> radius] .= NaN

    @unpack ring_radius, ring_tht, surf_area_p, surf_area_h, relevant_pix = canrad
    for rix = 1:1:size(ring_radius,1)-1
        relevant_pix[:,rix] = (ring_radius[rix] .< vec(g_rad) .< ring_radius[rix+1])
        surf_area_p[rix] = sin(ring_tht[rix+1]/360*2*pi)^2 - sin(ring_tht[rix]/360*2*pi)^2
        surf_area_h[rix] = 2*pi*(cos(ring_tht[rix]/360*2*pi) - cos(ring_tht[rix+1]/360*2*pi))/2/pi
    end

    @simd for crx = 1:size(pts_x,1)

        sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],st.time_zone)
        @unpack trans_for = solar

        if (st.phenology in ("leafon","both")) && (pts_m[crx] .== 1)

            mat2ev = Int64.(shi_ds["SHI_leafon"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_leafon"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_leafon"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["tvt_leafon"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if st.calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,st.calc_swr)
                dataset["swr_total_leafon"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_leafon"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if (st.phenology in ("leafoff","both"))  && (pts_m[crx] .== 1)

            mat2ev = Int64.(shi_ds["SHI_leafoff"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_leafoff"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_leafoff"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["tvt_leafoff"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if st.calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,st.calc_swr)
                dataset["swr_total_leafoff"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_leafoff"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if st.calc_terrain

            mat2ev = Int64.(shi_ds["SHI_terrain"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_terrain"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_terrain"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["tvt_terrain"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if st.calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,st.calc_swr)
                dataset["swr_total_terrain"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_terrain"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if (st.forest_type == "evergreen") && (pts_m[crx] .== 1)

            mat2ev = Int64.(shi_ds["SHI_evergreen"][:,:,crx])

            svf_p, svf_h = calc_svf(canrad,mat2ev)
            dataset["svf_planar_evergreen"][crx] = Int8(round(svf_p*100))
            dataset["svf_hemi_evergreen"][crx]   = Int8(round(svf_h*100))

            fill!(trans_for,0)
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            dataset["tvt_evergreen"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

            if st.calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,st.calc_swr)
                dataset["swr_total_evergreen"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_evergreen"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
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