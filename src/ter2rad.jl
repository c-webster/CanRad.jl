function ter2rad!(pts::Matrix{Float64},dat_in::Dict{String, String},par_in::Dict{String, Any},
    exdir::String,taskID="task")

    ################################################################################
    # Initialise

    # run compatability check then extract settings
    update_deprecated_settings!(par_in,dat_in)
    terrain_defaults!(par_in) # sets forest settings to "none"

    fp = FILEPATHS(; (Symbol.(keys(dat_in)) .=> values(dat_in))... )
    st = SETTINGS(; (Symbol.(keys(par_in)) .=> values(par_in))... )

    check_conflicts(st,fp,"t2r")

    st.step_progress && (start = time())

    # separate the points to vectors
    pts_x, pts_y = [pts[:,ptdx] for ptdx in 1:size(pts,2)]

    canrad  = CANRAD()
    fileio  = FILEIO(time_zone=st.time_zone,calc_swr=st.calc_swr)
    ter2rad = TER2RAD(pts_sz = size(pts_x,1))

    ################################################################################
    # > Get constants, organise the output and initiate progress reporting

    outdir, outstr  = organise_outf(taskID,exdir,st.batch)
    global outtext = "Processing "*cfmt.("%.$(0)f", 0)*"% ... "*string(0)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    ################################################################################
    # > Import surface data

    if st.terrain_highres

        limits_highres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,st.highres_peri)
        hrdtm_x, hrdtm_y, hrdtm_z, hrdtm_cellsize = read_griddata_window(fp.hrdtmf,limits_highres,true, false);
        
        if !isempty(hrdtm_x) || !isnan(sum(hrdtm_z))
            pts_e = findelev!(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),pts_x,pts_y,limits_highres,10.0,Vector{Float64}(undef,size(pts_x)))
            rbins_hrdtm = collect(2*hrdtm_cellsize:sqrt(2).*hrdtm_cellsize:st.highres_peri)
            rows = findall(isnan,hrdtm_z); deleteat!(hrdtm_x,rows); deleteat!(hrdtm_y,rows); deleteat!(hrdtm_z,rows)
        # silly workaround to dealing with terrain data that doesn't have the same boundaries as grids.
        # this workaround loads the data including nans so that the above doesn't fail, and a horizon line
        #   can be calculated where high-res data exists. however, it will fill pts_e with nans, so this is 
        #   fixed below where the low-res terrain information is used. 
        # This fix is intended for the swissalti3d data which is clipped to the swiss border. 
        end

    else

        hrdtm_x = [] # define an empty variable to run with terrain_highres only (fix in the future to something less blunt)

    end

    if st.terrain_lowres

        limits_lowres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,st.lowres_peri)
        lrdtm_x, lrdtm_y, lrdtm_z, lrdtm_cellsize = read_griddata_window(fp.lrdtmf,limits_lowres,true,true)
        pts_e_lrdtm = findelev!(copy(lrdtm_x),copy(lrdtm_y),copy(lrdtm_z),pts_x,pts_y,limits_lowres,50,Vector{Float64}(undef,size(pts_x)))
        rbins_lrdtm = collect(2*lrdtm_cellsize:sqrt(2).*lrdtm_cellsize:st.lowres_peri)

    end

    if isempty(hrdtm_x) || isnan(sum(pts_e))
        pts_e = copy(pts_e_lrdtm)
    end

    ###############################################################################
    # > Image matrix  preparation

    # create an empty image matrix
    @unpack radius, mat2ev, g_rad, g_coorcrt = canrad
    g_rad[g_rad .> radius] .= NaN
    outside_img = isnan.(g_rad)
    g_coorcrt .= ((g_coorcrt .- radius) ./ radius) .* 90

    # make g_coorcrt a KDtree for easy look up
    kdtree = scipyspat.cKDTree(g_coorcrt)

    @unpack ring_radius, ring_tht, surf_area_p, surf_area_h, relevant_pix = canrad
    for rix = 1:1:size(ring_radius,1)-1
        relevant_pix[:,rix] = (ring_radius[rix] .< vec(g_rad) .< ring_radius[rix+1])
        surf_area_p[rix] = sin(ring_tht[rix+1]/360*2*pi)^2 - sin(ring_tht[rix]/360*2*pi)^2
        surf_area_h[rix] = 2*pi*(cos(ring_tht[rix]/360*2*pi) - cos(ring_tht[rix+1]/360*2*pi))/2/pi
    end

    ###############################################################################
    # > Initialise variables for transmissivity and swr calculations

    if st.calc_trans
        loc_time     = collect(Dates.DateTime(st.t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(st.t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(st.t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(st.tstep):Dates.DateTime(st.t2,"dd.mm.yyyy HH:MM:SS"))
        solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = st.tstep, radius = canrad.radius, time_zone = st.time_zone)

        @unpack trans_for = solar

    end

    if st.calc_swr > 0
        radiation = RADIATION(loc_time = loc_time)
        if st.calc_swr == 2
            swrdat   = readdlm(fp.swrf)
            radiation.swr_open = float.(swrdat[:,3])
        end
    end

    ###############################################################################
    # > Create the output file structure

    if st.calc_trans
        dataset = createfiles(outdir,outstr,pts,st,fileio,loc_time_agg)
        pts_lat, pts_lon = get_latlon(pts_x,pts_y,st.epsg_code)
    else    
        dataset = createfiles(outdir,outstr,pts,st,fileio)
    end

    st.save_images && (images = create_exmat(outdir,outstr,pts,canrad.mat2ev,st))

    st.save_horizon && (hlm = create_exhlm(outdir,outstr,pts,ter2rad))

    ###############################################################################
    # start detailed progress reporting
    if st.step_progress
        progdir = joinpath(outdir,"ProgressLastPoint")
        progtextinit = "0. Pre-calc took "*cfmt.("%.$(2)f",  time() - start)*" seconds"
        ispath(progdir) ? (rm(progdir,recursive=true); mkdir(progdir)) : mkdir(progdir)
        writedlm(joinpath(progdir,progtextinit),NaN)            
    end

    ###############################################################################
    # > Loop through the points

    @simd for crx = 1:size(pts_x,1)

        st.step_progress && (start = time())

        # get the high-res local terrain
        if !isempty(hrdtm_x) && st.terrain_highres

            pt_hrdtm_x, pt_hrdtm_y, pt_hrdtm_z = getsurfdat(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),pts_x[crx],pts_y[crx],pts_e[crx],st.highres_peri);
            pt_hrdtm_x, pt_hrdtm_y = pcd2pol2cart!(ter2rad,pt_hrdtm_x, pt_hrdtm_y, pt_hrdtm_z,pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_hrdtm,st.image_height)

            if st.save_horizon 
                hrdtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
            end

        end

        # get the low-res regional terrain
        if st.terrain_lowres

            pt_lrdtm_x, pt_lrdtm_y, pt_lrdtm_z = getsurfdat(copy(lrdtm_x),copy(lrdtm_y),copy(lrdtm_z),pts_x[crx],pts_y[crx],pts_e_lrdtm[crx],st.lowres_peri);
            pt_lrdtm_x, pt_lrdtm_y = pcd2pol2cart!(ter2rad,pt_lrdtm_x, pt_lrdtm_y, pt_lrdtm_z,pts_x[crx],pts_y[crx],pts_e_lrdtm[crx],"terrain",rbins_lrdtm,st.image_height);

            if st.save_horizon
                if !isempty(hrdtm_x) && st.terrain_highres
                    lrdtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
                    hlm["tht"][:,crx] = Int.(round.(minimum(hcat(hrdtm_mintht,lrdtm_mintht),dims=2) .*100))
                else
                    hlm["tht"][:,crx] = Int.(round.(copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1]) .* 100))
                end
            end
        
        end
        
        # combine the datasets and occupy the image matrix
        fill!(mat2ev,1);
        if !isempty(hrdtm_x) && st.terrain_highres
            prepterdat!(append!(pt_hrdtm_x,pt_lrdtm_x),append!(pt_hrdtm_y,pt_lrdtm_y));
            fillmat!(canrad,kdtree,hcat(pt_hrdtm_x,pt_hrdtm_y),13,mat2ev);
        else
            prepterdat!(pt_lrdtm_x,pt_lrdtm_y);
            fillmat!(canrad,kdtree,hcat(pt_lrdtm_x,pt_lrdtm_y),13,mat2ev);
        end


        # create the image matrix
        mat2ev[outside_img] .= 1;

        st.save_images && (images["SHI_terrain"][:,:,crx] = mat2ev;)

        # calculate svf and transmissivity
        svf_p, svf_h = calc_svf(canrad,mat2ev)

        dataset["svf_planar_terrain"][crx] = Int(round(svf_p*100));
        dataset["svf_hemi_terrain"][crx]   = Int(round(svf_h*100));

        if st.calc_trans

            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],st.time_zone)

            @unpack trans_for = solar
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)

            dataset["tvt_terrain"][:,crx] = Int.(round.((vec(aggregate_data(solar,trans_for))).*100));

            if st.calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,st.calc_swr)
                dataset["swr_total_terrain"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_terrain"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end
        
        end

        if st.step_progress
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "Process took "*cfmt.("%.$(2)f", time() - start)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        # save the progress
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*cfmt.("%.$(0)f", Int(floor((crx / size(pts,1)) * 100)))*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end

    close(dataset)
    st.save_images && close(images)
    st.save_horizon && close(hlm)

    (st.save_images && st.make_pngs) && make_SHIs(outdir)
    
    # Save settings and data information to text file if not running batch
    if !st.batch 
        write_metadata(outdir,outstr,st,fp)
    end

end