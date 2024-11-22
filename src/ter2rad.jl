function ter2rad!(pts::Matrix{Float64},dat_in::Dict{String, String},par_in::Dict{String, Any},
    exdir::String,taskID="task")

    ################################################################################
    # Initialise

    # run compatability check then extract settings
    compatability_check!(par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    progress && (start = time())

    # separate the points to vectors
    pts_x, pts_y = [pts[:,ptdx] for ptdx in 1:size(pts,2)]

    canrad  = CANRAD()
    ter2rad = TER2RAD(pts_sz = size(pts_x,1))

    ################################################################################
    # > Get constants, organise the output and initiate progress reporting

    outdir, outstr  = organise_outf(taskID,exdir,batch)
    global outtext = "Processing "*cfmt.("%.$(0)f", 0)*"% ... "*string(0)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    if calc_trans
        loc_time     = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
    
        dataset = createfiles_terrain(outdir,outstr,pts,calc_trans,calc_swr,loc_time_agg,time_zone)
        pts_lat, pts_lon = calc_latlon(pts_x,pts_y,coor_system)
        solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = tstep, radius = canrad.radius, time_zone = time_zone)
    else    
        dataset = createfiles_terrain(outdir,outstr,pts,calc_trans,calc_swr)
    end

    save_images && (images = create_exmat_terrain(outdir,outstr,pts,canrad.mat2ev))

    save_horizon && (hlm = create_exhlm(outdir,outstr,pts,ter2rad))

    ################################################################################
    # > Import surface data

    if terrain_highres

        limits_highres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,highres_peri)
        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits_highres,true, false);
        
        if !isempty(dtm_x) || !isnan(sum(dtm_z))
            pts_e = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y,limits_highres,10.0,Vector{Float64}(undef,size(pts_x)))
            rbins_dtm = collect(2*dtm_cellsize:sqrt(2).*dtm_cellsize:highres_peri)
            rows = findall(isnan,dtm_z); deleteat!(dtm_x,rows); deleteat!(dtm_y,rows); deleteat!(dtm_z,rows)
        # silly workaround to dealing with terrain data that doesn't have the same boundaries as grids.
        # this workaround loads the data including nans so that the above doesn't fail, and a horizon line
        #   can be calculated where high-res data exists. however, it will fill pts_e with nans, so this is 
        #   fixed below where the low-res terrain information is used. 
        # This fix is intended for the swissalti3d data which is clipped to the swiss border. 
        end

    end

    if terrain_lowres

        limits_lowres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,lowres_peri)
        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits_lowres,true,true)
        pts_e_dem = findelev!(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,limits_lowres,50,Vector{Float64}(undef,size(pts_x)))
        rbins_dem = collect(2*dem_cellsize:sqrt(2).*dem_cellsize:lowres_peri)

    end

    if isempty(dtm_x) || isnan(sum(pts_e))
        pts_e = copy(pts_e_dem)
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
    # > SWR calculations

    if calc_swr > 0
        radiation = RADIATION(loc_time = loc_time)
        if calc_swr == 2
            swrdat   = readdlm(swrf)
            radiation.swr_open = float.(swrdat[:,3])
        end
    end

    ###############################################################################
    # start detailed progress reporting
    if progress
        elapsed = time() - start
        progdir = joinpath(outdir,"ProgressLastPoint")
        progtextinit = "0. Pre-calc took "*cfmt.("%.$(2)f", elapsed)*" seconds"
        if ispath(progdir)
            rm(progdir,recursive=true); mkdir(progdir)
        else
            mkdir(progdir)
        end
        writedlm(joinpath(progdir,progtextinit),NaN)            
    end

    ###############################################################################
    # > Loop through the points

    @simd for crx = 1:size(pts_x,1)

        progress && (start = time())

        # get the high-res local terrain
        if !isempty(dtm_x) && terrain_highres

            pt_dtm_x, pt_dtm_y, pt_dtm_z = getsurfdat(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],highres_peri);
            pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,pt_dtm_x, pt_dtm_y, pt_dtm_z,pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)

            if save_horizon 
                dtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
            end

        end

        # get the low-res regional terrain
        if terrain_lowres

            pt_dem_x, pt_dem_y, pt_dem_z = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e_dem[crx],lowres_peri);
            pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,pt_dem_x, pt_dem_y, pt_dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],"terrain",rbins_dem,image_height);

            if save_horizon
                if !isempty(dtm_x) && terrain_highres
                    dem_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
                    hlm["tht"][:,crx] = Int.(round.(minimum(hcat(dtm_mintht,dem_mintht),dims=2) .*100))
                else
                    hlm["tht"][:,crx] = Int.(round.(copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1]) .* 100))
                end
            end
        
        end
        
        # combine the datasets and occupy the image matrix
        fill!(mat2ev,1);
        if !isempty(dtm_x) && terrain_highres
            prepterdat!(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
            fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),13,mat2ev);
        else
            prepterdat!(pt_dem_x,pt_dem_y);
            fillmat!(canrad,kdtree,hcat(pt_dem_x,pt_dem_y),13,mat2ev);
        end


        # create the image matrix
        mat2ev[outside_img] .= 1;

        save_images && (images["SHI_terrain"][:,:,crx] = mat2ev;)

        # calculate svf and transmissivity
        svf_p, svf_h = calc_svf(canrad,mat2ev)

        dataset["svf_planar_t"][crx] = Int(round(svf_p*100));
        dataset["svf_hemi_t"][crx]   = Int(round(svf_h*100));

        if calc_trans

            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],time_zone)

            @unpack trans_for = solar
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)

            dataset["trans_t"][:,crx] = Int.(round.((vec(aggregate_data(solar,trans_for))).*100));

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,calc_swr)
                dataset["swr_total_t"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct_t"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end
        
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "Process took "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        # save the progress
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*cfmt.("%.$(0)f", Int(floor((crx / size(pts,1)) * 100)))*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end

    close(dataset)
    save_images && close(images)
    save_horizon && close(hlm)

    (save_images && make_pngs) && make_SHIs(outdir,"none","none",true)

    println("done with "*taskID)

end