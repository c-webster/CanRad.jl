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

    outdir, outstr, crxstart, append_file, percentdone = organise_outf(taskID,exdir,batch,size(pts_x,1))
    crxstart = 1; append_file = false; percentdone = 0;     # force restart the tile
    global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crxstart-1)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    if calc_trans
        loc_time     = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file,loc_time_agg,time_zone)
        pts_lat, pts_lon = calc_latlon(pts_x,pts_y,coor_system)
        solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = tstep, radius = canrad.radius, time_zone = time_zone)
    else    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file)
    end

    if save_images 
        images = create_exmat(outdir,outstr,pts,canrad.mat2ev,append_file)
    end

    if save_hlm
        hlm = create_exhlm(outdir,outstr,pts,ter2rad)
    end

    ################################################################################
    # > Import surface data

    if terrain_highres

        limits_highres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,dtm_peri)
        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits_highres,true, false);
        
        if !isempty(dtm_x) || !isnan(sum(dtm_z))
            pts_e = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y,limits_highres,10.0,Vector{Float64}(undef,size(pts_x)))
            rbins_dtm = collect(2*dtm_cellsize:sqrt(2).*dtm_cellsize:dtm_peri)
            rows = findall(isnan,dtm_z); deleteat!(dtm_x,rows); deleteat!(dtm_y,rows); deleteat!(dtm_z,rows)
        # silly workaround to dealing with terrain data that doesn't have the same boundaries as grids.
        # this workaround loads the data including nans so that the above doesn't fail, and a horizon line
        #   can be calculated where high-res data exists. however, it will fill pts_e with nans, so this is 
        #   fixed below where the low-res terrain information is used. 
        # This fix is intended for the swissalti3d data which is clipped to the swiss border. 
        end

    end

    if terrain_lowres

        limits_lowres = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,terrain_peri)
        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits_lowres,true,true)
        pts_e_dem = findelev!(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,limits_lowres,50,Vector{Float64}(undef,size(pts_x)))
        rbins_dem = collect(2*dem_cellsize:sqrt(2).*dem_cellsize:terrain_peri)

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
        progtextinit = "0. Pre-calc took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
        if ispath(progdir)
            rm(progdir,recursive=true); mkdir(progdir)
        else
            mkdir(progdir)
        end
        writedlm(joinpath(progdir,progtextinit),NaN)            
    end

    dtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
    dem_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])

    ###############################################################################
    # > Loop through the points

    @simd for crx = crxstart:size(pts_x,1)

        if progress; start = time(); end

        # get the high-res local terrain
        if !isempty(dtm_x) && terrain_highres

            pt_dtm_x, pt_dtm_y, pt_dtm_z = getsurfdat(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],dtm_peri);
            pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,pt_dtm_x, pt_dtm_y, pt_dtm_z,pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)

            if save_hlm 
                if sum(isnan.(ter2rad.mintht)) .> 0
                    temp_mintht = replace(ter2rad.mintht,NaN=>missing)
                    dtm_mintht = Impute.interp(temp_mintht)[ter2rad.dx1:ter2rad.dx2-1]
                else
                    copy!(dtm_mintht,ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
                end
            end

        end

        # get the low-res regional terrain
        if terrain_lowres

            pt_dem_x, pt_dem_y, pt_dem_z = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e_dem[crx],terrain_peri);
            pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,pt_dem_x, pt_dem_y, pt_dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],"terrain",rbins_dem,image_height);

            if save_hlm
                if !isempty(dtm_x) && terrain_highres
                    if sum(isnan.(ter2rad.mintht)) .> 0
                        temp_mintht = replace(ter2rad.mintht,NaN=>missing)
                        dem_mintht = Impute.interp(temp_mintht)[ter2rad.dx1:ter2rad.dx2-1]
                    else
                        copy!(dem_mintht,ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
                    end
                    hlm["tht"][:,crx] = Int8.(round.(minimum(hcat(dtm_mintht,dem_mintht),dims=2)))
                else
                    hlm["tht"][:,crx] = Int8.(round.(copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])))
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

        # calculate Vf and transmissivity
        Vf_p, Vf_h = calculateVf(canrad,mat2ev)

        dataset["Vf_planar_terrain"][crx] = Int8(round(Vf_p*100));
        dataset["Vf_hemi_terrain"][crx]   = Int8(round(Vf_h*100));

        if calc_trans

            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],time_zone)

            @unpack trans_for = solar
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)

            dataset["Forest_Transmissivity_terrain"][:,crx] = Int8.(round.(vec(aggregate_data(solar,trans_for))).*100);

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,Vf_p,calc_swr)
                dataset["SWR_total_terrain"][:,crx]  = Int16.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["SWR_direct_terrain"][:,crx] = Int16.(round.(vec(aggregate_data(solar,swrdir))))
            end
        
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "Process took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        # save the progress
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*sprintf1.("%.$(0)f", Int(floor((crx / size(pts,1)) * 100)))*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end

    close(dataset)
    save_images && close(images)
    save_hlm && close(hlm)

    println("done with "*taskID)

end