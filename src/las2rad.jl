function las2rad!(pts::Matrix{Float64},dat_in::Dict{String, String},par_in::Dict{String, Any},
    exdir::String,taskID="task")

    ################################################################################
    # > Initialise

    # run compatability check then extract settings
    compatability_check!(par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    progress && (start = time())

    # separate the points to vectors
    pts_x, pts_y =[pts[:,x] for x in 1:size(pts,2)]

    # get model number (1 = canopy, 2 = terrain only)
    size(pts,2) > 2 ? pts_m = pts[:,3] : pts_m = ones(size(pts_x))

    canrad  = CANRAD()
    las2rad = LAS2RAD(point_size=point_size,forest_peri=forest_peri)

    ################################################################################
    # > Get constants, organise the output and initiate progress reporting

    outdir, outstr  = organise_outf(taskID,exdir,batch)
    crxstart = 1; percentdone = 0; # initialise the tile
    global outtext = "Processing "*cfmt.("%.$(0)f", percentdone)*"% ... "*string(crxstart-1)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    if calc_trans
        loc_time     = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,loc_time_agg,time_zone)
        pts_lat, pts_lon = get_latlon(pts_x,pts_y,epsg_code)
        solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = tstep, radius = canrad.radius, time_zone = time_zone)
    else    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr)
    end

    save_images && (images = create_exmat(outdir,outstr,pts,canrad.mat2ev))

    ################################################################################
    # > Import surface data

    @unpack limits_canopy = canrad
    limits_canopy = getlimits!(limits_canopy,pts_x,pts_y,forest_peri)

    #  "canopy", "ground+canopy", "all"
    if keep_points == "canopy"
        keep_ground = false
        keep_all = false
    elseif keep_points == "ground+canopy"
        keep_ground = true
        keep_all = false
    elseif keep_points == "all"
        keep_ground = false
        keep_all = true
    end

    dsm_x, dsm_y, dsm_z = readlas(dsmf,limits_canopy,keep_ground,keep_all)

    ################################################################################
    # > Import and prepare terrain data

    if terrainmask_precalc
        terrain_mask = getterrainmask(canrad,terf,pts_x,pts_y)
    elseif horizon_line
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
        hlm_tht = load_hlm(ter2rad,hlmf,pts_x,pts_y)
    end

    if terrain_highres || terrain_lowres || horizon_line
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
    end

    if terrain_highres

        @unpack limits_highres, pts_e = ter2rad
        getlimits!(limits_highres,pts_x,pts_y,highres_peri)
        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits_highres,true, true)
        rbins_dtm = collect(2*dtm_cellsize:sqrt(2).*dtm_cellsize:highres_peri)

        pts_e = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y,limits_highres,10.0,pts_e)

    end

    # Get the low-res terrain data
    if terrain_lowres && !horizon_line

        @unpack limits_lowres, pts_e_dem = ter2rad
        getlimits!(limits_lowres,pts_x,pts_y,lowres_peri)
        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits_lowres,true,true)
        rbins_dem = collect(2*dem_cellsize:sqrt(2).*dem_cellsize:lowres_peri)
        findelev!(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,limits_lowres,dem_cellsize*3,pts_e_dem)

    end

    # determine ground elevation of points
    if size(pts,2) > 2
        pts_e = pts[:,3]
    elseif terrain_highres
        pts_e     = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y,limits_highres,10,Vector{Float64}(undef,size(pts_x)))
    elseif terrain_lowres
        pts_e     = pts_e_dem
    else # if there's no terrain input, just assume the ground is flat
        pts_e     = zeros(size(pts_x))
    end

    # determine true elevatation of las points if lidar is normalised
        # (if the absolute difference of mode(lidar) and mode(terrain) is greater
        #   than a realistic canopy height (60 m)
    if abs(mode(dsm_z) - mode(dtm_z)) > 60 && terrain_highres
        dsm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),dsm_x,dsm_y)
    end

    # > Calculate slope and aspect from dtm data
    if tilt # currently disabled in Preparatory_Functions.jl
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x))
    end

    # load the buildings
    if buildings
        bhm_x, bhm_y, bhm_z, bhm_cellsize = read_griddata_window(bhmf,limits_canopy,true,true)
        rbins_bhm = collect(2*bhm_cellsize:sqrt(2).*bhm_cellsize:forest_peri)
        if !isempty(bhm_x)
            build = true
            if !(abs(mode(bhm_z) - mode(dtm_z)) < 60) # check if data is normalised, if yes, fit to terrain
                bhm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),bhm_x,bhm_y)
            end
    	else; build = false; end
    else; build = false; end

    ###############################################################################
    # > Generate extra canopy elements

    limits_trees =  getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,forest_peri)

    # load the dbh
    if trunks
        dbh_x, dbh_y, dbh_z, dbh_r, lastc = loaddbh(dbhf,limits_trees)
        if !isempty(dbh_x)
            dbh_e = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),dbh_x,dbh_y,limits_trees,10.0,Vector{Float64}(undef,size(dbh_x))            )
            tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,30,0.1,dbh_e)
            include_trunks = true
        else; include_trunks = false; end
    else; include_trunks = false;
    end

    # load the ltc
    if branches
        if extension(ltcf) == ".txt"
            ltc = loadltc_txt(ltcf,limits_trees,0)
        elseif extension(ltcf) == ".laz"
             ltc = loadltc_laz(ltcf,limits_trees,dbh_x,dbh_y,lastc,season)  
            if abs(mode(ltc.z) - mode(dtm_z)) < 60 # if the data's not normalised, it needs to be normalised)
                ltc.z .-= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),ltc.x,ltc.y)
            end
            # remove the points < 1m from the ground
            ltc = ltc[setdiff(1:end, findall(ltc.z.<1)), :]
        end
        bsm_x, bsm_y, bsm_z, _  = make_branches(ltc.x,ltc.y,ltc.z,
                    ltc.tx,ltc.ty,ltc.hd,ltc.ang,ltc.cls,branch_spacing,season)

        bsm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),bsm_x,bsm_y)
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

    calc_trans && (@unpack trans_for = solar)

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
    # > get a couple more constants

    @unpack rbins, knum, knum_t = las2rad

    ###############################################################################
    # > Loop through the points

    @simd for crx = 1:size(pts_x,1)

        progress && (start = time())

        #### transfer point clouds to polar coordinates
        if branches
            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(copy(dsm_x),copy(dsm_y),copy(dsm_z),copy(bsm_x),copy(bsm_y),copy(bsm_z),pts_x[crx],pts_y[crx],pts_e[crx],forest_peri)
        else
            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(copy(dsm_x),copy(dsm_y),copy(dsm_z),pts_x[crx],pts_y[crx],pts_e[crx],forest_peri);
        end
        pt_dsm_x, pt_dsm_y, pt_dsm_z = pcd2pol2cart!(pt_dsm_x,pt_dsm_y,pt_dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],"surface",image_height)

        prepsurfdat!(pt_dsm_x, pt_dsm_y, pt_dsm_z) 

        if include_trunks
            pt_tsm_x, pt_tsm_y, pt_tsm_z = getsurfdat(copy(tsm_x),copy(tsm_y),copy(tsm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(forest_peri))
            tidx = findall(dist(dbh_x,dbh_y,pts_x[crx],pts_y[crx]) .< 5)
            if size(tidx,1) > 0
                hdt  = dist(dbh_x[tidx],dbh_y[tidx],pts_x[crx],pts_y[crx])
                npt  = fill(NaN,(size(tidx)))
                hint = fill(NaN,(size(tidx)))
                        for tixt = 1:1:size(tidx,1)
                    if hdt[tixt] < 1
                        npt[tixt] = Int.(150); hint[tixt] = 0.005
                    else
                        npt[tixt] = Int.(100); hint[tixt] = 0.01
                    end
                end
                tsm_tmp = calculate_trunks(dbh_x[tidx],dbh_y[tidx],dbh_z[tidx],dbh_r[tidx],npt,hint,dbh_e[tidx])
                pt_tsm_x, pt_tsm_y, pt_tsm_z = pcd2pol2cart!(append!(pt_tsm_x,tsm_tmp[1]),append!(pt_tsm_y,tsm_tmp[2]),append!(pt_tsm_z,tsm_tmp[3]),
                                                        pts_x[crx],pts_y[crx],pts_e[crx],"trunks",image_height)

            else
                pcd2pol2cart!(pt_tsm_x,pt_tsm_y,pt_tsm_z,pts_x[crx],pts_y[crx],pts_e[crx],"trunks",image_height)
            end

            if season == "winter"
                prepsurfdat!(pt_tsm_x,pt_tsm_y,pt_tsm_z)
            else
                prepsurfdat!(append!(pt_dsm_x,pt_tsm_x),append!(pt_dsm_y,pt_tsm_y),append!(pt_dsm_z,pt_tsm_z))
            end
        end

        if !terrainmask_precalc

            if terrain_highres
                # get the high-res local terrain
                pt_dtm_x, pt_dtm_y, pt_dtm_z = getsurfdat(copy(dtm_x),copy(dtm_y),copy(dtm_z),
                                                pts_x[crx],pts_y[crx],pts_e[crx],highres_peri);
                pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,pt_dtm_x, pt_dtm_y, pt_dtm_z,pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)
        
            end

            if terrain_lowres 
                # get the low-res regional terrain
                pt_dem_x, pt_dem_y, pt_dem_z = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e[crx],lowres_peri);
                pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,pt_dem_x, pt_dem_y, pt_dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],"terrain",rbins_dem,image_height);
        
            end

            if horizon_line && !oshd_flag
                pt_dtm_x,pt_dtm_y = hlm2cart(ter2rad,hlm_tht[:,crx])
            elseif terrain_highres && (terrain_lowres || (horizon_line && oshd_flag))
                prepterdat!(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
            elseif terrain_highres
                prepterdat!(pt_dtm_x,pt_dtm_y)
            elseif (terrain_lowres && !terrain_highres)
                pt_dtm_x, pt_dtm_y = prepterdat(pt_dem_x,pt_dem_y)
            end

        end

        if build

            pt_bhm_x, pt_bhm_y, pt_bhm_z = getsurfdat(copy(bhm_x),copy(bhm_y),copy(bhm_z),pts_x[crx],pts_y[crx],pts_e[crx],300);    
            pt_bhm_x, pt_bhm_y =  pcd2pol2cart!(ter2rad,pt_bhm_x,pt_bhm_y,pt_bhm_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                "buildings",rbins_bhm,image_height)
            if !terrainmask_precalc
                prepterdat!(append!(pt_dtm_x,pt_bhm_x),append!(pt_dtm_y,pt_bhm_y));
            else
                pt_dtm_x, pt_dtm_y = prepterdat(pt_bhm_x,pt_bhm_y)
            end

        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "1. Transferring to polar took "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        ######################
        ### Classify image
        progress && (start = time())

        terrainmask_precalc ? (mat2ev .= terrain_mask[:,:,crx]) : (fill!(mat2ev,1))

        # occupy matrix with surface points
        for zdx = 1:1:size(rbins,1)-1
            ridx = findall(rbins[zdx] .<= pt_dsm_z .< rbins[zdx+1])
            fillmat!(canrad,kdtree,hcat(pt_dsm_x[ridx],pt_dsm_y[ridx]),knum[zdx],mat2ev)
            
            if season == "winter" && trunks
                tridx = findall(rbins[zdx] .<= pt_tsm_z .< rbins[zdx+1])
                fillmat!(canrad,kdtree,hcat(pt_tsm_x[tridx],pt_tsm_y[tridx]),knum_t[zdx],mat2ev)
            end
        end

        # add terrain
        !terrainmask_precalc && (fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),10,mat2ev))

        (season == "winter" && trunks) && (fillmat!(canrad,kdtree,hcat(pt_tsm_x,pt_tsm_y),10,mat2ev))

        mat2ev[outside_img] .= 1

        save_images && (images["SHI"][:,:,crx] = mat2ev)

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext2)); catch; end; end
            global progtext2 = "2. Classifying image took "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext2),NaN)
        end

        ########################
        ### Perform calculations
        progress && (start = time())

        ##### Calculate svf and save
        svf_p, svf_h = calc_svf(canrad,mat2ev)

        dataset["svf_planar"][crx] = Int(round(svf_p*100))
        dataset["svf_hemi"][crx]   = Int(round(svf_h*100))

        ##### CalcSWR
        if calc_trans
            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],time_zone)

            @unpack trans_for = solar
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)

            dataset["for_trans"][:,crx] = Int.(round.((vec(aggregate_data(solar,trans_for))).*100));

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,calc_swr)
                dataset["swr_total"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext3)); catch; end; end
            global progtext3 = "3. Calculations and export took "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext3),NaN)
        end

        # save the progress
        percentdone = Int(floor((crx / size(pts,1)) * 100))
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*cfmt.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end # end crx

    close(dataset)
    save_images && close(images)

    (save_images && make_pngs) && make_SHIs(outdir)

end
