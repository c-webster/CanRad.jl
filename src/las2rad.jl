function las2rad!(pts::Matrix{Float64},dat_in::Dict{String, String},par_in::Dict{String, Any},
    exdir::String,taskID="task")

    ################################################################################
    # > Initialise

    # run compatability check, extract settings and check for basic conflicts
    update_deprecated_settings!(par_in,dat_in)

    fp = FILEPATHS(; (Symbol.(keys(dat_in)) .=> values(dat_in))... )
    st = SETTINGS(; (Symbol.(keys(par_in)) .=> values(par_in))... )

    check_conflicts(st,fp,"l2r")

    # load the required constants
    canrad  = CANRAD()
    fileio  = FILEIO(time_zone=st.time_zone,calc_swr=st.calc_swr)
    las2rad = LAS2RAD(point_size=st.point_size,forest_peri=st.forest_peri)

    ################################################################################
    # > Get constants, organise the output and initiate progress reporting

    st.step_progress && (start = time())

    outdir, outstr  = organise_outf(taskID,exdir,st.batch)
    global outtext = "Processing "*cfmt.("%.$(0)f", 0)*"% ... "*string(0)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    ################################################################################
    # > Organise the point vectors depending on input into the function

    # separate the points to vectors
    pts_x, pts_y =[pts[:,x] for x in 1:size(pts,2)]

    # get model number (1 = canopy, 2 = terrain only)
    size(pts,2) > 2 ? pts_m = pts[:,3] : pts_m = ones(size(pts_x))

    ################################################################################
    # > Import surface data

    @unpack limits_canopy = canrad
    limits_canopy = getlimits!(limits_canopy,pts_x,pts_y,st.forest_peri)

    st.keep_points in ("canopy","all") ? keep_ground = false : keep_ground = true
    st.keep_points in ("canopy","ground+canopy") ? keep_all = false : keep_all = true

    dsm_x, dsm_y, dsm_z = readlas(fp.lasf,limits_canopy,keep_ground,keep_all)

    ################################################################################
    # > Import and prepare terrain data

    if st.terrain_highres || st.terrain_lowres || st.hlm_precalc || st.terrainmask_precalc
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
        flat_terrain_flag = false
        if st.hlm_precalc
            if st.special_implementation == "oshd"
            oshd_flag = true
            pt_lrdtm_x, pt_lrdtm_y = load_hlm_oshd(fp.hlmf)
            else
                oshd_flag = false
                hlm_tht = load_hlm(ter2rad,fp.hlmf,pts_x,pts_y)
            end
        elseif st.terrainmask_precalc
            terrain_mask = getterrainmask(canrad,fp.termf,pts_x,pts_y)
        end
    else
        flat_terrain_flag = true
    end

    if st.terrain_highres #&& !(st.hlm_precalc || st.terrainmask_precalc)

        @unpack limits_highres, pts_e = ter2rad
        getlimits!(limits_highres,pts_x,pts_y,st.highres_peri)
        hrdtm_x, hrdtm_y, hrdtm_z, hrdtm_cellsize = read_griddata_window(fp.hrdtmf,limits_highres,true, true)
        rbins_hrdtm = collect(2*hrdtm_cellsize:sqrt(2).*hrdtm_cellsize:st.highres_peri)

    end

    # Get the low-res terrain data
    if st.terrain_lowres && !(st.hlm_precalc || st.terrainmask_precalc)

        @unpack limits_lowres, pts_e_lrdtm = ter2rad
        getlimits!(limits_lowres,pts_x,pts_y,st.lowres_peri)
        lrdtm_x, lrdtm_y, lrdtm_z, lrdtm_cellsize = read_griddata_window(fp.lrdtmf,limits_lowres,true,true)
        rbins_lrdtm = collect(lrdtm_cellsize:sqrt(2).*lrdtm_cellsize:st.lowres_peri)
        findelev!(copy(lrdtm_x),copy(lrdtm_y),copy(lrdtm_z),pts_x,pts_y,limits_lowres,lrdtm_cellsize*3,pts_e_lrdtm)

    end

    # determine ground elevation of points
    if size(pts,2) > 2
        pts_e = pts[:,3]
    elseif st.terrain_highres
        pts_e = findelev!(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),pts_x,pts_y,limits_highres,10.0,pts_e)
    elseif st.terrain_lowres
        pts_e     = pts_e_lrdtm
    else # if there's no terrain input, just assume the ground is flat
        pts_e     = zeros(size(pts_x))
    end

    # determine true elevatation of las points if lidar is normalised
        # (if the absolute difference of mode(lidar) and mode(terrain) is greater
        #   than a realistic canopy height (60 m)
    if (st.terrain_highres || st.terrain_highres) && abs(mode(dsm_z) - mode(hrdtm_z)) > 60
        dsm_z .+= findelev(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),dsm_x,dsm_y)
    end

    # > Calculate slope and aspect from dtm data
    if st.tilt # currently disabled in Preparatory_Functions.jl
        pts_slp = findelev(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_s),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x))
    end

    # load the buildings
    if st.buildings
        bhm_x, bhm_y, bhm_z, bhm_cellsize = read_griddata_window(fp.bhmf,limits_canopy,true,true)
        rbins_bhm = collect(2*bhm_cellsize:sqrt(2).*bhm_cellsize:st.forest_peri)
        if !isempty(bhm_x)
            buildings_flag = true
            if !(abs(mode(bhm_z) - mode(hrdtm_z)) < 60) # check if data is normalised, if yes, fit to terrain
                bhm_z .+= findelev(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),bhm_x,bhm_y)
            end
    	else; buildings_flag = false; end
    else; buildings_flag = false; end

    ###############################################################################
    # > Generate extra canopy elements

    limits_trees =  getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,st.forest_peri)

    # load the dbh
    if st.trunks
        dbh_x, dbh_y, dbh_z, dbh_r, lastc = loaddbh(fp.dbhf,limits_trees)
        if !isempty(dbh_x)
            dbh_e = findelev!(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),dbh_x,dbh_y,limits_trees,10.0,Vector{Float64}(undef,size(dbh_x))            )
            tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,30,0.1,dbh_e)
            include_trunks = true
        else; include_trunks = false; end
    else; include_trunks = false;
    end

    # load the ltc
    if st.branches
        if splitext(fp.ltcf)[2] == ".txt"
            ltc = loadltc_txt(fp.ltcf,limits_trees,0)
        elseif splitext(fp.ltcf)[2] == ".laz"
             ltc = loadltc_laz(fp.ltcf,limits_trees,dbh_x,dbh_y,lastc,st.phenology)  
            if abs(mode(ltc.z) - mode(hrdtm_z)) < 60 # if the data's not normalised, it needs to be normalised)
                ltc.z .-= findelev(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),ltc.x,ltc.y)
            end
            # remove the points < 1m from the ground
            ltc = ltc[setdiff(1:end, findall(ltc.z.<1)), :]
        end
        bsm_x, bsm_y, bsm_z, _  = make_branches(ltc.x,ltc.y,ltc.z,
                    ltc.tx,ltc.ty,ltc.hd,ltc.ang,ltc.cls,st.branch_spacing,st.phenology)

        bsm_z .+= findelev(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),bsm_x,bsm_y)
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

    if st.calc_trans
        dataset = createfiles(outdir,outstr,pts,st,fileio,loc_time_agg)
        pts_lat, pts_lon = get_latlon(pts_x,pts_y,st.epsg_code)
    else    
        dataset = createfiles(outdir,outstr,pts,st,fileio)
    end

    st.save_images && (images = create_exmat(outdir,outstr,pts,canrad.mat2ev,st))

    envstrings = environments_flag(st)
    outsuf = "_$(envstrings[1])"

    ###############################################################################
    # start detailed progress reporting
    if st.step_progress
        progdir = joinpath(outdir,"ProgressLastPoint")
        progtextinit = "0. Pre-calc took "*cfmt.("%.$(2)f", time() - start)*" seconds"
        ispath(progdir) ? (rm(progdir,recursive=true); mkdir(progdir)) : mkdir(progdir)
        writedlm(joinpath(progdir,progtextinit),NaN)            
    end

    ###############################################################################
    # > get a couple more constants

    @unpack rbins, knum, knum_t = las2rad

    ###############################################################################
    # > Loop through the points

    @simd for crx = 1:size(pts_x,1)

        st.step_progress && (start = time())

        #### transfer point clouds to polar coordinates
        if st.branches
            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(copy(dsm_x),copy(dsm_y),copy(dsm_z),copy(bsm_x),copy(bsm_y),copy(bsm_z),pts_x[crx],
                                                pts_y[crx],pts_e[crx],st.forest_peri)
        else
            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(copy(dsm_x),copy(dsm_y),copy(dsm_z),pts_x[crx],pts_y[crx],pts_e[crx],st.forest_peri);
        end
        pt_dsm_x, pt_dsm_y, pt_dsm_z = pcd2pol2cart!(pt_dsm_x,pt_dsm_y,pt_dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],"surface",st.image_height)

        prepsurfdat!(pt_dsm_x, pt_dsm_y, pt_dsm_z) 

        if include_trunks
            pt_tsm_x, pt_tsm_y, pt_tsm_z = getsurfdat(copy(tsm_x),copy(tsm_y),copy(tsm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(st.forest_peri))
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
                                                        pts_x[crx],pts_y[crx],pts_e[crx],"trunks",st.image_height)

            else
                pcd2pol2cart!(pt_tsm_x,pt_tsm_y,pt_tsm_z,pts_x[crx],pts_y[crx],pts_e[crx],"trunks",st.image_height)
            end

            if st.phenology == "leafoff"
                prepsurfdat!(pt_tsm_x,pt_tsm_y,pt_tsm_z)
            else
                prepsurfdat!(append!(pt_dsm_x,pt_tsm_x),append!(pt_dsm_y,pt_tsm_y),append!(pt_dsm_z,pt_tsm_z))
            end
        end

        if !st.terrainmask_precalc

            if st.terrain_highres
                # get the high-res local terrain
                pt_hrdtm_x, pt_hrdtm_y, pt_hrdtm_z = getsurfdat(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),
                                                        pts_x[crx],pts_y[crx],pts_e[crx],st.highres_peri);
                pt_hrdtm_x, pt_hrdtm_y = pcd2pol2cart!(ter2rad,pt_hrdtm_x, pt_hrdtm_y, pt_hrdtm_z,pts_x[crx],
                                            pts_y[crx],pts_e[crx],"terrain",rbins_hrdtm,st.image_height)
        
            end

            if st.terrain_lowres 
                # get the low-res regional terrain
                pt_lrdtm_x, pt_lrdtm_y, pt_lrdtm_z = getsurfdat(copy(lrdtm_x),copy(lrdtm_y),copy(lrdtm_z),pts_x[crx],pts_y[crx],pts_e[crx],st.lowres_peri);
                pt_lrdtm_x, pt_lrdtm_y = pcd2pol2cart!(ter2rad,pt_lrdtm_x, pt_lrdtm_y, pt_lrdtm_z,pts_x[crx],pts_y[crx],pts_e_lrdtm[crx],
                                            "terrain",rbins_lrdtm,st.image_height);
        
            end

            if st.horizon_line && !oshd_flag
                pt_hrdtm_x,pt_hrdtm_y = hlm2cart(ter2rad,hlm_tht[:,crx])
            elseif st.terrain_highres && (st.terrain_lowres || (st.horizon_line && oshd_flag))
                prepterdat!(append!(pt_hrdtm_x,pt_lrdtm_x),append!(pt_hrdtm_y,pt_lrdtm_y));
            elseif st.terrain_highres
                prepterdat!(pt_hrdtm_x,pt_hrdtm_y)
            elseif (st.terrain_lowres && !st.terrain_highres)
                pt_hrdtm_x, pt_hrdtm_y = prepterdat(pt_lrdtm_x,pt_lrdtm_y)
            end

        end

        if buildings_flag

            pt_bhm_x, pt_bhm_y, pt_bhm_z = getsurfdat(copy(bhm_x),copy(bhm_y),copy(bhm_z),pts_x[crx],pts_y[crx],pts_e[crx],300);    
            pt_bhm_x, pt_bhm_y =  pcd2pol2cart!(ter2rad,pt_bhm_x,pt_bhm_y,pt_bhm_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                "buildings",rbins_bhm,st.image_height)
            if !terrainmask_precalc
                prepterdat!(append!(pt_hrdtm_x,pt_bhm_x),append!(pt_hrdtm_y,pt_bhm_y));
            else
                pt_hrdtm_x, pt_hrdtm_y = prepterdat(pt_bhm_x,pt_bhm_y)
            end

        end

        if st.step_progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "1. Transferring to polar took "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        ######################
        ### Classify image
        st.step_progress && (start = time())

        st.terrainmask_precalc ? (mat2ev .= terrain_mask[:,:,crx]) : (fill!(mat2ev,1))

        # occupy matrix with surface points
        for zdx = 1:1:size(rbins,1)-1
            ridx = findall(rbins[zdx] .<= pt_dsm_z .< rbins[zdx+1])
            fillmat!(canrad,kdtree,hcat(pt_dsm_x[ridx],pt_dsm_y[ridx]),knum[zdx],mat2ev)
            
            if st.phenology == "leafoff" && include_trunks
                tridx = findall(rbins[zdx] .<= pt_tsm_z .< rbins[zdx+1])
                fillmat!(canrad,kdtree,hcat(pt_tsm_x[tridx],pt_tsm_y[tridx]),knum_t[zdx],mat2ev)
            end
        end

        # add terrain
        !st.terrainmask_precalc && (fillmat!(canrad,kdtree,hcat(pt_hrdtm_x,pt_hrdtm_y),10,mat2ev))

        (st.phenology == "leafoff" && include_trunks) && (fillmat!(canrad,kdtree,hcat(pt_tsm_x,pt_tsm_y),10,mat2ev))

        mat2ev[outside_img] .= 1

        st.save_images && (images["SHI$(outsuf)"][:,:,crx] = mat2ev)

        if st.step_progress
            if crx != 1; try; rm(joinpath(progdir,progtext2)); catch; end; end
            global progtext2 = "2. Classifying image took "*cfmt.("%.$(2)f", time() - start)*" seconds"
            writedlm(joinpath(progdir,progtext2),NaN)
        end

        ########################
        ### Perform calculations
        st.step_progress && (start = time())

        ##### Calculate svf and save
        svf_p, svf_h = calc_svf(canrad,mat2ev)

        dataset["svf_planar$(outsuf)"][crx] = Int(round(svf_p*100))
        dataset["svf_hemi$(outsuf)"][crx]   = Int(round(svf_h*100))

        ##### CalcSWR
        if st.calc_trans
            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],st.time_zone)

            @unpack trans_for = solar
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)

            dataset["tvt$(outsuf)"][:,crx] = Int.(round.((vec(aggregate_data(solar,trans_for))).*100));

            if st.calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p,st.calc_swr)
                dataset["swr_total$(outsuf)"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["swr_direct$(outsuf)"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if st.step_progress
            if crx != 1; try; rm(joinpath(progdir,progtext3)); catch; end; end
            global progtext3 = "3. Calculations and export took "*cfmt.("%.$(2)f", time() - start)*" seconds"
            writedlm(joinpath(progdir,progtext3),NaN)
        end

        # save the progress
        percentdone = Int(floor((crx / size(pts,1)) * 100))
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*cfmt.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end # end crx

    close(dataset)
    st.save_images && close(images)

    (st.save_images && st.make_pngs) && make_SHIs(outdir)

    # Save settings and data information to text file if not running batch
    if !st.batch 
        write_metadata(outdir,outstr,st,fp)
    end

end
