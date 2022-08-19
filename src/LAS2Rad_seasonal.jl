function LAS2Rad_seasonal(pts::Matrix{Float64},dat_in::Dict{String, String},par_in::Dict{String, Any},
    exdir::String,taskID="task")

    ################################################################################
    # > Initialise

    # run compatability check then extract settings
    dat_in, par_in = compatability_check(dat_in,par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    # separate the points to vectors
    pts_x = float(pts[:,1])
    pts_y = float(pts[:,2])

    if progress; start = time(); end

    ################################################################################
    # > Organise the progress reporting

    outdir, outstr, crxstart, append_file, percentdone = organise_outf(taskID,exdir,batch,size(pts_x,1))
    global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crxstart-1)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    ################################################################################
    # > Import surface data

    limits_canopy = hcat((floor(minimum(pts_x))-surf_peri),(ceil(maximum(pts_x))+surf_peri),
                    (floor(minimum(pts_y))-surf_peri),(ceil(maximum(pts_y))+surf_peri))

    dsm_x, dsm_y, dsm_z = readlas(dsmf,limits_canopy)

    ################################################################################
    # > Import and prepare terrain data

    if terrain_highres || terrain_lowres || horizon_line
        terrain = true
    else
        terrain = false
    end

    # load and clip the dtm
    if terrain_highres

        limits_highres = hcat((floor(minimum(pts_x))-surf_peri*2),(ceil(maximum(pts_x))+surf_peri*2),
                        (floor(minimum(pts_y))-surf_peri*2),(ceil(maximum(pts_y))+surf_peri*2))

        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits_highres,true, true)

    end

    # load and clip the dem
    if terrain_lowres

        limits_lowres = hcat((floor(minimum(pts_x))-terrain_peri),(ceil(maximum(pts_x))+terrain_peri),
                        (floor(minimum(pts_y))-terrain_peri),(ceil(maximum(pts_y))+terrain_peri))

        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits_lowres,true,true)
        pts_e_dem = findelev(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,100)

    end

    # determine ground elevation of points
    if size(pts,2) > 2
        pts_e = pts[:,3]
    elseif terrain_highres
        pts_e     = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y)
    elseif terrain_lowres
        pts_e     = pts_e_dem
    else # if there's no terrain input, just assume the ground is flat
        pts_e     = fill(0.0,size(pts_x))
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
        pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x))
        imcX = 0.0; imcY = 0.0
    end

    # get the tile horizon line
    if horizon_line  # get pre-calculated horizon line
        pt_dem_x, pt_dem_y = load_hlm(hlmf,taskID)
    elseif terrain_tile && !horizon_line
        pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,mean(pts_x),mean(pts_y),mean(pts_e_dem),
                                            terrain_peri,"terrain",image_height,0.0,dem_cellsize);
    end

    # load the buildings
    if buildings
    	bhm_x, bhm_y, bhm_z, bhm_cellsize = read_griddata_window(bhmf,limits_canopy,true,true,true)
    	if !isempty(bhm_x)
    		if abs(mode(bhm_z) - mode(dtm_z)) > 40 && terrain_highres # if normalsed
    			bhm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),bhm_x,bhm_y)
    		end
            build = true
    	else; build = false; end
    else; build = false;
    end

    ###############################################################################
    # > Generate extra canopy elements

    limits_trees =  hcat((floor(minimum(pts_x))-surf_peri/2),(ceil(maximum(pts_x))+surf_peri/2),
    				(floor(minimum(pts_y))-surf_peri/2),(ceil(maximum(pts_y))+surf_peri/2))

    # load the dbh
    if trunks
        dbh_x, dbh_y, dbh_z, dbh_r, lastc = loaddbh(dbhf,limits_trees)
        if !isempty(dbh_x)
            dbh_e = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),dbh_x,dbh_y)
            tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,40,0.1,dbh_e)
            trunk = true
        else; trunk = false; end
    else; trunk = false;
    end

    # load the ltc
    if branches || season == "complete"

        # load all the lidar and branch/tree information
        ltc = loadltc_laz(ltcf,limits_canopy,dbh_x,dbh_y,lastc,season)  
 
        if season == "complete"
            # replace the dsm with the ltc points so that they're classified as deciduous or not
            dsm_evergreen, dsm_deciduous = split_points(ltc.x,ltc.y,ltc.z,Bool.(ltc.cls))
            dsm_x = []; dsm_y = []; dsm_z = []
        end

        # get the points that correspond to trees in the tsm data
        inbounds = indexin(ltc.num, lastc) .!= nothing;
        ltc = clipdat(ltc[inbounds,:],limits_trees);

        # remove the points < 1m from the ground
        if abs(mode(ltc.z) - mode(dtm_z)) < 60 # if the data's not normalised, it needs to be normalised)
            ltc.z .-= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),ltc.x,ltc.y)
        end
        ltc = ltc[setdiff(1:end, findall(ltc.z.<1)), :]

        # make the branch/canopy points and classify as deciduous or not
        bsm_x, bsm_y, bsm_z, bsm_d  = make_branches(ltc.x,ltc.y,ltc.z,
                                        ltc.tx,ltc.ty,ltc.hd,ltc.ang,ltc.cls,branch_spacing,season)

        bsm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),bsm_x,bsm_y)

        # split all the data
        bsm_evergreen, bsm_deciduous = split_points(bsm_x,bsm_y,bsm_z,bsm_d)

        # clear up some memory
        ltc = [];
        bsm_x = []; bsm_y = []; bsm_z = [];

    end

    ###############################################################################
    # > Load meteorological data

    if calc_swr == 2
        swrdat   = readdlm(swrf)
        t_step   = Dates.value(Dates.Minute(Dates.DateTime(swrdat[3,1].*" ".*swrdat[3,2],"dd.mm.yyyy HH:MM:SS")-Dates.DateTime(swrdat[2,1].*" ".*swrdat[2,2],"dd.mm.yyyy HH:MM:SS")))
        t_start  = swrdat[1,1].*" ".*swrdat[1,2]
        t_end    = swrdat[end,1].*" ".*swrdat[end,2]
        swr_open = float.(swrdat[:,3])
    else
        t_start = t1; t_end = t2; t_step = tstep;
    end

    ###############################################################################
    # > Image matrix preparation

    # create the empty matrix
    radius = 500
    g_rad, g_coorpol, g_coorcrt, g_img = create_mat(radius)
    g_coorcrt = ((g_coorcrt .- radius) ./ radius) .* 90
    g_img[isnan.(g_rad)] .= 1

    # make g_coorcrt a KDtree for easy look up
    if ~tilt
        kdtree = scipyspat.cKDTree(g_coorcrt)
        kdtreedims = size(g_coorcrt,1)
    end

    ###############################################################################
    # > get the time stamp variable

    loc_time     = collect(Dates.DateTime(t_start,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t_end,"dd.mm.yyyy HH:MM:SS"))
    loc_time_agg = collect(Dates.DateTime(t_start,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t_end,"dd.mm.yyyy HH:MM:SS"))

    if season == "complete"
        # index loc_time for the start/end of the seasonal transition periods
        sdx = Vector{Int32}() # create a seasonal index of the loc_time variable

        midnights = findall((hour.(loc_time) .== 00) .&& (minute.(loc_time) .== 00) .&& (second.(loc_time) .== 00))

        spring1 = findall(loc_time .== Dates.DateTime(spring_start,"dd.mm.yyyy HH:MM:SS"))[1]
        spring2 = findall(loc_time .== (Dates.DateTime(spring_start,"dd.mm.yyyy HH:MM:SS")+Dates.Day(spring_duration)))[1]

        append!(sdx,[1,spring1]) # day 1 + start of spring
        append!(sdx,midnights[spring1 .< midnights .< spring2])
        append!(sdx,spring2) # first day of summer

        autumn1 = findall(loc_time .== Dates.DateTime(autumn_start,"dd.mm.yyyy HH:MM:SS"))[1]
        autumn2 = findall(loc_time .== (Dates.DateTime(autumn_start,"dd.mm.yyyy HH:MM:SS")+Dates.Day(autumn_duration)))[1]

        append!(sdx,autumn1) # start of autumn
        append!(sdx,midnights[autumn1 .< midnights .< autumn2])
        append!(sdx,autumn2) #start of winter

        sdx_longer = vcat(sdx,size(loc_time,1))

        # # do the time vector for the different sky-view fraction data
        # loc_time_Vf  = Vector{DateTime}()
        # loc_time_Vf = [loc_time[1],loc_time[spring1]]
        # append!(loc_time_Vf,loc_time[midnights[spring1 .< midnights .<= spring2]])
        # append!(loc_time_Vf,loc_time[midnights[autumn1 .<= midnights .<= autumn2]])
        
    end

    ###############################################################################
    # > organise the output folder

    # create the output files
    if calc_trans
        dataset      = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file,loc_time_agg,time_zone,season)
    else
        dataset      = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file)
    end

    if save_images
        images = create_exmat(outdir,outstr,pts,g_img,append_file,season)
    end

    # initialise progress reporting
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

    ###############################################################################
    # > get constants

    rbins = collect(0:(surf_peri-0)/5:surf_peri)

    if calc_trans
        drad, im_centre, lens_profile_tht, lens_profile_rpix, trans_for = get_constants(g_img,loc_time)
    end

    ###############################################################################
    # > Loop through the points

    @simd for crx = crxstart:size(pts_x,1)
        if progress; start = time(); end

        #### transfer point clouds to polar c
        if branches || season == "complete"
            
            # keep deciduous points separate
            pt_bsm_d_x, pt_bsm_d_y, pt_bsm_d_z = getsurfdat(bsm_deciduous.x,bsm_deciduous.y,bsm_deciduous.z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri)

            pt_dsm_d_x, pt_dsm_d_y, pt_dsm_d_z = getsurfdat(dsm_deciduous.x,dsm_deciduous.y,dsm_deciduous.z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri)

            pt_bsm_d_x, pt_bsm_d_y, pt_bsm_d_z = pcd2pol2cart(pt_bsm_d_x,pt_bsm_d_y,pt_bsm_d_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                        surf_peri,"surface",image_height,pts_slp[crx],0);

            pt_dsm_d_x, pt_dsm_d_y, pt_dsm_d_z = pcd2pol2cart(pt_dsm_d_x,pt_dsm_d_y,pt_dsm_d_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                        surf_peri,"surface",image_height,pts_slp[crx],0);

            # combined dsm points are all evergreen
            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_evergreen.x,dsm_evergreen.y,dsm_evergreen.z,bsm_evergreen.x,bsm_evergreen.y,bsm_evergreen.z,
                                                        pts_x[crx],pts_y[crx],pts_e[crx],surf_peri)

        elseif branches # combine into one dataset for easier plotting

            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_x,dsm_y,dsm_z,bsm_x,bsm_y,bsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri)
        
        else

            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_x,dsm_y,dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri);
        
        end

        pt_dsm_x, pt_dsm_y, pt_dsm_z = pcd2pol2cart(pt_dsm_x,pt_dsm_y,pt_dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                        surf_peri,"surface",image_height,pts_slp[crx],0);

        if trunk
            pt_tsm_x, pt_tsm_y, pt_tsm_z = getsurfdat(tsm_x,tsm_y,tsm_z,pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5))
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
                pt_tsm_x, pt_tsm_y, _ = pcd2pol2cart(append!(pt_tsm_x,tsm_tmp[1]),append!(pt_tsm_y,tsm_tmp[2]),append!(pt_tsm_z,tsm_tmp[3]),
                                                        pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5),"surface",image_height,pts_slp[crx],0);

            else
                pt_tsm_x, pt_tsm_y, _ = pcd2pol2cart(pt_tsm_x,pt_tsm_y,pt_tsm_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                        Int.(surf_peri*0.5),"surface",image_height,pts_slp[crx],0);
            end
        end

        if terrain
            if terrain_highres
                pt_dtm_x, pt_dtm_y =  pcd2pol2cart(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],
                                                    Int.(300),"terrain",image_height,pts_slp[crx],dtm_cellsize);
            end

            if terrain_lowres && !terrain_tile
                pt_dem_x, pt_dem_y = pcd2pol2cart(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e_dem[crx],
                                                    terrain_peri,"terrain",image_height,pts_slp[crx],dem_cellsize);
            end

            if terrain_highres && (terrain_lowres || horizon_line)
                pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
            elseif terrain_highres
                pt_dtm_x, pt_dtm_y = prepterdat(pt_dtm_x,pt_dtm_y)
            elseif terrain_lowres && !terrain_highres
                pt_dtm_x, pt_dtm_y = prepterdat(pt_dem_x,pt_dem_y)
            end
        end

        if build
            pt_bhm_x, pt_bhm_y =  pcd2pol2cart(copy(bhm_x),copy(bhm_y),copy(bhm_z),pts_x[crx],pts_y[crx],pts_e[crx],
                                                Int.(50),"terrain",image_height,pts_slp[crx],bhm_cellsize);
            if terrain
                pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_bhm_x),append!(pt_dtm_y,pt_bhm_y));
            else
                pt_dtm_x, pt_dtm_y = prepterdat(pt_bhm_x,pt_bhm_y)
            end
        end

        if tilt # disabled in Preparatory_Functions.jl
            # transform image matrix
            imcX, imcY = getimagecentre(pts_slp[crx],pts_asp[crx])
            gcrcrt = fill(NaN,size(g_coorcrt))
            gcrcrt[:,1] = g_coorcrt[:,1] .+ imcX
            gcrcrt[:,2] = g_coorcrt[:,2] .+ imcY
            kdtree = scipyspat.cKDTree(gcrcrt)
            kdtreedims = size(gcrcrt,1)
        else
            imcX = 0.0; imcY = 0.0;
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "1. Transferring to polar took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        ######################
        ### Classify image
        if progress; start = time(); end

        sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts_x[crx],pts_y[crx],loc_time,time_zone,coor_system)

        mat2ev  = copy(g_img);

        # create the base winter map
        # add terrain
        if terrain
            mat2ev = fillmat(kdtree,hcat(pt_dtm_x,pt_dtm_y),1.0,kdtreedims,10,radius,mat2ev);
        end

        # add trunks
        if trunk
            mat2ev = fillmat(kdtree,hcat(pt_tsm_x,pt_tsm_y),1.0,kdtreedims,15,radius,mat2ev);
        end

        # add the evergreen points
        tol   = tolerance[2].*collect(reverse(0.5:0.1:1))
        for zdx = 1:1:size(rbins,1)-1
            ridx = findall(rbins[zdx] .<= pt_dsm_z .< rbins[zdx+1])
            imdx = reshape(findpairs(kdtree,hcat(pt_dsm_x[ridx],pt_dsm_y[ridx]),tol[zdx],kdtreedims,40),(radius*2,radius*2))
            mat2ev[imdx.==1] .= 0
        end

        # add the deciduous points
        tol   = tolerance[1].*collect(reverse(0.75:0.05:1))
        for zdx = 1:1:size(rbins,1)-1
            ridx = findall(rbins[zdx] .<= pt_dsm_d_z .< rbins[zdx+1])
            mat2ev = fillmat(kdtree,hcat(pt_dsm_d_x[ridx],pt_dsm_d_y[ridx]),tol[zdx],kdtreedims,20,radius,mat2ev);
        end

        mat2ev[isnan.(g_rad)] .= 1;

        # calculate transmissivity for the start and end of the season
        Vf_p = Vector{Float64}(undef,size(sdx,1))
        Vf_h = Vector{Float64}(undef,size(sdx,1))
        Vf_p[1], Vf_h[1] = calculateVf(mat2ev,g_rad,radius)
        Vf_p[end] = Vf_p[1]
        Vf_h[end] = Vf_h[1]
    
        transfor = Vector{Float64}(undef,size(loc_time,1))
        transfor[sdx[1]:sdx[2]] = calc_transmissivity(float.(mat2ev),loc_time[sdx[1]:sdx[2]],tstep,radius,sol_phi[sdx[1]:sdx[2]],sol_tht[sdx[1]:sdx[2]],g_coorpol,0.0,0.0,
                        im_centre,trans_for[sdx[1]:sdx[2]],lens_profile_tht,lens_profile_rpix)
        transfor[sdx[end]:end] = calc_transmissivity(float.(mat2ev),loc_time[sdx[end]:end],tstep,radius,sol_phi[sdx[end]:end],sol_tht[sdx[end]:end],g_coorpol,0.0,0.0,
                        im_centre,trans_for[sdx[end]:end],lens_profile_tht,lens_profile_rpix)

        temp_img = Array{Int8}(undef,(1000,1000,size(sdx,1)))
        temp_img .= 0
        temp_img[:,:,1] = Bool.(copy(mat2ev))

        tol_range = (tolerance[2] - tolerance[1])/2
        tol = tolerance[1]:(tol_range/(spring_duration/2)):tolerance[2]

        # sub-sample the bsm and dsm points for randomly allocating to the image matrix
        bx = Int.(floor.(rand(Uniform(1,ceil(spring_duration/2)+1),size(pt_bsm_d_x,1))))
        dx = Int.(floor.(rand(Uniform(1,ceil(spring_duration/2)+1),size(pt_dsm_d_x,1))))

        thres = Int.(ceil(spring_duration/2)) # when to ramp up the point increases
        # start the loop for spring and plot pt_bsm_d_x
        start_spring = findall(sdx .== spring1)[1]
        start_summer = findall(sdx .== spring2)[1]
        for tdx = start_spring:start_summer ### what to loop through?

            if tdx-1 <= thres
                # add the deciduous points
                tol2   = tol[tdx].*collect(reverse(0.5:0.1:1))
                for zdx = 1:1:size(rbins,1)-1
                    ridx = (rbins[zdx] .<= pt_bsm_d_z .< rbins[zdx+1]) .* (bx.==tdx-1)
                    mat2ev = fillmat(kdtree,hcat(pt_bsm_d_x[ridx],pt_bsm_d_y[ridx]),tol2[zdx],kdtreedims,40,radius,mat2ev);
                end
            else
                # start making them bigger
                tol2   = tolerance[2].*collect(reverse(0.5:0.1:1))
                for zdx = 1:1:size(rbins,1)-1
                    ridx = (rbins[zdx] .<= pt_bsm_d_z .< rbins[zdx+1]) .* (bx.==tdx-1-thres)
                    mat2ev = fillmat(kdtree,hcat(pt_bsm_d_x[ridx],pt_bsm_d_y[ridx]),tol2[zdx],kdtreedims,40,radius,mat2ev);
                end
                mat2ev = fillmat(kdtree,hcat(pt_dsm_d_x[dx.==(tdx-1-thres)],pt_dsm_d_y[dx.==(tdx-1-thres)]),tol2[1],kdtreedims,20,radius,mat2ev);
            end

            mat2ev[isnan.(g_rad)] .= 1;

            # Calculate Vf
            Vf_p[tdx], Vf_h[tdx] = calculateVf(mat2ev,g_rad,radius)

            # Calculate SWR/forest transmissivity
            if calc_trans
                transfor[sdx[tdx]:sdx[tdx+1]] = calc_transmissivity(float.(mat2ev),loc_time[sdx[tdx]:sdx[tdx+1]],tstep,radius,sol_phi[sdx[tdx]:sdx[tdx+1]],
                                        sol_tht[sdx[tdx]:sdx[tdx+1]],g_coorpol,0.0,0.0,
                                        im_centre,trans_for[sdx[tdx]:sdx[tdx+1]],lens_profile_tht,lens_profile_rpix)
            end

            temp_img[:,:,tdx] = Bool.(copy(mat2ev))

        end

        # do the loop for autumn (inversely)
        tol = tolerance[1]:(tol_range/(autumn_duration/2)):tolerance[2]
        bx = Int.(floor.(rand(Uniform(1,ceil(autumn_duration/2)+1),size(pt_bsm_d_x,1))))
        dx = Int.(floor.(rand(Uniform(1,ceil(autumn_duration/2)+1),size(pt_dsm_d_x,1))))

        # start the loop for autumn and plot pt_bsm_d_x
        mat2ev = copy(Int64.(temp_img[:,:,1]))
        start_autumn = findall(sdx .== autumn1)[1]
        start_winter = findall(sdx .== autumn2)[1]
        thres = Int.(ceil(autumn_duration/2)) # when to ramp up the point increases

        adx = start_winter-1:-1:start_autumn # get an inverse index

        for tdx = 1:1:size(adx,1)

            if tdx < thres
                # add the deciduous points
                tol2   = tol[tdx].*collect(reverse(0.5:0.1:1))
                for zdx = 1:1:size(rbins,1)-1
                    ridx = (rbins[zdx] .<= pt_bsm_d_z .< rbins[zdx+1]) .* (bx.==tdx)
                    mat2ev = fillmat(kdtree,hcat(pt_bsm_d_x[ridx],pt_bsm_d_y[ridx]),tol2[zdx],kdtreedims,40,radius,mat2ev);
                end
            else
                # start making them bigger
                tol2   = tolerance[2].*collect(reverse(0.5:0.1:1))
                for zdx = 1:1:size(rbins,1)-1
                    ridx = (rbins[zdx] .<= pt_bsm_d_z .< rbins[zdx+1]) .* (bx.==tdx-thres)
                    mat2ev = fillmat(kdtree,hcat(pt_bsm_d_x[ridx],pt_bsm_d_y[ridx]),tol2[zdx],kdtreedims,40,radius,mat2ev);
                end
                mat2ev = fillmat(kdtree,hcat(pt_dsm_d_x[dx.==(tdx+1-thres)],pt_dsm_d_y[dx.==(tdx+1-thres)]),tol2[1],kdtreedims,20,radius,mat2ev);
            end

            mat2ev[isnan.(g_rad)] .= 1;

            ndx = adx[tdx] # new index for tranisitioning from winter back through autumn
            # Calculate Vf
            Vf_p[ndx], Vf_h[ndx] = calculateVf(mat2ev,g_rad,radius)

            # Calculate SWR/forest transmissivity
            if calc_trans
                transfor[sdx[ndx-1]:sdx[ndx]] = calc_transmissivity(float.(mat2ev),loc_time[sdx[ndx-1]:sdx[ndx]],tstep,radius,sol_phi[sdx[ndx-1]:sdx[ndx]],
                                        sol_tht[sdx[ndx-1]:sdx[ndx]],g_coorpol,0.0,0.0,
                                        im_centre,trans_for[sdx[ndx-1]:sdx[ndx]],lens_profile_tht,lens_profile_rpix)
            end

            temp_img[:,:,ndx] = Bool.(copy(mat2ev))

        end

        # put the winter image at the end
        temp_img[:,:,start_winter] = Bool.(copy(temp_img[:,:,1]))

        # Vf same size as loc_time and transfor
        Vf_p_all = Vector{Float64}(undef,size(loc_time,1))
        Vf_h_all = Vector{Float64}(undef,size(loc_time,1))
        for nx = 1:1:size(sdx_longer,1)-1
            Vf_p_all[sdx_longer[nx]:sdx_longer[nx+1]] .= Vf_p[nx]
            Vf_h_all[sdx_longer[nx]:sdx_longer[nx+1]] .= Vf_h[nx]
        end

        if calc_swr == 1
            swrtot, swrdir, _ = calculateSWR_seasonal(transfor,sol_sinelev,sol_tht,sol_phi,loc_time,max.(1367*sol_sinelev,0),Vf_p_all)
        elseif calc_swr == 2
            swrtot, swrdir, _ = calculateSWR(transfor,sol_sinelev,sol_tht,sol_phi,loc_time,swr_open,Vf_p)
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext3)); catch; end; end
            global progtext3 = "2. Classifying image and transmissivity calculations took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext3),NaN)
        end

        ######################
        # export the data [append to netcdf]
        if progress; start = time(); end

        if season == "complete"
            dataset["Vf_planar"][:,crx] = Int.(round.((vec(aggregate_data(loc_time,loc_time_agg,Vf_p_all,tstep))).*100))
            dataset["Vf_hemi"][:,crx]   = Int.(round.((vec(aggregate_data(loc_time,loc_time_agg,Vf_h_all,tstep))).*100))
        else
            dataset["Vf_planar"][crx] = Int(round(Vf_p*100))
            dataset["Vf_hemi"][crx]   = Int(round(Vf_h*100))
        end

        if calc_trans
            dataset["Forest_Transmissivity"][:,crx] = Int.(round.((vec(aggregate_data(loc_time,loc_time_agg,transfor,tstep))).*100))
            if calc_swr > 0
                if calc_swr == 2
                    swrtot[isnan.(swrtot)] .= -9999
                    swrdir[isnan.(swrdir)] .= -9999
                    # note this accommodates data gaps in measured above canopy radiation, and resultant NaN values in output modelled values.
                    # if calc_swr = 2, loc_time_agg and loc_time are the same, so there is no implication in the next lines if setting
                    # values to -9999
                end
                dataset["SWR_total"][:,crx] = Int.(round.(vec(aggregate_data(loc_time,loc_time_agg,swrtot,tstep))))
                dataset["SWR_direct"][:,crx] = Int.(round.(vec(aggregate_data(loc_time,loc_time_agg,swrdir,tstep))))
            end
        end

        if save_images
            if season == "complete"

                 # save the gif of the seasonal change
                if !ispath(joinpath(outdir,"gifs")); mkdir(joinpath(outdir,"gifs")); end
                outf = joinpath(outdir,"gifs","SHI_"*sprintf1.("%07.$(2)f",pts_x[crx])*"_"*sprintf1.("%07.$(2)f",pts_y[crx])*".gif")
                save(outf,colorview(Gray,float.(temp_img)))

                # save the winter and summer end members
                images["SHI"][:,:,crx,1] = temp_img[:,:,start_winter]
                images["SHI"][:,:,crx,1] = temp_img[:,:,start_summer]

            else
                images["SHI"][:,:,crx] = mat2ev
            end
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext4)); catch; end; end
            global progtext4 = "3. Exporting data took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext4),NaN)
        end

        # save the progress
        percentdone = Int(floor((crx / size(pts,1)) * 100))
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end # end crx

    close(dataset)
    if save_images; close(images); end

    println("done with "*taskID)

	return dat_in, par_in

end
