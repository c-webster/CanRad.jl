function CHM2Rad(pts,dat_in,par_in,exdir,taskID="task")

    ################################################################################
    # Initialise

    # run compatability check then extract settings
    dat_in, par_in = compatability_check(dat_in,par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    # separate the points to vectors
    pts_x = float(pts[:,1])
    pts_y = float(pts[:,2])

    if progress; start = time(); end

    ################################################################################
    # > Import surface data and clip

    # get analysis limits
    limits_canopy = hcat((floor(minimum(pts_x))-surf_peri),(ceil(maximum(pts_x))+surf_peri),
                    (floor(minimum(pts_y))-surf_peri),(ceil(maximum(pts_y))+surf_peri))

    chm_x, chm_y, chm_z, chm_cellsize = read_griddata_window(chmf,limits_canopy,true,true)

    # if @isdefined(cbhf)
        # chb_x, chb_y, chm_b, _ = read_griddata(cbhf,true)
        # _, _, chm_b, _    = clipdat(copy(chb_x),copy(chb_y),chm_b,limits_canopy)
    # else
        chm_b = fill(2.0,size(chm_z))
        # chm_b[chm_z .>= 2] .= 2.0
    # end

    # load the trunk data
    if trunks
        dbh_x, dbh_y, dbh_z, dbh_r = loaddbh(dbhf,limits_canopy,50)
        if !isempty(dbh_x)
            dbh_e = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),dbh_x,dbh_y)
            tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,30,0.1,dbh_e)
            trunk = true
        else; trunk = false; end
    else; trunk = false; end

    ################################################################################
    # > Import, clip and prepare terrain data

    # set terrain = true if any terrain is to be plotted in the final image matrix
    if terrain_highres || terrain_lowres || horizon_line
        terrain = true
    else
        terrain = false
    end

    if terrain_highres

        limits_highres = hcat((floor(minimum(pts_x))-surf_peri*2),(ceil(maximum(pts_x))+surf_peri*2),
                        (floor(minimum(pts_y))-surf_peri*2),(ceil(maximum(pts_y))+surf_peri*2))

        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits_highres,true, true)

    end

    # Get the low-res terrain data
    if terrain_lowres || !horizon_line

        limits_lowres = hcat((floor(minimum(pts_x))-terrain_peri),(ceil(maximum(pts_x))+terrain_peri),
                        (floor(minimum(pts_y))-terrain_peri),(ceil(maximum(pts_y))+terrain_peri))

        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits_lowres,true,true)
        pts_e_dem = findelev(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,100)

    end

    # Get elevation of surface and evaluation points
    if terrain_highres
        pts_e     = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y)
        chm_e     = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),chm_x,chm_y)
    elseif terrain_lowres
        pts_e     = pts_e_dem
        chm_e     = findelev(copy(dem_x),copy(dem_y),copy(dem_z),chm_x,chm_y)
    else # if there's no terrain input, just assume the ground is flat
        pts_e     = fill(0.0,size(pts_x))
        chm_e     = fill(0.0,size(chm_x))
    end

    chm_b .+= chm_e
    chm_z .+= chm_e

    # Calculate slope and aspect from dtm data
    if tilt # currently disabled in Preparatory_Functions.jl
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
        pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = fill(0.0,size(pts_x,1))
        imcX = 0.0; imcY = 0.0
    end

    # get the tile horizon line
    if horizon_line # get pre-calculated horizon line
        pt_dem_x, pt_dem_y = load_hlm(hlmf,taskID)
    elseif terrain_tile && !horizon_line
        pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,mean(pts_x),mean(pts_y),mean(pts_e_dem),
                                            terrain_peri,"terrain",image_height,0.0,dem_cellsize);
    end

    # load the buildings
    if buildings
        bhm_x, bhm_y, bhm_z, bhm_cellsize = read_griddata_window(bhmf,limits_canopy,true,true)
        if !isempty(bhm_x)
            if abs(mode(bhm_z) - mode(dtm_z)) > 40 && terrain_highres # if normalsed
                bhm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),bhm_x,bhm_y)
            end
            build = true
    	else; build = false; end
    else; build = false; end

    ###############################################################################
    # Get crown volume density information

    if OSHD

        # load forest mix ratio and correct for lavd
        mr_x, mr_y, mr_val, _ = read_griddata_window(mrdf,limits_canopy,true,true)
        mr_val[mr_val .== 0] .= 1.0

        limits_tile = hcat((floor(minimum(pts_x))),(ceil(maximum(pts_x))),
                        (floor(minimum(pts_y))),(ceil(maximum(pts_y))))
        for_type = mode(read_griddata_window(fcdf,limits_tile,true,true)[3])

        # make sure areas in Jura and central Switzerland become alpine forests
        if median(pts_e) > 1000
            for_type = 2
        end

        if for_type == 1
            temp = reverse(collect(0.02:(0.2-0.02)/9999:0.2))
        elseif for_type == 2
            temp = reverse(collect(0.05:(0.67-0.05)/9999:0.67))
        end
        lavd_val = fill(0.0,size(mr_val))

        for vx in eachindex(mr_val)
            lavd_val[vx] = temp[Int.(mr_val[vx])]
        end

        chm_lavd = findelev(copy(mr_x),copy(mr_y),copy(lavd_val),chm_x,chm_y,10,"cubic")

    else

        ch_x, ch_y, chm_lavd, _ = read_griddata(lavdf,true,true)
        chm_lavd = clipdat(copy(ch_x),copy(ch_y),chm_lavd,limits_canopy)[3]

    end

    ###############################################################################
    # Load meteorological data

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
    # > Tile preparation

    # create an empty image matrix
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
    # > create the output files

    if OSHD
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3]
        global outdir = exdir*"/"*outstr
    elseif batch
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3][1:end-4]
        global outdir = exdir*"/"*outstr
    else
        outstr = String(split(exdir,"/")[end-1])
        global outdir = exdir
    end

    if !ispath(outdir)
        mkpath(outdir)
    end

    # set start point within tile
    fx = readdir(outdir)[findall(startswith.(readdir(outdir),"Processing"))]
    if isempty(fx)
        crxstart = 1; append_file = false
    else
        crxstart = parse(Int,split(fx[1])[4])
        if crxstart == size(pts_x,1) # if tile is complete, restart tile
            crxstart = 1
            append_file = false
        else
            append_file = true
        end
    end

    # create the output files
    if calc_trans
        loc_time     = collect(Dates.DateTime(t_start,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t_end,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t_start,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t_end,"dd.mm.yyyy HH:MM:SS"))
        dataset      = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file,loc_time_agg)
    else
         dataset     = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file)
    end

    if save_images
        images = create_exmat(outdir,outstr,pts,g_img,append_file)
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
    # get constants
    rbins = collect(0:(surf_peri-0)/5:surf_peri)
    tol   = tolerance.*collect(reverse(0.75:0.05:1))

    if calc_trans
        drad, im_centre, lens_profile_tht, lens_profile_rpix, trans_for = get_constants(g_img,loc_time)
    end

    ###############################################################################
    # > Loop through the points

    @simd for crx = crxstart:size(pts_x,1)

        ######################
        #### calculate horizon lines

        if progress; start = time(); end

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
            elseif (terrain_lowres && !terrain_highres)
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

        if pt_corr
            pt_chm_x, pt_chm_y, pt_chm_r, pt_chm_x_thick, pt_chm_y_thick = calcCHM_Ptrans(copy(chm_x),copy(chm_y),copy(chm_z),
                                copy(chm_b),copy(chm_lavd),pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,image_height,chm_cellsize) # calculated points

            pt_chm_x_pts, pt_chm_y_pts, pt_chm_r_pts = pcd2pol2cart(copy(chm_x),copy(chm_y),copy(chm_z),pts_x[crx],pts_y[crx],
                                pts_e[crx],surf_peri,"chm",image_height,pts_slp[crx],chm_cellsize) # pts from the CHM

            if trunk
                pt_tsm_x, pt_tsm_y, pt_tsm_z = getsurfdat(tsm_x,tsm_y,tsm_z,pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5))
                tidx = findall(dist(dbh_x,dbh_y,pts[crx,1],pts[crx,2]) .< 4)
                if size(tidx,1) > 0
                    hdt  = dist(dbh_x[tidx],dbh_y[tidx],pts[crx,1],pts[crx,2])
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
                    pt_tsm_x, pt_tsm_y, _ = pcd2pol2cart(pt_tsm_x,pt_tsm_y,pt_tsm_z,
                                                        pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5),"surface",image_height,pts_slp[crx],0);
                end
            end

        else
            #  100% opaque canopy:
            pt_chm_x, pt_chm_y = pcd2pol2cart(copy(chm_x),copy(chm_y),copy(chm_z),pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,"terrain",image_height,pts_slp[crx],chm_cellsize)
            if terrain # merge if using terrain
                pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_chm_x,pt_dtm_x),append!(pt_chm_y,pt_dtm_y));
            else
                pt_dtm_x, pt_dtm_y = prepterdat(pt_chm_x,pt_chm_y);
            end
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

        mat2ev  = copy(g_img);

        # occupy matrix
        if pt_corr
            for zdx = 1:1:size(rbins,1)-1
                ridx   = findall(rbins[zdx] .<= pt_chm_r .< rbins[zdx+1])
                mat2ev = fillmat(kdtree,hcat(pt_chm_x[ridx],pt_chm_y[ridx]),tol[zdx],kdtreedims,30,radius,mat2ev);
            end
            # include canopy surface points
            # mat2ev = fillmat(kdtree,hcat(pt_chm_x_pts[pt_chm_r_pts .> 10],pt_chm_y_pts[pt_chm_r_pts .> 10]),4.0,kdtreedims,30,radius,mat2ev); # include canopy surface points
            if terrain
                # mat2ev = fillmat(kdtree,hcat(vcat(pt_chm_x_thick,pt_dtm_x),vcat(pt_chm_y_thick,pt_dtm_y)),1.5,kdtreedims,10,radius,mat2ev); # distance canopy is opaque and treated with terrain
                mat2ev = fillmat(kdtree,hcat(pt_dtm_x,pt_dtm_y),1.5,kdtreedims,10,radius,mat2ev); # distance canopy is opaque and treated with terrain
            else
                mat2ev = fillmat(kdtree,hcat(pt_chm_x_thick,pt_chm_y_thick),1.5,kdtreedims,10,radius,mat2ev); # distance canopy is opaque and treated with terrain
            end
        else # treat all canopy as opaque (like terrain) -> all points come in as terrain points from above section
            mat2ev = fillmat(kdtree,hcat(pt_dtm_x,pt_dtm_y),1.0,kdtreedims,10,radius,mat2ev); # use this line if plotting opaque canpoy
        end

        # add trunks
        if trunk
            mat2ev = fillmat(kdtree,hcat(pt_tsm_x,pt_tsm_y),2.0,kdtreedims,15,radius,mat2ev)
        end

        mat2ev[isnan.(g_rad)] .= 1;

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext2)); catch; end; end
            global progtext2 = "2. Classifying image took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext2),NaN)
        end

        ########################
        ### Perform calculations

        if progress; start = time(); end

        ##### Calculate Vf
        Vf_p, Vf_h = calculateVf(mat2ev,g_rad,radius)

        ##### CalcSWR
        if calc_trans
            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts_x[crx],pts_y[crx],loc_time,time_zone,coor_system)
            transfor = calc_transmissivity(float.(mat2ev),loc_time,tstep,radius,sol_phi,sol_tht,g_coorpol,0.0,0.0,
                                    im_centre,trans_for,lens_profile_tht,lens_profile_rpix)
            if calc_swr == 1
                swrtot, swrdir, _ = calculateSWR(transfor,sol_sinelev,sol_tht,sol_phi,loc_time,max.(1367*sol_sinelev,0),Vf_p)
            elseif calc_swr == 2
                swrtot, swrdir, _ = calculateSWR(transfor,sol_sinelev,sol_tht,sol_phi,loc_time,swr_open,Vf_p)
            end
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext3)); catch; end; end
            global progtext3 = "3. Vf and SWR calculations took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext3),NaN)
        end

        ######################
        #export the data [append to netcdf]
        if progress; start = time(); end

        dataset["Vf_planar"][crx] = Int(round(Vf_p*100))
        dataset["Vf_hemi"][crx]   = Int(round(Vf_h*100))

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
            images["SHI"][:,:,crx] = mat2ev
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext4)); catch; end; end
            global progtext4 = "4. Exporting data took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext4),NaN)
        end

        # save the progress
        percentdone = Int(floor((crx / size(pts,1)) * 100))
        try
            rm(outdir*"/"*outtext)
        catch
            for f in readdir(outdir)
                    if startswith.(f,"Processing")
                            rm(outdir*"/"*f)
                    end
            end
        end
        global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end # end crx

    close(dataset)
    if save_images; close(images); end

    println("done with "*taskID)

    return dat_in, par_in

end
