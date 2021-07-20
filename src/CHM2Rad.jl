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

    # set start point within tile
    crxstart = 1; append_file  = false
    # currently hard coded to always start at 1
    # this should be updated eventually to allow appending to a file and starting part way through a tile
    # if size(readdir(outdir),1) == 5
    #     crxstart = parse(Int,(split(reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])))[4])+1
    #     global outtext  = reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])
    #     append   = true
    # else; crxstart = 1; append_file  = false
    # end

    ################################################################################
    # > Import surface data and clip

    chm_x, chm_y, chm_z, chm_cellsize = read_ascii(chmf,true,true)

    _, _, chm_lavd, _ = read_ascii(lavdf,true,true)

    if @isdefined(cbhf)
        _, _, chm_b, _ = read_ascii(cbhf,true)
    else
        chm_b = fill(0.0,size(chm_z))
        chm_b[chm_z .>= 2] .= 2.0
    end

    # get analysis limits
    limits = Int.(floor.(hcat(vcat(minimum(pts_x),maximum(pts_x)),vcat(minimum(pts_y),maximum(pts_y)))))

    # clip the CHM data to the analysis limits + surf_peri
    _, _, chm_lavd, _ = clipdat(copy(chm_x),copy(chm_y),chm_lavd,limits,surf_peri)
    _, _, chm_b, _ = clipdat(copy(chm_x),copy(chm_y),chm_b,limits,surf_peri)
    chm_x, chm_y, chm_z, _ = clipdat(chm_x,chm_y,chm_z,limits,surf_peri)

    # load the trunk data
    if trunks
        dbh_x, dbh_y, dbh_z, dbh_r = loaddbh(dbhf,limits,50)
        if !isempty(dbh_x)
            dbh_e = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),dbh_x,dbh_y)
            tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,30,0.1,dbh_e)
            trunks_2 = true
        else
            trunks_2 = false
        end
    else
        trunks_2 = false
    end

    ################################################################################
    # > Import, clip and prepare terrain daa

    # set terrain = true if any terrain is to be plotted in the final image matrix
    if terrain_highres || terrain_lowres || horizon_line
        terrain = true
    else
        terrain = false
    end

    if terrain_highres
        if tilt
            dtm_x, dtm_y, dtm_z, dtm_s, dtm_a, dtm_cellsize = importdtm(dtmf,tilt)
        else
            dtm_x, dtm_y, dtm_z, dtm_cellsize = importdtm(dtmf,tilt)
        end
        dtm_x, dtm_y, dtm_z, _ = clipdat(dtm_x,dtm_y,dtm_z,limits,surf_peri*4)
    end

    if terrain_lowres
        dem_x, dem_y, dem_z, dem_cellsize = read_ascii(demf,true)
        dem_x, dem_y, dem_z, _ = clipdat(dem_x,dem_y,dem_z,limits,terrain_peri)
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

    chm_b     = chm_b .+ chm_e
    chm_z     = chm_z .+ chm_e

    # Calculate slope and aspect from dtm data
    if tilt
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
        pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = fill(0.0,size(pts_x,1))
        imcX = 0.0; imcY = 0.0
    end

    ###############################################################################
    # Load meteorological data

    if calc_swr == 2
        swrdat   = readdlm(swrf)
        t_step    = Dates.value(Dates.Minute(Dates.DateTime(swrdat[3,1].*" ".*swrdat[3,2],"dd.mm.yyyy HH:MM:SS")-Dates.DateTime(swrdat[2,1].*" ".*swrdat[2,2],"dd.mm.yyyy HH:MM:SS")))
        t_start       = swrdat[1,1].*" ".*swrdat[1,2]
        t_end       = swrdat[end,1].*" ".*swrdat[end,2]
        swr_open = float.(swrdat[:,3])
    else
        t_start = t1; t_end = t2; t_step = tstep;
    end

    ###############################################################################
    # > Tile preparation

    if horizon_line
        pt_dem_x = readdlm(hldir*"/HL_"*taskID)[:,1]
        pt_dem_y = readdlm(hldir*"/HL_"*taskID)[:,2]
    elseif terrain_tile && !horizon_line
        pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,mean(pts_x),mean(pts_y),mean(pts_e_dem),terrain_peri,"terrain",ch,0.0,dem_cellsize);
    end

    # create an empty matrix for Vf calculation
    g_rad, g_coorpol, g_coorcrt, g_img = create_mat(radius)
    g_coorcrt = ((g_coorcrt .- radius) ./ radius) .* 90
    g_img[isnan.(g_rad)] .= 1

    # make g_coorcrt a KDtree for easy look up
    if ~tilt
        kdtree = scipyspat.cKDTree(g_coorcrt)
        kdtreedims = size(g_coorcrt,1)
    end

    # get constants
    rbins = collect(0:(surf_peri-0)/5:surf_peri)
    tol   = tolerance.*collect(reverse(0.75:0.05:1))

    if calc_trans
        drad, im_centre, lens_profile_tht, lens_profile_rpix, trans_for = get_constants(g_img,loc_time)
    end

    ###############################################################################
    # > create the output files

    if batch
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3][1:end-4]
        global outdir = exdir*"/"*outstr
    else
        outstr = String(split(exdir,"/")[end-1])
        global outdir = exdir
    end

    if !ispath(outdir)
        mkpath(outdir)
    end

    if calc_trans
        loc_time = collect(Dates.DateTime(t_start,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t_end,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t_start,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t_end,"dd.mm.yyyy HH:MM:SS"))
        if calc_swr > 0
            swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file,loc_time_agg)
        else
            for_tau, Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file,loc_time_agg)
        end
    else
         Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file)
    end

    if save_images
        SHIs, images = create_exmat(outdir,outstr,pts,g_img,append_file)
    end

    if progress
        elapsed = time() - start
        progtextinit = "0. Pre-calc took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
        if ispath(outdir*"/"*"ProgressLastPoint/")
            rm(outdir*"/"*"ProgressLastPoint/",recursive=true)
            mkdir(outdir*"/"*"ProgressLastPoint/")
        else
            mkdir(outdir*"/"*"ProgressLastPoint/")
        end
        writedlm(outdir*"/"*"ProgressLastPoint/"*progtextinit,NaN)
    end

    ###############################################################################
    # > Loop through the points

        @simd for crx = crxstart:size(pts_x,1)

            ######################
            #### calculate horizon lines

            if progress; start = time(); end

            if terrain
                if terrain_highres
                    pt_dtm_x, pt_dtm_y =  pcd2pol2cart(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(300),"terrain",ch,pts_slp[crx],dtm_cellsize);
                end

                if terrain_lowres && !terrain_tile
                    pt_dem_x, pt_dem_y = pcd2pol2cart(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e_dem[crx],terrain_peri,"terrain",ch,pts_slp[crx],dem_cellsize);
                end

                if terrain_highres && (terrain_lowres || horizon_line)
                    pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
                elseif terrain_highres
                    pt_dtm_x, pt_dtm_y = prepterdat(pt_dtm_x,pt_dtm_y)
                elseif terrain_lowres && !terrain_highres
                    pt_dtm_x, pt_dtm_y = prepterdat(pt_dem_x,pt_dem_y)
                end
            end

            if pt_corr
                pt_chm_x, pt_chm_y, pt_chm_r, pt_chm_x_thick, pt_chm_y_thick = calcCHM_Ptrans(copy(chm_x),copy(chm_y),copy(chm_z),
                                    copy(chm_b),copy(chm_lavd),pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,ch,chm_cellsize) # calculated points

                pt_chm_x_pts, pt_chm_y_pts, pt_chm_r_pts = pcd2pol2cart(copy(chm_x),copy(chm_y),copy(chm_z),pts_x[crx],pts_y[crx],
                                    pts_e[crx],surf_peri,"chm",ch,pts_slp[crx],chm_cellsize) # pts from the CHM

                if trunks_2
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
                                                            pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5),"surface",ch,pts_slp[crx],0);

                    else
                        pt_tsm_x, pt_tsm_y, _ = pcd2pol2cart(pt_tsm_x,pt_tsm_y,pt_tsm_z,
                                                            pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5),"surface",ch,pts_slp[crx],0);
                    end
                end

            else
                #  100% opaque canopy:
                pt_chm_x, pt_chm_y = pcd2pol2cart(copy(chm_x),copy(chm_y),copy(chm_z),pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,"terrain",ch,pts_slp[crx],chm_cellsize)
                if terrain # merge if using terrain
                    pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_chm_x,pt_dtm_x),append!(pt_chm_y,pt_dtm_y));
                else
                    pt_dtm_x, pt_dtm_y = prepterdat(pt_chm_x,pt_chm_y);
                end
            end

            if progress
                elapsed = time() - start
                if crx != 1; try; rm(outdir*"/ProgressLastPoint/"*progtext1); catch; end; end
                global progtext1 = "1. Transferring to polar took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                writedlm(outdir*"/ProgressLastPoint/"*progtext1,NaN)
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
                mat2ev = fillmat(kdtree,hcat(pt_chm_x_pts[pt_chm_r_pts .> 10],pt_chm_y_pts[pt_chm_r_pts .> 10]),4.0,kdtreedims,30,radius,mat2ev); # include canopy surface points
                if terrain
                    mat2ev = fillmat(kdtree,hcat(vcat(pt_chm_x_thick,pt_dtm_x),vcat(pt_chm_y_thick,pt_dtm_y)),1.5,kdtreedims,10,radius,mat2ev); # distance canopy is opaque and treated with terrain
                else
                    mat2ev = fillmat(kdtree,hcat(pt_chm_x_thick,pt_chm_y_thick),1.5,kdtreedims,10,radius,mat2ev); # distance canopy is opaque and treated with terrain
                end
            else # treat all canopy as opaque (like terrain) -> all points come in as terrain points from above section
                mat2ev = fillmat(kdtree,hcat(pt_dtm_x,pt_dtm_y),1.0,kdtreedims,10,radius,mat2ev); # use this line if plotting opaque canpoy
            end

            # add trunks
            if trunks_2
                mat2ev = fillmat(kdtree,hcat(pt_tsm_x,pt_tsm_y),2.0,kdtreedims,15,radius,mat2ev)
            end

            mat2ev[isnan.(g_rad)] .= 1;

            if progress
                elapsed = time() - start
                if crx != 1; try; rm(outdir*"/ProgressLastPoint/"*progtext2); catch; end; end
                global progtext2 = "2. Classifying image took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                writedlm(outdir*"/ProgressLastPoint/"*progtext2,NaN)
            end

            ########################
            ### Perform calculations

            if progress; start = time(); end

            ##### Calculate Vf
            Vf_w, Vf_f = calculateVf(mat2ev,g_rad,radius)

            ##### CalcSWR
            if calc_trans
                sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts_x[crx],pts_y[crx],loc_time,time_zone,coor_system,utm_zone)
                transfor = calc_transmissivity(float.(mat2ev),loc_time,tstep,radius,sol_phi,sol_tht,g_coorpol,0.0,0.0,
                                        im_centre,trans_for,lens_profile_tht,lens_profile_rpix)
                if calc_swr == 1
                    swrtot, swrdir, _ = calculateSWR(transfor,sol_sinelev,sol_tht,sol_phi,loc_time,max.(1367*sol_sinelev,0),Vf_w)
                elseif calc_swr == 2
                    swrtot, swrdir, _ = calculateSWR(transfor,sol_sinelev,sol_tht,sol_phi,loc_time,swr_open,Vf_w)
                end
            end

            if progress
                elapsed = time() - start
                if crx != 1; try; rm(outdir*"/ProgressLastPoint/"*progtext3); catch; end; end
                global progtext3 = "3. Vf and SWR calculations took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                writedlm(outdir*"/ProgressLastPoint/"*progtext3,NaN)
            end

            ######################
            #export the data [append to netcdf]
            if progress; start = time(); end

            Vf_weighted[crx] = np.array(Vf_w)
            Vf_flat[crx]     = np.array(Vf_f)

            if calc_trans
                for_tau[crx] = np.array(vec(aggregate_data(loc_time,loc_time_agg,transfor,tstep)))
                if calc_swr > 0
                    swr_tot[crx] = np.array(vec(aggregate_data(loc_time,loc_time_agg,swrtot,tstep)))
                    swr_dir[crx] = np.array(vec(aggregate_data(loc_time,loc_time_agg,swrdir,tstep)))
                end
            end

            if save_images
                SHIs[crx] = np.array(transpose(mat2ev))
            end

            if progress
                elapsed = time() - start
                if crx != 1; try; rm(outdir*"/ProgressLastPoint/"*progtext4); catch; end; end
                global progtext4 = "4. Exporting data took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                writedlm(outdir*"/ProgressLastPoint/"*progtext4,NaN)
            end


            # save the progress
            percentdone = Int(floor((crx / size(pts,1)) * 100))
            if crx == 1
                # delete any other processing text files in the output folder
                for f in readdir(outdir)
                        if startswith.(f,"Processing")
                                rm(outdir*"/"*f)
                        end
                end
            else
                rm(outdir*"/"*outtext)
            end
            global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
            writedlm(outdir*"/"*outtext,NaN)

        end # end crx

        dataset.close()
        if save_images; images.close(); end

        println("done with "*taskID)

end
