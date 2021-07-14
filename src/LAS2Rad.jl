function LAS2Rad(pts,dat_in,par_in,exdir,taskID="task")


    ################################################################################
    # Initialise

    dat_in, par_in = compatability_check(dat_in,par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    pts = float.(pts)
    pts_x = float(pts[:,1])
    pts_y = float(pts[:,2])

    if progress; start = time(); end

    if batch
        # outstr = string(Int(floor(pts[1,1])))*"_"*string(Int(floor(pts[1,2])))
        outstr = split(taskID,"_")[2]*"_"*split(taskID,"_")[3][1:end-4]
        global outdir = exdir*"/"*outstr
    else
        outstr = String(split(exdir,"/")[end-1])
        global outdir = exdir
    end
    if !ispath(outdir)
        mkpath(outdir)
    end

    # writedlm(outdir*"/"*basename(taskID),taskID)
    # if size(readdir(outdir),1) == 5
    #     crxstart = parse(Int,(split(reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])))[4])+1
    #     global outtext  = reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])
    #     append   = true
    # else; crxstart = 1; append_file  = false
    # end

    crxstart = 1; append_file  = false

    ################################################################################

    # load the las
    dsm_x, dsm_y, dsm_z = readlas(dsmf)

    # load the dtm
    if tershad < 3
        if tilt; dtm_x, dtm_y, dtm_z, dtm_s, dtm_a, dtm_cellsize = importdtm(dtmf,tilt)
        else; dtm_x, dtm_y, dtm_z, dtm_cellsize = importdtm(dtmf,tilt); end
    end

    # load the dem
    if tershad < 2
        dem_x, dem_y, dem_z, dem_cellsize = read_ascii(demf)
    end

    # load the buildings
    if buildings
        rsm_x, rsm_y, rsm_z, rsm_cellsize = read_ascii(rsmf)
    end

    # Load SWR data
    # if calc_trans && calc_swr == 2
    #     swrdat   = readdlm(swrf)
    #     tstep    = Dates.value(Dates.Minute(Dates.DateTime(swrdat[3,1].*" ".*swrdat[3,2],"dd.mm.yyyy HH:MM:SS")-Dates.DateTime(swrdat[2,1].*" ".*swrdat[2,2],"dd.mm.yyyy HH:MM:SS")))
    #     t1       = swrdat[1,1].*" ".*swrdat[1,2]
    #     t2       = swrdat[end,1].*" ".*swrdat[end,2]
    #     swr_open = float.(swrdat[:,3])
    # end

    ###############################################################################
    # > Preprocess DSM/DTM

    # clip dsm within eval-peri of min/max pts
    dsm_x, dsm_y, dsm_z, _ = clipdat(dsm_x,dsm_y,dsm_z,pts[:,1:2],surf_peri)
    limits = Int.(floor.(hcat(vcat(minimum(dsm_x),maximum(dsm_x)),vcat(minimum(dsm_y),maximum(dsm_y)))))

    if tershad < 3; dtm_x, dtm_y, dtm_z, _ = clipdat(dtm_x,dtm_y,dtm_z,limits,surf_peri*4); end
    if tershad < 2; dem_x, dem_y, dem_z, _ = clipdat(dem_x,dem_y,dem_z,limits,terrain_peri); end

    # determine ground elevation of points
    if tershad == 3;
        if size(pts,2) > 2; pts_e = pts[:,3]; else;  pts_e = fill(0.0,size(pts_x)); end
    else; pts_e = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y); end

    if tershad < 2; pts_e_dem = findelev(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,100); end

    ###############################################################################
    # > Generate extra canopy elements

    # load the dbh
    if trunks
        dbh_x, dbh_y, dbh_z, dbh_r, lastc = loaddbh(dbhf,limits,-50)
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

    # load the ltc
    if branches
        if extension(ltcf) == ".txt"
            ltc = loadltc_txt(ltcf,limits,0)
        elseif extension(ltcf) == ".laz"
            ltc = loadltc_laz(ltcf,limits,-50,dbh_x,dbh_y,dbh_e,lastc)
            ltc[:,3] = ltc[:,3] .- findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),ltc[:,1],ltc[:,2])
            ltc = ltc[setdiff(1:end, findall(ltc[:,3].<1)), :]
        end

        bsm_x, bsm_y, bsm_z = make_branches(ltc,b_space)
        bsm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),bsm_x,bsm_y)
    end

    # Load SWR data
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
    # > Calculate slope and aspect from dtm data
    if tilt
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
        pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x))
        imcX = 0.0; imcY = 0.0
    end

    ###############################################################################
    # > Start creating the 2D canopy/terrain map

    # create the empty matrix
    g_rad, g_coorpol, g_coorcrt, g_img = create_mat(radius)
    g_coorcrt = ((g_coorcrt .- radius) ./ radius) .* 90
    g_img[isnan.(g_rad)] .= 1

    # make g_coorcrt a KDtree for easy look up
    if ~tilt
        kdtree = scipyspat.cKDTree(g_coorcrt)
        kdtreedims = size(g_coorcrt,1)
    end

    # generate the output files
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

    # calculate horizon matrix for grid-cell
    if tershad < 2 && terrain == 2
        pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,mean(pts_x),mean(pts_y),mean(pts_e_dem),terrain_peri,"terrain",ch,0.0,dem_cellsize);
    end


    if progress
        elapsed = time() - start
        progtextinit = "0. Pre-calc took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
        if ispath(outdir*"/ProgressLastPoint/")
            rm(outdir*"/ProgressLastPoint/",recursive=true)
            mkdir(outdir*"/ProgressLastPoint/")
        else
            mkdir(outdir*"/ProgressLastPoint/")
        end
        writedlm(outdir*"/ProgressLastPoint/"*progtextinit,NaN)
    end

    ###############################################################################
    # > Loop through the points
    # try

        rbins = collect(0:(surf_peri-0)/5:surf_peri)
        tol   = tolerance.*collect(reverse(0.75:0.05:1))

        if calc_trans
            drad, im_centre, lens_profile_tht, lens_profile_rpix, trans_for = get_constants(g_img,loc_time)
        end

        @simd for crx = crxstart:size(pts_x,1)

            # try # catch statement to close dataset if error

                if progress; start = time(); end


                #### transfer point clouds to polar coordinates
                if branches
                    pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_x,dsm_y,dsm_z,bsm_x,bsm_y,bsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri)
                else
                    pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_x,dsm_y,dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri);
                end
                pt_dsm_x, pt_dsm_y, pt_dsm_z = pcd2pol2cart(pt_dsm_x,pt_dsm_y,pt_dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,"surface",ch,pts_slp[crx],0);

                if trunks_2
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
                                                            pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5),"surface",ch,pts_slp[crx],0);

                    else
                        pt_tsm_x, pt_tsm_y, _ = pcd2pol2cart(pt_tsm_x,pt_tsm_y,pt_tsm_z,
                                                            pts_x[crx],pts_y[crx],pts_e[crx],Int.(surf_peri*0.5),"surface",ch,pts_slp[crx],0);
                    end
                end

                if tershad < 3
                    pt_dtm_x, pt_dtm_y =  pcd2pol2cart(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(300),"terrain",ch,pts_slp[crx],dtm_cellsize);
                    if tershad < 2
                        if terrain == 1
                        pt_dem_x, pt_dem_y = pcd2pol2cart(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e_dem[crx],terrain_peri,"terrain",ch,pts_slp[crx],dem_cellsize);
                        end
                        pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
                    else
                        pt_dtm_x, pt_dtm_y = prepterdat(pt_dtm_x,pt_dtm_y)
                    end
                end

                if buildings
                    pt_rsm_x, pt_rsm_y =  pcd2pol2cart(copy(rsm_x),copy(rsm_y),copy(rsm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(50),"terrain",ch,pts_slp[crx],rsm_cellsize);
                    if tershad < 3
                        pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_rsm_x),append!(pt_dtm_y,pt_rsm_y));
                    else
                        pt_dtm_x, pt_dtm_y = prepterdat(pt_rsm_x,pt_rsm_y)
                    end
                end

                if tilt
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
                    if crx != 1; try; rm(outdir*"/ProgressLastPoint/"*progtext1); catch; end; end
                    global progtext1 = "1. Transferring to polar took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                    writedlm(outdir*"/ProgressLastPoint/"*progtext1,NaN)
                end

                ######################
                ### Classify image
                if progress; start = time(); end

                mat2ev  = copy(g_img);

                # occupy matrix with surface points
                for zdx = 1:1:size(rbins,1)-1
                    ridx = findall(rbins[zdx] .<= pt_dsm_z .< rbins[zdx+1])
                    imdx = reshape(findpairs(kdtree,hcat(pt_dsm_x[ridx],pt_dsm_y[ridx]),tol[zdx],kdtreedims,40),(radius*2,radius*2))
                    mat2ev[imdx.==1] .= 0
                end

                # add trunks
                if trunks_2
                    mat2ev = fillmat(kdtree,hcat(pt_tsm_x,pt_tsm_y),2.0,kdtreedims,15,radius,mat2ev)
                end

                #occupy matrix with dtm
                if @isdefined(pt_dtm_x)
                    mat2ev = fillmat(kdtree,hcat(pt_dtm_x,pt_dtm_y),1.0,kdtreedims,10,radius,mat2ev);
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

                ##### Calculate SWR/forest transmissivity
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

                percentdone = Int(floor((crx / size(pts,1)) * 100))
                if crx == 1
                    # delete any other processing text files in the output folder
                    for f in readdir(outdir)
                            if startswith.(f,"Processing") && !startswith.(f,"Processing 100");
                                    rm(outdir*"/"*f)
                            end
                    end
                else
                    rm(outdir*"/"*outtext)
                end
                if crx != 1; rm(outdir*"/"*outtext); end
                global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
                writedlm(outdir*"/"*outtext,NaN)

            # catch
            #     writedlm(outdir*"/Error_"*string(Int(pts[crx,1]))*"_"*string(Int(pts[crx,2]))*".txt",NaN)
            # end

        end #end crx

        dataset.close()
        if save_images; images.close(); end

        println("done with "*taskID)

    # catch
    #     dataset.close()
    #     if save_images; images.close(); end
    # end
end
