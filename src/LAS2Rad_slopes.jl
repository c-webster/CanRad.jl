function LAS2Rad_slopes(pts,dat_in,par_in,exdir,taskID="task")


    ################################################################################
    # Initialise
    eval(extract(dat_in))
    eval(extract(par_in))

    pts = float.(pts)
    pts_x = float(pts[:,1])
    pts_y = float(pts[:,2])

    if progress; start = time(); end

    if batch
        outstr = string(Int(floor(pts[1,1])))*"_"*string(Int(floor(pts[1,2])))
        global outdir = exdir*"/"*outstr
    else
        outstr = "Output_"*String(split(exdir,"/")[end-1])
        global outdir = exdir
    end
    if !ispath(outdir)
        mkpath(outdir)
    end

    crxstart = 1; append_file  = false

    ################################################################################

    # load the las
    dsm_x, dsm_y, dsm_z = readlas(dsmf)

    # load the dtm
    if tershad < 3
        if tilt
            dtm_x, dtm_y, dtm_z, dtm_s, dtm_a, dtm_cellsize = importdtm(dtmf,tilt)
        else
            dtm_x, dtm_y, dtm_z, dtm_cellsize = importdtm(dtmf,tilt)
        end

    # load the dem
        if tershad < 2
            dem_x, dem_y, dem_z, dem_cellsize = read_ascii(demf)
        end
    end
    ###############################################################################
    # > Preprocess DSM/DTM

    # clip dsm within eval-peri of min/max pts
    dsm_x, dsm_y, dsm_z, _ = clipdat(dsm_x,dsm_y,dsm_z,Int.(floor.(pts)),surf_peri)
    limits = Int.(floor.(hcat(vcat(minimum(dsm_x),maximum(dsm_x)),vcat(minimum(dsm_y),maximum(dsm_y)))))
    dtm_x, dtm_y, dtm_z, _ = clipdat(dtm_x,dtm_y,dtm_z,limits,surf_peri)
    if tershad < 2; dem_x, dem_y, dem_z, _ = clipdat(dem_x,dem_y,dem_z,limits,terrain_peri); end

    # determine ground elevation of points
    pts_e     = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y)
    if tershad < 2; pts_e_dem = findelev(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,100); end

    ###############################################################################
    # > Generate extra canopy elements

    # load the dbh
    if trunks
        dbh_x, dbh_y, dbh_z, dbh_r = loaddbh(dbhf,limits)
        dbh_e = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),dbh_x,dbh_y)

        tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,30,0.1,dbh_e)
    end

    # load the ltc
    if branches
        ltc = loadltc(ltcf,limits)

        bsm_x, bsm_y, bsm_z = make_branches(ltc)
        bsm_z .+= findelev((dtm_x),(dtm_y),(dtm_z),bsm_x,bsm_y)
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
    loc_time = collect(Dates.DateTime(t1,"yyyy.mm.dd HH:MM:SS"):Dates.Minute(int):Dates.DateTime(t2,"yyyy.mm.dd HH:MM:SS"))

        swr_tot_s, swr_dir_s, for_tau_s, _, _, dataset_s = createfiles(outdir,outstr*"_SOUTH",pts,calc_trans,false,loc_time)
        swr_tot_w, swr_dir_w, for_tau_w, _, _, dataset_w = createfiles(outdir,outstr*"_WEST",pts,calc_trans,false,loc_time)
        swr_tot_n, swr_dir_n, for_tau_n, _, _, dataset_n = createfiles(outdir,outstr*"_NORTH",pts,calc_trans,false,loc_time)
        swr_tot_e, swr_dir_e, for_tau_e, _, _, dataset_e = createfiles(outdir,outstr*"_EAST",pts,calc_trans,false,loc_time)

        Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,calc_trans,false)

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
    # try

        rbins = collect(0:(surf_peri-0)/5:surf_peri)
        tol   = tolerance.*collect(reverse(0.5:0.1:1))

        @simd for crx = crxstart:size(pts_x,1)

            # try # catch statement to close dataset if error

                if progress; start = time(); end


                #### transfer point clouds to polar coordinates
                if branches
                    pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_x,dsm_y,dsm_z,bsm_x,bsm_y,bsm_z,pts_x[crx],pts_y[crx],surf_peri)
                else
                    pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_x,dsm_y,dsm_z,pts_x[crx],pts_y[crx],surf_peri);
                end
                pt_dsm_x, pt_dsm_y, pt_dsm_z = pcd2pol2cart(pt_dsm_x,pt_dsm_y,pt_dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,"surface",ch,pts_slp[crx],0);

                if trunks
                    pt_tsm_x, pt_tsm_y, pt_tsm_z = getsurfdat(tsm_x,tsm_y,tsm_z,pts_x[crx],pts_y[crx],Int.(surf_peri*0.3))
                    tidx = findall(dist(dbh_x,dbh_y,pts[crx,1],pts[crx,2]) .< 3)
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

                if tershad < 3
                    pt_dtm_x, pt_dtm_y =  pcd2pol2cart(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(300),"terrain",ch,pts_slp[crx],dtm_cellsize);
                    if tershad < 2 && terrain == 1
                        pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],terrain_peri,"terrain",ch,pts_slp[crx],dem_cellsize);
                        pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
                    else
                        pt_dtm_x, pt_dtm_y = prepterdat(pt_dtm_x,pt_dtm_y)
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
                if trunks
                    mat2ev = fillmat(kdtree,hcat(pt_tsm_x,pt_tsm_y),2.0,kdtreedims,15,radius,mat2ev)
                end

                #occupy matrix with dtm
                if tershad < 3
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
                sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts_x[crx],pts_y[crx],loc_time,time_zone,coor_system,utm_zone)
                swrtot_s, swrdir_s, _, transfor_s = calculateSWR(float.(mat2ev),loc_time,radius,sol_phi,sol_tht,Vf_w,g_coorpol,sol_sinelev,imcX,imcY)
                swrtot_w, swrdir_w, _, transfor_w = calculateSWR(float.(mat2ev),loc_time,radius,sol_phi.-270,sol_tht,Vf_w,g_coorpol,sol_sinelev,imcX,imcY)
                swrtot_n, swrdir_n, _, transfor_n = calculateSWR(float.(mat2ev),loc_time,radius,sol_phi.-180,sol_tht,Vf_w,g_coorpol,sol_sinelev,imcX,imcY)
                swrtot_e, swrdir_e, _, transfor_e = calculateSWR(float.(mat2ev),loc_time,radius,sol_phi.-90,sol_tht,Vf_w,g_coorpol,sol_sinelev,imcX,imcY)

                if progress
                    elapsed = time() - start
                    if crx != 1; try; rm(outdir*"/ProgressLastPoint/"*progtext3); catch; end; end
                    global progtext3 = "3. Vf and SWR calculations took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                    writedlm(outdir*"/ProgressLastPoint/"*progtext3,NaN)
                end

                ######################
                #export the data [append to netcdf]
                if progress; start = time(); end

                if ~tilt & calc_swr
                    swr_tot_s[crx] = np.array(swrtot_s); swr_dir_s[crx] = np.array(swrdir_s); for_tau_s[crx] = np.array(transfor_s)
                    swr_tot_w[crx] = np.array(swrtot_w); swr_dir_w[crx] = np.array(swrdir_w); for_tau_w[crx] = np.array(transfor_w)
                    swr_tot_n[crx] = np.array(swrtot_n); swr_dir_n[crx] = np.array(swrdir_n); for_tau_n[crx] = np.array(transfor_n)
                    swr_tot_e[crx] = np.array(swrtot_e); swr_dir_e[crx] = np.array(swrdir_e); for_tau_e[crx] = np.array(transfor_e)
                end
                Vf_weighted[crx] = np.array(Vf_w)
                Vf_flat[crx]     = np.array(Vf_f)

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
                if crx != 1; rm(outdir*"/"*outtext); end
                global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
                writedlm(outdir*"/"*outtext,NaN)

            # catch
            #     writedlm(outdir*"/Error_"*string(Int(pts[crx,1]))*"_"*string(Int(pts[crx,2]))*".txt",NaN)
            # end

        end #end crx

        dataset.close()
        dataset_s.close()
        dataset_w.close()
        dataset_n.close()
        dataset_e.close()
        if save_images; images.close(); end

        println("done with: "*taskID)

    # catch
    #     dataset.close()
    #     if save_images; images.close(); end
    # end
end