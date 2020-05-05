function LAS2Rad(pts,dat_in,par_in,exdir,taskID="task")


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
        outstr = String(split(exdir,"/")[end])
        global outdir = exdir
    end
    if !ispath(outdir)
        mkpath(outdir)
    end
    # writedlm(outdir*"/"*basename(taskID),taskID)
    if size(readdir(outdir),1) == 5
        crxstart = parse(Int,(split(reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])))[4])+1
        global outtext  = reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])
        append   = true
    else; crxstart = 1; append  = false
    end

    ################################################################################

    # load the las
    dsm_x, dsm_y, dsm_z = readlas(dsmf)

    # load the dtm
    if tilt
        dtm_x, dtm_y, dtm_z, dtm_s, dtm_a, dtm_cellsize = importdtm(dtmf,tilt)
    else
        dtm_x, dtm_y, dtm_z, dtm_cellsize = importdtm(dtmf,tilt)
    end

    # load the dem
     dem_x, dem_y, dem_z, dem_cellsize = read_ascii(demf)

    ###############################################################################
    # > Preprocess DSM/DTM

    # clip dsm within eval-peri of min/max pts
    dsm_x, dsm_y, dsm_z, _ = clipdat(dsm_x,dsm_y,dsm_z,Int.(floor.(pts)),surf_peri)
    bounds = Int.(floor.(hcat(vcat(minimum(dsm_x),maximum(dsm_x)),vcat(minimum(dsm_y),maximum(dsm_y)))))
    dtm_x, dtm_y, dtm_z, _ = clipdat(dtm_x,dtm_y,dtm_z,bounds,surf_peri*4)
    dem_x, dem_y, dem_z, _ = clipdat(dem_x,dem_y,dem_z,bounds,terrain_peri)

    # determine ground elevation of points
    pts_e     = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y)
    pts_e_dem = findelev(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,100)
    # dsm[:,3] .-= findelev(dtm,dsm[:,1],dsm[:,2])
    # dsm[:,3] .-=  dsm[:,4]

    ###############################################################################
    # > Generate extra canopy elements

    # load the dbh
    if trunks
        dbh_x, dbh_y, dbh_z, dbh_r = loaddbh(dbhf,bounds)
        dbh_e = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),dbh_x,dbh_y)

        tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,30,0.1,dbh_e)
    end

    # load the ltc
    if branches
        ltc = loadltc(ltcf,bounds)

        bsm_x, bsm_y, bsm_z = make_branches(ltc)
        bsm_z .+= findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),bsm_x,bsm_y)
    end

    ###############################################################################
    # > Calculate slope and aspect from dtm data
    if tilt
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
        pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x))
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

    loc_time = collect(Dates.DateTime(t1,"yyyy.mm.dd HH:MM:SS"):Dates.Minute(int):
                                Dates.DateTime(t2,"yyyy.mm.dd HH:MM:SS"))

    # generate the output files

    if calc_swr
        swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,loc_time,t1,t2,int,calc_swr,append)
    else
         Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,loc_time,t1,t2,int,calc_swr,append)
    end

    if save_images
        SHIs, images = create_exmat(outdir,outstr,pts,g_img,append)
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
    try

        rbins = collect(0:(surf_peri-0)/5:surf_peri)
        tol   = tolerance.*collect(reverse(0.75:0.05:1))

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
                    pt_tsm_x, pt_tsm_y, pt_tsm_z = getsurfdat(tsm_x,tsm_y,tsm_z,pts_x[crx],pts_y[crx],Int.(surf_peri*0.5))
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

                if tershad < 3
                    pt_dtm_x, pt_dtm_y =  pcd2pol2cart(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(3*dem_cellsize),"terrain",ch,pts_slp[crx],dtm_cellsize);
                    if tershad < 2 && terrain == 1
                        pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],terrain_peri,"terrain",ch,pts_slp[crx],dem_cellsize);
                        pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
                    end
                end

                if tilt
                    imcX, imcY = getimagecentre(pts_slp[crx],pts_asp[crx])
                    gcrcrt = zeros(size(g_coorcrt)) .* NaN
                    gcrcrt[:,1] = g_coorcrt[:,1] .+ imcX
                    gcrcrt[:,2] = g_coorcrt[:,2] .+ imcY
                    kdtree = scipyspat.cKDTree(gcrcrt)
                    kdtreedims = size(gcrcrt,1)
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
                if ~tilt & calc_swr
                    sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts_x[crx],pts_y[crx],loc_time,time_zone,coor_system,utm_zone)
                    swrtot, swrdir, swr_dif, transfor = calculateSWR(float.(mat2ev),loc_time,radius,sol_phi,
                                                                sol_tht,Vf_w,g_coorpol,sol_sinelev)
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

                if ~tilt & calc_swr
                    swr_tot[crx]     = np.array(swrtot)
                    swr_dir[crx]     = np.array(swrdir)
                    for_tau[crx]     = np.array(transfor)
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

    catch
        dataset.close()
        if save_images; images.close(); end

    finally
        dataset.close()
        if save_images; images.close(); end

        println("done with: "*taskID)
    end
end
