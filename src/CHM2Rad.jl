function CHM2Rad(pts,dat_in,par_in,exdir,taskID="task")

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

    crxstart = 1; append  = false

    ################################################################################

    # Import terrain/surface data and clip
    chm_x, chm_y, chm_z, chm_cellsize = read_ascii(chmf)

    _, _, chm_lavd, _ = read_ascii(lavdf)

    dtm_x, dtm_y, dtm_z, dtm_cellsize = importdtm(dtmf,tilt)

    dem_x, dem_y, dem_z, dem_cellsize = read_ascii(demf)

    limits = Int.(floor.(hcat(vcat(minimum(pts_x),maximum(pts_x)),vcat(minimum(pts_y),maximum(pts_y)))))
     chm_x, chm_y, chm_z, _ = clipdat(chm_x,chm_y,chm_z,limits,surf_peri)
    dtm_x, dtm_y, dtm_z, _ = clipdat(dtm_x,dtm_y,dtm_z,limits,surf_peri*4)
    dem_x, dem_y, dem_z, _ = clipdat(dem_x,dem_y,dem_z,limits,terrain_peri)

    pts_e     = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y)
    pts_e_dem = findelev(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,100)
    chm_z     = chm_z + (findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),chm_x,chm_y))

    ###############################################################################
    # > Calculate slope and aspect from dtm data
    if tilt
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
        pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x))
    end

    # create an empty matrix for Vf calculation
    g_rad, g_coorpol, g_coorcrt, g_img = create_mat(radius)
    g_coorcrt = ((g_coorcrt .- radius) ./ radius) .* 90
    g_img[isnan.(g_rad)] .= 1

    # make g_coorcrt a KDtree for easy look up
    kdtree = scipyspat.cKDTree(g_coorcrt)
    kdtreedims = size(g_coorcrt,1)

    # create the output files
    if calc_trans
        loc_time = collect(Dates.DateTime(t1,"yyyy.mm.dd HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"yyyy.mm.dd HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t1,"yyyy.mm.dd HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"yyyy.mm.dd HH:MM:SS")-Dates.Minute(tstep))
        for_tau, Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,calc_trans,false,loc_time_agg)
    else
         Vf_weighted, Vf_flat, dataset = createfiles(outdir,outstr,pts,calc_trans,false)
    end
    if save_images
        SHIs, images = create_exmat(outdir,outstr,pts,g_img,append)
    end

    if tershad < 2 && terrain == 2
        pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,mean.(pts_x),mean.(pts_y),mean.(pts_e_dem),terrain_peri,"terrain",ch,0.0,dem_cellsize);
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

    rbins = collect(0:(surf_peri-0)/5:surf_peri)
    tol   = tolerance.*collect(reverse(0.75:0.05:1))

    if calc_trans
        drad, im_centre, lens_profile_tht, lens_profile_rpix, trans_for = get_constants(g_img,loc_time)
    end

    try
        @simd for crx = crxstart:size(pts_x,1)

            ######################
            #### calculate horizon lines

            if progress; start = time(); end

            # pt_chm_x, pt_chm_y = pcd2pol2cart(copy(chm_x),copy(chm_y),copy(chm_z),pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,"terrain",ch,pts_slp[crx],dtm_cellsize);
            pt_chm_x, pt_chm_y, pt_chm_r, pt_chm_x_thick, pt_chm_y_thick = calcCHM_Ptrans(copy(chm_x),copy(chm_y),copy(chm_z),chm_lavd,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,ch,chm_cellsize)

            pt_chm_x_pts, pt_chm_y_pts, pt_chm_r_pts = pcd2pol2cart(copy(chm_x),copy(chm_y),copy(chm_z),pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,"chm",ch,pts_slp[crx],chm_cellsize)
            pt_dtm_x, pt_dtm_y = pcd2pol2cart(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(3*dem_cellsize),"terrain",ch,pts_slp[crx],dtm_cellsize);

            if terrain == 1
                pt_dem_x, pt_dem_y = pcd2pol2cart(dem_x,dem_y,dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],terrain_peri,"terrain",ch,pts_slp[crx],dem_cellsize);
            end

            # prep the data for occupying the matrix
            pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));

            elapsed = time() - start

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

            #occupy matrix with dtm
            for zdx = 1:1:size(rbins,1)-1
                # global mat2ev
                ridx   = findall(rbins[zdx] .<= pt_chm_r .< rbins[zdx+1])
                ridx2  = findall(rbins[zdx] .<= pt_chm_r_pts .< rbins[zdx+1])
                mat2ev = fillmat(kdtree,hcat(vcat(pt_chm_x[ridx],pt_chm_x_pts[ridx2]),vcat(pt_chm_y[ridx],pt_chm_y_pts[ridx2])),tol[zdx],kdtreedims,50,radius,mat2ev);
            end
            # mat2ev = fillmat(kdtree,hcat(vcat(pt_chm_x,pt_chm_x_pts),vcat(pt_chm_y,pt_chm_y_pts)),0.5,kdtreedims,50,radius,mat2ev);
            mat2ev = fillmat(kdtree,hcat(vcat(pt_chm_x_thick,pt_dtm_x),vcat(pt_chm_y_thick,pt_dtm_y)),1.0,kdtreedims,10,radius,mat2ev);
            # mat2ev = fillmat(kdtree,hcat(pt_dtm_x,pt_dtm_y),1.0,kdtreedims,10,radius,mat2ev);
            # mat2ev = fillmat(kdtree,hcat(pt_chm_x,pt_chm_y),1.0,kdtreedims,10,radius,mat2ev);

            mat2ev[isnan.(g_rad)] .= 1;
            writedlm("C:/Users/webster/GDrive/Temp/LAS2Rad/temp.txt",mat2ev)

            elapsed = time() - start

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

            ##### CalcSWR to be added later
            if calc_trans
                sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts_x[crx],pts_y[crx],loc_time,time_zone,coor_system,utm_zone)
                transfor = calc_transmissivity(float.(mat2ev),loc_time,tstep,radius,sol_phi,sol_tht,g_coorpol,0.0,0.0,drad,
                                        im_centre,trans_for,lens_profile_tht,lens_profile_rpix)
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
            if crx != 1; rm(outdir*"/"*outtext); end
            global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
            writedlm(outdir*"/"*outtext,NaN)

        end

        dataset.close()
        if save_images; images.close(); end

        println("done with: "*taskID)

    catch
        dataset.close()
        if save_images; images.close(); end
        println(taskID*" failed")

    end

end
