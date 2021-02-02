function TLS2Rad(pts,dat_in,par_in,exdir,taskID="task")


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
        outstr = String(split(exdir,"/")[end-1])
        global outdir = exdir
    end

    if !ispath(outdir)
        mkpath(outdir)
    end

    crxstart = 1; append_file  = false

    ###############################################################################
    # > Preprocess DSM/DTM

    # load the las
    dsm_x, dsm_y, dsm_z = readlas(dsmf)

    # load the dtm
    dtm_x, dtm_y, dtm_z, dtm_cellsize = importdtm(dtmf,tilt)

    # clip dsm within eval-peri of min/max pts
    dsm_x, dsm_y, dsm_z, _ = clipdat(dsm_x,dsm_y,dsm_z,Int.(floor.(pts)),surf_peri)
    limits = Int.(floor.(hcat(vcat(minimum(dsm_x),maximum(dsm_x)),vcat(minimum(dsm_y),maximum(dsm_y)))))
    dtm_x, dtm_y, dtm_z, _ = clipdat(dtm_x,dtm_y,dtm_z,limits,surf_peri)

    pts_e     = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y)

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
        loc_time = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
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

        rbins = collect(0:(surf_peri-0)/10:surf_peri)
        tol   = tolerance.*collect(reverse(0.5:0.05:1))

        if calc_trans
            drad, im_centre, lens_profile_tht, lens_profile_rpix, trans_for = get_constants(g_img,loc_time)
        end

        @simd for crx = 1:size(pts_x,1)

            if progress; start = time(); end

            pt_dsm_x, pt_dsm_y, pt_dsm_z = getsurfdat(dsm_x,dsm_y,dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri);

            pt_dsm_x, pt_dsm_y, pt_dsm_z = pcd2pol2cart(pt_dsm_x,pt_dsm_y,pt_dsm_z,pts_x[crx],pts_y[crx],pts_e[crx],surf_peri,"surface");


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
            if crx != 1; rm(outdir*"/"*outtext); end
            global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
            writedlm(outdir*"/"*outtext,NaN)


        end #end crx

        dataset.close()
        if save_images; images.close(); end

        println("done with: "*taskID)

    # catch
    #     dataset.close()
    #     if save_images; images.close(); end
    # end
end
