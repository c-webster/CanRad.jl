function LAS2Rad(pts,dat_in,par_in,exdir,taskID="task")

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
    # > Import surface/terrain data and clip

    dsm_x, dsm_y, dsm_z = readlas(dsmf)

    # clip dsm within eval-peri of min/max pts
    limits = hcat((floor(minimum(pts_x))-surf_peri),(ceil(maximum(pts_x))+surf_peri),
                    (floor(minimum(pts_y))-surf_peri),(ceil(maximum(pts_y))+surf_peri))

    dsm_x, dsm_y, dsm_z, _ = clipdat(dsm_x,dsm_y,dsm_z,limits,surf_peri)

    # import, clip and prepare terrain daa
    if terrain_highres || terrain_lowres || horizon_line
        terrain = true
    else
        terrain = false
    end

    # load and clip the dtm
    if terrain_highres

        limits = hcat((floor(minimum(pts_x))-surf_peri*2),(ceil(maximum(pts_x))+surf_peri*2),
                        (floor(minimum(pts_y))-surf_peri*2),(ceil(maximum(pts_y))+surf_peri*2))

        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits,true, true)

    end

    # load and clip the dem
    if terrain_lowres

        limits = hcat((floor(minimum(pts_x))-terrain_peri),(ceil(maximum(pts_x))+terrain_peri),
                        (floor(minimum(pts_y))-terrain_peri),(ceil(maximum(pts_y))+terrain_peri))

        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits,true,true)
        pts_e_dem = findelev(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,100)

    end

    # load the buildings
    if buildings
        rsm_x, rsm_y, rsm_z, rsm_cellsize = read_griddata(rsmf,true,true)
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

    # > Calculate slope and aspect from dtm data
    if tilt # disabled in Preparatory_Functions.jl
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
        pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x))
        imcX = 0.0; imcY = 0.0
    end

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

    # create the empty matrix
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
        if ispath(outdir*"/ProgressLastPoint/")
            rm(outdir*"/ProgressLastPoint/",recursive=true)
            mkdir(outdir*"/ProgressLastPoint/")
        else
            mkdir(outdir*"/ProgressLastPoint/")
        end
        writedlm(outdir*"/ProgressLastPoint/"*progtextinit,NaN)
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

                if buildings
                    pt_rsm_x, pt_rsm_y =  pcd2pol2cart(copy(rsm_x),copy(rsm_y),copy(rsm_z),pts_x[crx],pts_y[crx],pts_e[crx],Int.(50),"terrain",ch,pts_slp[crx],rsm_cellsize);
                    if tershad < 3
                        pt_dtm_x, pt_dtm_y = prepterdat(append!(pt_dtm_x,pt_rsm_x),append!(pt_dtm_y,pt_rsm_y));
                    else
                        pt_dtm_x, pt_dtm_y = prepterdat(pt_rsm_x,pt_rsm_y)
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

                # add terrain
                if terrain
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

                Vf_weighted[crx] = Int(round(Vf_w*100))
                Vf_flat[crx]     = Int(round(Vf_f*100))

                if calc_trans
                    for_tau[:,crx] = Int.(round.((vec(aggregate_data(loc_time,loc_time_agg,transfor,tstep))).*100))
                    if calc_swr > 0
                        if calc_swr == 2
                            swr_tot[isnan.(swr_tot)] .= -9999
                            swr_dir[isnan.(swr_dir)] .= -9999
                        end
                        swr_tot[:,crx] = Int.(round.(vec(aggregate_data(loc_time,loc_time_agg,swrtot,tstep))))
                        swr_dir[:,crx] = Int.(round.(vec(aggregate_data(loc_time,loc_time_agg,swrdir,tstep))))
                    end
                end

                if save_images
                    SHIs[:,:,crx] = mat2ev
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

        close(dataset)
        if save_images; close(images); end

        println("done with "*taskID)

end
