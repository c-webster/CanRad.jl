function LAS2Rad(pts,dat_in,par_in,exdir,taskID="task")


    ################################################################################
    # Initialise
    eval(extract(dat_in))
    eval(extract(par_in))

    pts = float.(pts)

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
        crxstart = parse(Int,(split(reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])))[4])
        append   = true
    else; crxstart = 1; append  = false
    end

    ################################################################################

    # load the las
    dsm = readlas(dsmf);

    # load the dtm
    dtm, dtm_cellsize = importdtm(dtmf,tilt)

    # load the dem
    dem, dem_cellsize = read_ascii(demf)

    ###############################################################################
    # > Preprocess DSM/DTM

    # clip dsm within eval-peri of min/max pts
    dsm = clipdat(dsm,Int.(floor.(pts)),surf_peri)
    bounds = Int.(floor.(hcat(vcat(minimum(dsm[:,1]),maximum(dsm[:,1])),vcat(minimum(dsm[:,2]),maximum(dsm[:,2])))))
    dtm = clipdat(dtm,bounds,surf_peri*4)
    dem = clipdat(dem,bounds,terrain_peri)

    # determine ground elevation of points
    pts_e    = findelev(dtm,pts[:,1],pts[:,2])
    # dsm[:,3] .-= findelev(dtm,dsm[:,1],dsm[:,2])
    # dsm[:,3] .-=  dsm[:,4]

    ###############################################################################
    # > Generate extra canopy elements

    # load the dbh
    if trunks
        dbh, header = readdlm(dbhf,'\t',header=true)
        dbh[:,4]/=2 # diameter to radius
        dbh[:,4]/=100 # cm to m

        dbh   = clipdat(dbh,bounds,0)
        dbh_e = findelev(dtm,dbh[:,1],dbh[:,2])
        tsm   = calculate_trunks(dbh,30,0.05,dbh_e)
        tsm[:,3] .+= tsm[:,4]
        tsm = tsm[:, 1:end .!= 4]

    end

    # load the ltc
    if branches
        ltc, _ = readdlm(ltcf,'\t',header=true)
        replace!(ltc, -9999=>NaN)
        ltc = clipdat(ltc,bounds,0)

        bsm = make_branches(ltc)
        bsm[:,3] .+= findelev(dtm,bsm[:,1],bsm[:,2])
        bsm = bsm[:, 1:end .!= 4]
    end


    ###############################################################################
    # > Calculate slope and aspect from dtm data
    if tilt
        pts_slp = findelev(hcat(dtm[:,1:2],dtm[:,4]),pts[:,1],pts[:,2])
        pts_asp = findelev(hcat(dtm[:,1:2],dtm[:,5]),pts[:,1],pts[:,2])
    else
        pts_slp = zeros(size(pts[:,1]))
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

    # dem terrain mask
    if tershad < 2 && terrain == 2
        demcrt = dem2pol(dem,pts[:,1],pts[:,2],ch,terrain_peri,dem_cellsize,0.0)

        # add dem into empty matrix
        demcrt = pol2cart(dempol[:,1],dempol[:,2])
        idx    = findpairs(kdtree,demcrt,1,kdtreedims)
        imdx   = float(reshape(idx,(radius*2,radius*2)))
        g_img[imdx.==1] .= 0

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

        @simd for crx = crxstart:size(pts,1)

            # try # catch statement to close dataset if error

                if progress; start = time(); end


                #### transfer point clouds to polar coordinates
                if branches
                    pt_dsm = getsurfdat(dsm,bsm,pts[crx,1],pts[crx,2],surf_peri);
                else
                    pt_dsm = dsm[dist(dsm,pts[crx,1],pts[crx,2]).< surf_peri,:];
                end
                pt_dsm = pcd2pol2cart(pt_dsm,pts[crx,1],pts[crx,2],pts_e[crx,1],surf_peri,"surface",ch,pts_slp[crx],0);
                pt_dsm = prepcrtdat(pt_dsm,3);

                if trunks
                    pt_tsm = pcd2pol2cart(tsm[dist(tsm,pts[crx,1],pts[crx,2]).< surf_peri*0.5,:],pts[crx,1],pts[crx,2],pts_e[crx,1],Int.(surf_peri*0.5),"surface",ch,pts_slp[crx],0);
                    tidx = findall(dist(dbh,pts[crx,1],pts[crx,2]) .< 4)
                    if size(tidx,1) > 0
                        hdt = dist(dbh[tidx,:],pts[crx,1],pts[crx,2])
                        npt  = zeros(size(tidx)) .* NaN
                        hint = zeros(size(tidx)) .* NaN
                        for tixt = 1:1:size(tidx,1)
                            if hdt[tixt] < 1
                                npt[tixt] = Int.(150); hint[tixt] = 0.005
                            else
                                npt[tixt] = Int.(100); hint[tixt] = 0.01
                            end
                        end
                        pt_tsm = vcat(pt_tsm,pcd2pol2cart(calculate_trunks(dbh[tidx,:],npt,hint,dbh_e[tidx,:]),pts[crx,1],pts[crx,2],pts_e[crx,1],Int.(surf_peri*0.5),"surface",ch,0))
                    end
                    pt_tsm = prepcrtdat(pt_tsm,3)
                end

                if tershad < 3
                    pt_dtm =  pcd2pol2cart(dtm,pts[crx,1],pts[crx,2],pts_e[crx,1],Int.(3*dem_cellsize),"terrain",ch,pts_slp[crx],dtm_cellsize);
                    if tershad < 2 && terrain == 1
                        # pt_dem = Array{Float64,2}(undef,size(dem,1),3)
                        pt_dem = dem2pol(dem,pts[crx,1],pts[crx,2],ch,terrain_peri,dem_cellsize,pts_slp[crx])
                        pt_dtm = prepterdat(vcat(pt_dtm,pt_dem));
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
                    if crx != 1; rm(outdir*"/ProgressLastPoint/"*progtext1); end
                    global progtext1 = "1. Transferring to polar took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                    writedlm(outdir*"/ProgressLastPoint/"*progtext1,NaN)
                end

                ######################
                ### Classify image
                if progress; start = time(); end

                mat2ev  = copy(g_img);

                # occupy matrix with surface points
                for zdx = 1:1:size(rbins,1)-1
                    ridx = findall(rbins[zdx] .<= pt_dsm[:,3] .< rbins[zdx+1])
                    imdx = reshape(findpairs(kdtree,pt_dsm[ridx,1:2],tol[zdx],kdtreedims,40),(radius*2,radius*2))
                    mat2ev[imdx.==1] .= 0
                end

                # add trunks
                if trunks
                    mat2ev = fillmat(kdtree,pt_tsm[:,1:2],2.0,kdtreedims,15,radius,mat2ev)
                end

                #occupy matrix with dtm
                if tershad < 3
                    mat2ev = fillmat(kdtree,pt_dtm,1.0,kdtreedims,10,radius,mat2ev);
                end

                mat2ev[isnan.(g_rad)] .= 1;

                if progress
                    elapsed = time() - start
                    if crx != 1; rm(outdir*"/ProgressLastPoint/"*progtext2); end
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
                    sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts[crx,1],pts[crx,2],loc_time,time_zone,coor_system,utm_zone)
                    swrtot, swrdir, swr_dif, transfor = calculateSWR(float.(mat2ev),loc_time,radius,sol_phi,
                                                                sol_tht,Vf_w,g_coorpol,sol_sinelev)
                end

                if progress
                    elapsed = time() - start
                    if crx != 1; rm(outdir*"/ProgressLastPoint/"*progtext3); end
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
                    if crx != 1; rm(outdir*"/ProgressLastPoint/"*progtext4); end
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

    finally
        dataset.close()
        if save_images
            images.close()
        end
        println("done with: "*taskID)
    end
end
