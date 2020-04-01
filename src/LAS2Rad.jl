function LAS2Rad(pts,dat_in,par_in,exdir,taskID=empty)


    ################################################################################
    # Initialise
    eval(extract(dat_in))
    eval(extract(par_in))

    if progress; start = time(); end

    global outdir = exdir*"/"*string(Int(pts[1,1]))*"_"*string(Int(pts[1,2]))
    if !ispath(outdir)
        mkdir(outdir)
    end
    # writedlm(outdir*"/"*basename(taskID),taskID)

    ################################################################################

    # load the las
    dsm = readlas(dsmf)

    # load the dtm
    dtm, dtm_cellsize = importdtm(dtmf,tilt)

    # load the dem
    dem, dem_cellsize = read_ascii(demf)

    ###############################################################################
    # > Preprocess DSM/DTM

    # clip dsm within eval-peri of min/max pts
    dsm = clipdat(dsm,pts,surf_peri)
    dtm = clipdat(dtm,dsm,surf_peri*4)
    dem = clipdat(dem,dsm,terrain_peri)

    # determine ground elevation of points
    pts_e    = findelev(dtm,pts[:,1],pts[:,2])
    dsm[:,4] = findelev(dtm,dsm[:,1],dsm[:,2])
    dsm[:,3] .-=  dsm[:,4]

    ###############################################################################
    # > Generate extra canopy elements

    # load the dbh
    if trunks
        dbh, header = readdlm(dbhf,'\t',header=true)
        dbh[:,4]/=2 # diameter to radius
        dbh[:,4]/=100 # cm to m

        dbh   = clipdat(dbh,dsm,0)
        dbh_e = findelev(dtm,dbh[:,1],dbh[:,2])
        tsm   = calculate_trunks(dbh,30,0.05,dbh_e)
    end

    # load the ltc
    if branches
        ltc, _ = readdlm(ltcf,'\t',header=true)
        replace!(ltc, -9999=>NaN)
        ltc = clipdat(ltc,dsm,0)

        bsm = make_branches(ltc)
        bsm[:,4] = findelev(dtm,bsm[:,1],bsm[:,2])
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
        dempol = dem2pol(dem,pts[:,1],pts[:,2],ch,terrain_peri,dem_cellsize,0)

        # add dem into empty matrix
        demcrt = pol2cart(dempol[:,1],dempol[:,2])
        idx    = findpairs(kdtree,demcrt,1,kdtreedims)
        imdx   = float(reshape(idx,(radius*2,radius*2)))
        g_img[imdx.==1] .= 0
    end

    # calculate solar track
    loc_time = collect(Dates.DateTime(t1,"yyyy.mm.dd HH:MM:SS"):Dates.Minute(int):
                                Dates.DateTime(t2,"yyyy.mm.dd HH:MM:SS"))

    if calc_swr
        sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts,loc_time,time_zone,coor_system,utm_zone)
        swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset = createfiles(outdir,pts,loc_time,t1,t2,int,calc_swr)
    else
         Vf_weighted, Vf_flat, dataset = createfiles(outdir,pts,loc_time,t1,t2,int,calc_swr)
    end

    if save_images
        SHIs, images = create_exmat(outdir,pts,g_img)
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

        @simd for crx = 1:size(pts,1)

            try # catch statement to close dataset if error

                if progress; start = time(); end

                #### transfer point clouds to polar coordinates
                dsmpol = pcd2pol(dsm,pts[crx,1],pts[crx,2],pts_e[crx,1],surf_peri,ch)
                dsmcrt = pol2cart(dsmpol[:,1], dsmpol[:,2])

                if trunks
                    tidx = findall((sqrt.(((dbh[:,1].-pts[crx,1]).^2) .+ ((dbh[:,2].-pts[crx,2]).^2))) .< 4)
                    tsmpol = pcd2pol(tsm,pts[crx,1],pts[crx,2],pts_e[crx,1],surf_peri*0.5,ch)
                    if size(tidx,1) > 0
                        hdt = sqrt.(((dbh[tidx,1].-pts[crx,1]).^2) .+ ((dbh[tidx,2].-pts[crx,2]).^2))
                        npt  = zeros(size(tidx)) .* NaN
                        hint = zeros(size(tidx)) .* NaN
                        for tixt = 1:1:size(tidx,1)
                            if hdt[tixt] < 1
                                npt[tixt] = 150; hint[tixt] = 0.005
                            else
                                npt[tixt] = 100; hint[tixt] = 0.01
                            end
                        end

                        tsmadd = calculate_trunks(dbh[tidx,:],npt,hint,dbh_e[tidx,:])

                        tsmaddpol = pcd2pol(tsmadd,pts[crx,1],pts[crx,2],pts_e[crx,1],surf_peri*0.5,ch)

                        tsmpol = vcat(tsmpol,tsmaddpol)
                    end
                    tsmcrt = pol2cart(tsmpol[:,1], tsmpol[:,2])
                    tsmcrt, tsmrad = prepcrtdat(tsmcrt,tsmpol[:,3],3)
                end

                if branches
                    bsmpol = pcd2pol(bsm,pts[crx,1],pts[crx,2],pts_e[crx,1],surf_peri*0.5,ch)
                    bsmcrt = pol2cart(bsmpol[:,1], bsmpol[:,2])
                end

                if tershad < 3
                    dtmpol = pcd2pol(dtm,pts[crx,1],pts[crx,2],pts_e[crx,1],4*surf_peri,ch)
                    dtmpol = calc_horizon_lines(1,3*surf_peri,dtmpol,pts_slp[crx])
                    tercrt = pol2cart(dtmpol[:,1],dtmpol[:,2])

                    if tershad < 2 && terrain == 1
                        dempol = dem2pol(dem,pts[crx,1],pts[crx,2],ch,terrain_peri,dem_cellsize,pts_slp[crx])
                        demcrt = pol2cart(dempol[:,1],dempol[:,2])
                        tercrt = prepterdat(tercrt,demcrt)
                    end
                end

                #### Prepare surface data
                if trunks && branches
                    datcrt = vcat(dsmcrt,bsmcrt,tsmcrt)
                    datrad = vcat(dsmpol[:,3],bsmpol[:,3],tsmrad)
                elseif trunks && ~branches
                    datcrt = vcat(dsmcrt,tsmcrt)
                    datrad = vcat(dsmpol[:,3],tsmrad)
                elseif ~trunks && branches
                    datcrt = vcat(dsmcrt,bsmcrt)
                    datrad = vcat(dsmpol[:,3],bsmpol[:,3])
                else
                    datpol = copy(dsmpol)
                end

                if tilt
                    slp = pts_slp[crx]
                    imcX, imcY = getimagecentre(slp,pts_asp[crx])
                    gcrcrt = zeros(size(g_coorcrt)) .* NaN
                    gcrcrt[:,1] = g_coorcrt[:,1] .+ imcX
                    gcrcrt[:,2] = g_coorcrt[:,2] .+ imcY
                    kdtree = scipyspat.cKDTree(gcrcrt)
                    kdtreedims = size(gcrcrt,1)
                else
                    datpol[datpol[:,2] .> 90,2] .= 90
                    slope = pts_slp[crx]
                end

                # datcrt, datrad = prepcrtdat(pol2cart(datpol[:,1],datpol[:,2]), datpol[:,3])

                if progress
                    elapsed = time() - start
                    if crx != 1; rm(outdir*"/ProgressLastPoint/"*progtext1); end
                    global progtext1 = "1. Transferring to polar took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
                    writedlm(outdir*"/ProgressLastPoint/"*progtext1,NaN)
                end

                ######################
                ### Classify image
                if progress; start = time(); end

                mat2ev  = copy(g_img)

                rbins = collect(0:(surf_peri-0)/5:surf_peri)
                tol   = tolerance.*collect(reverse(0.75:0.05:1))

                # occupy matrix with surface points
                for zdx = 1:1:size(rbins,1)-1
                    ridx = findall(rbins[zdx] .<= datrad[:,1] .< rbins[zdx+1])
                    ringcrt, _ = prepcrtdat(datcrt[ridx,:], datrad[ridx,:],2)
                    idx  = findpairs(kdtree,ringcrt,tol[zdx],kdtreedims,40)
                    imdx = float(reshape(idx,(radius*2,radius*2)))
                    mat2ev[imdx.==1] .= 0
                end

                #occupy matrix with dtm
                if tershad < 3
                    mat2ev = fillmat(kdtree,tercrt,1,kdtreedims,radius,mat2ev)
                end

                mat2ev[isnan.(g_rad)] .= 1

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
                Vf_w, Vf_f = calculateVf(mat2ev,g_rad)

                ##### Calculate SWR/forest transmissivity
                if ~tilt & calc_swr
                    swrtot, swrdir, swr_dif, transfor = calculateSWR(mat2ev,loc_time,radius,sol_phi,
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
                    SHIs[crx] = np.array(mat2ev)
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

            catch
                writedlm(outdir*"/Error_"*string(Int(pts[crx,1]))*"_"*string(Int(pts[crx,2]))*".txt",NaN)
            end

        end #end crx

    finally
        dataset.close()
        if save_images
            images.close()
        end
        println("done with: "*taskID)
    end
end
