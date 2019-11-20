function CHM2Rad(pts,dat_in,par_in,exdir,taskID)

    ################################################################################
    # Initialise
    extract(dat_in)
    extract(par_in)

    if progress; start = time(); end

    global outdir = exdir*"/"*string(Int(pts[1,1]))*"_"*string(Int(pts[1,2]))
    if !ispath(outdir)
        mkdir(outdir)
    end
    # writedlm(outdir*"/"*basename(taskID),taskID)

    # Import terrain/surface data and clip
    chm, chm_cellsize = read_ascii(chmf*segn*"_chm.asc")
    chm = clipdat(chm,pts,chm_peri)

    dtm, dtm_cellsize = importdtm(dtmf*"DTM_"*segn*".mat")
    dtm = clipdat(dtm,pts,dtm_peri)

    if !isempty(dem)
         dem, dem_cellsize = read_ascii(demf)
         dem = clipdat(dem,pts,terrain_peri)
    end

    # calculate tht
    if phi_int == 1
        div = 180
    elseif phi_int == 0.5
        div = 360
    elseif phi_int == 0.25
        div = 720
    else
        error("phi_int variable not 1, 0.5 or 0.25")
    end
    phi = collect(-pi:pi/div:pi)

    # create an empty matrix for Vf calculation
    g_rad, g_coorpol, g_coorcrt, g_img = create_mat(radius)
    g_coorcrt = ((g_coorcrt .- radius) ./ radius) .* 90
    g_img[isnan.(g_rad)] .= 1

    # make g_coorcrt a KDtree for easy look up
    kdtree = scipyspat.cKDTree(g_coorcrt)
    kdtreedims = size(g_coorcrt,1)

    # organise file structure for export
    if !ispath(exdir)
        mkdir(exdir)
    end
    if !ispath(exdir*segn)
        mkdir(exdir*segn)
    end

    # calculate solar track
    loc_time = collect(Dates.DateTime(t1,"yyyy.mm.dd HH:MM:SS"):Dates.Minute(int):
                                Dates.DateTime(t2,"yyyy.mm.dd HH:MM:SS"))

    sol_tht, sol_phi, sol_sinelev  = calc_solar_track(pts,loc_time,time_zone,coor_system)

    # create the output files

    swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset = createfiles(outdir,pts,loc_time,t1,t2,int)


    ###############################################################################
    # > Loop through the points
    # try

        @simd for crx = 1:size(pts,1)





        end










end
