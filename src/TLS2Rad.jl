function TLS2Rad(pts,dat_in,par_in,exdir,taskID=empty)


    ################################################################################
    # Initialise
    eval(extract(dat_in))
    eval(extract(par_in))

    ################################################################################

    # load the TLS high res
    tsm = readlas(tsmf)
    tsm = clipdat(tsm,pts,surf_peri)

    # load the decimated las
    dsm = readlas(dsmf)
    dsm = clipdat(dsm,pts,surf_peri)


    # load the DTM?
    dtm, dtm_cellsize = importdtm(dtmf,tilt)


    ###############################################################################
    # > Preprocess DSM/DTM

    # clip dsm within eval-peri of min/max pts
    dtm = clipdat(dtm,dsm,surf_peri*4)

    # determine ground elevation of points
    pts_e    = findelev(dtm,pts[:,1],pts[:,2])























end
