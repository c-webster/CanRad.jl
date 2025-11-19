function chm2rad!(pts::Matrix{Float64},dat_in::Dict{String, String},par_in::Dict{String, Any},
    exdir::String,taskID="task",pts_z=-1)

    ################################################################################
    # Initialise
    # run compatability check then extract settings
    update_deprecated_settings!(par_in,dat_in)

    fp = FILEPATHS(; (Symbol.(keys(dat_in)) .=> values(dat_in))... )
    st = SETTINGS(; (Symbol.(keys(par_in)) .=> values(par_in))... )

    check_conflicts(st,fp,"c2r")

    st.step_progress && (start = time())

    # load the required constants
    canrad  = CANRAD()
    fileio  = FILEIO(time_zone=st.time_zone,calc_swr=st.calc_swr)
    chm2rad = CHM2RAD()

    ################################################################################
    # > Organise the output file structre and initiate process reporting

    outdir, outstr  = organise_outf(taskID,exdir,st.batch)
    global outtext = "Processing "*cfmt.("%.$(0)f", 0)*"% ... "*string(0)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    ################################################################################
    # > Organise the point vectors depending on input into the function

    # separate the points to vectors
    pts_x, pts_y =[pts[:,x] for x in 1:size(pts,2)]

    # get model number (1 = canopy, 2 = terrain only)
    size(pts,2) > 2 ? pts_m = pts[:,3] : pts_m = ones(size(pts_x))

    # fill pts_z with image_height if not given as pts_z vector
    (pts_z == -1) && (pts_z = fill(st.image_height,size(pts_x,1)))

    ################################################################################
    # > Import surface data

    limits_canopy = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,st.forest_peri)

    chm_x, chm_y, chm_z, chm_cellsize = read_griddata_window(fp.chmf,limits_canopy,true,true)

    if isempty(chm_x) && st.batch # create a dataset of zeros if chm data is unavailable for tile

        warning("$(taskID): chm empty for requested area. running with no canopy")

        limits_chmfill = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,5)

        chm_cellsize = 2.0
        chm_x = vec((limits_chmfill[1]:chm_cellsize:limits_chmfill[2]-1)'  .* ones(Int(limits_chmfill[4]-limits_chmfill[3])))
        chm_y = vec(ones(Int(limits_chmfill[2]-limits_chmfill[1]))' .* (limits_chmfill[3]:chm_cellsize:limits_chmfill[4]-1))
        chm_z = zeros(size(chm_x))

    elseif isempty(chm_x) && !st.batch

        error("CHM file is empty or missing data over the requested area. Please check the input CHM file.")

    end

    if st.special_implementation == "oshd-alps"
        # force a resolution change from 10 to 2m 
        chm_x, chm_y, chm_z, chm_cellsize = change_res(chm_x,chm_y,chm_z)
        rbins_chm = collect(1:sqrt(2).*chm_cellsize:st.forest_peri)
        chm_cellsize *= 0.5
    else
        rbins_chm = collect(4*chm_cellsize:sqrt(2).*chm_cellsize:st.forest_peri)
    end

    # if special_implementation == "aso-sites"
    #     chm_z .= ifelse.(chm_z .< 0 .|| isinf.(chm_z), 0.0, chm_z)
    # end

    ################################################################################
    # > Import prepare terrain data

    if st.terrain_highres || st.terrain_lowres || st.hlm_precalc || st.terrainmask_precalc
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
        flat_terrain_flag = false
        if st.hlm_precalc
            hlm_tht = load_hlm(ter2rad,fp.hlmf,pts_x,pts_y)
        elseif st.terrainmask_precalc
            terrain_mask = getterrainmask(canrad,fp.termf,pts_x,pts_y)
        end
    else
        flat_terrain_flag = true
    end

    if st.save_horizon
        hrdtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
        lrdtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
    end

    if st.terrain_highres && !(st.hlm_precalc || st.terrainmask_precalc)

        @unpack limits_highres, pts_e = ter2rad
        getlimits!(limits_highres,pts_x,pts_y,st.highres_peri)
        hrdtm_x, hrdtm_y, hrdtm_z, hrdtm_cellsize = read_griddata_window(fp.hrdtmf,limits_highres,true, true)
        rbins_hrdtm = collect(2*hrdtm_cellsize:sqrt(2).*hrdtm_cellsize:st.highres_peri)
        pts_e = findelev!(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),pts_x,pts_y,limits_highres,10.0,pts_e,"nearest")

    end

    # Get the low-res terrain data
    if st.terrain_lowres && !(st.hlm_precalc || st.terrainmask_precalc)

        @unpack limits_lowres, pts_e_lrdtm = ter2rad
        getlimits!(limits_lowres,pts_x,pts_y,st.lowres_peri)
        lrdtm_x, lrdtm_y, lrdtm_z, lrdtm_cellsize = read_griddata_window(fp.lrdtmf,limits_lowres,true,true)
        
        if st.terrain_highres
            rbins_lrdtm = collect(st.highres_peri-lrdtm_cellsize:sqrt(2).*lrdtm_cellsize:st.lowres_peri)
        else # include lowres terrain under pts if not using highres terrain
            rbins_lrdtm = collect(lrdtm_cellsize:sqrt(2).*lrdtm_cellsize:st.lowres_peri)
        end

        findelev!(copy(lrdtm_x),copy(lrdtm_y),copy(lrdtm_z),pts_x,pts_y,limits_lowres,lrdtm_cellsize*3,pts_e_lrdtm,"nearest")

    end

    # Get elevation of surface and evaluation points
    if st.terrain_highres
        chm_e     = findelev!(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),chm_x,chm_y,limits_highres,10,Vector{Float64}(undef,size(chm_x)))
    elseif st.terrain_lowres # use the lower resolution dataset if a higres isn't available
        chm_e     = findelev!(copy(lrdtm_x),copy(lrdtm_y),copy(lrdtm_z),chm_x,chm_y,limits_lowres,lrdtm_cellsize*3,Vector{Float64}(undef,size(chm_x)))
        pts_e     = copy(pts_e_lrdtm)
    else # if there's no terrain input, just assume the ground is flat
        pts_e     = zeros(size(pts_x))
        chm_e     = zeros(size(chm_x))
    end

    # Calculate slope and aspect from dtm data
    if st.tilt # currently disabled in initialise.jl and constants.jl
        pts_slp = findelev(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_s),pts_x,pts_y)
        # pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x,1))
    end

    # load the buildings
    if st.buildings
        bhm_x, bhm_y, bhm_z, bhm_cellsize = read_griddata_window(fp.bhmf,limits_canopy,true,true)
        rbins_bhm = collect(2*bhm_cellsize:sqrt(2).*bhm_cellsize:st.forest_peri)
        if !isempty(bhm_x)
            buildings_flag= true
    	else; buildings_flag= false; end
    else; buildings_flag= false; end

    ###############################################################################
    # > Set up the canopy structure information

    if st.special_implementation in ("swissrad","oshd")

        limits_tile = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,0.0)

        # define forest type for the tile (broadleaf or needleleaf)
        # make sure higher elevation areas are larch, not broadleaf
        if median(pts_e) > 1500
            forest_flag = 2
        else # load the forest type data for Switzerland
            # get the dominant forest type for the tile (with check for incomplete
            #   tiles along the border)
            types_tile = read_griddata_window(fp.fcdf,limits_tile,true,true)[3]
            if !isempty(types_tile[types_tile .> 0])
                forest_flag = mode(types_tile[types_tile .> 0])
            else # empty data tiles on the border get to be broadleaf
                forest_flag = 1
            end
        end

        if (st.phenology == "leafoff") || (st.phenology == "both") # LA value varies with mix ratio of deciduous/evergreen

            # load forest mix ratio and correct for lavd
            # get the range of LA values for the tile
            if forest_flag == 1 # broadleaf (broadleaf -> conifer)
                mr_x, mr_y, mr_val, _ = read_griddata_window(fp.mrdf,limits_canopy,true,true)
            elseif forest_flag == 2 # needleleaf (larch -> conifer)
            # check incompatibilities between mix rate and forest type
                # if mix rate says evergreen (<50%), but copernicus says deciduous, force larch in mr_val 75%
                # but only above 1500m
                mr_x, mr_y, mr_val, _ = read_griddata_window(fp.mrdf,limits_canopy,true,false)
                if median(pts_e) > 1500
                    _, _, ft_val,_ = read_griddata_window(fp.ftdf,limits_canopy,true,false)
                    tmp_dx = findall((mr_val .< 5000) .& (ft_val .== 1))
                    if sum(tmp_dx) > 0
                        mr_val[tmp_dx] .= 7500
                    end
                    
                    # vectorise the mr data
                    mr_x = vec(mr_x); mr_y = vec(mr_y); mr_val = vec(reverse(mr_val,dims=1))
                    rows = findall(isnan,mr_val)
                    deleteat!(mr_x,rows); deleteat!(mr_y,rows); deleteat!(mr_val,rows)

                end
            end

            temp_LA_L = reverse(collect(0.04:(0.2-0.04)/9999:0.2)) # low elevation
            temp_LA_M = reverse(collect(0.02:(0.67-0.02)/9999:0.67)) # mid elevation
            temp_LA_H = reverse(collect(0.2:(0.67-0.2)/9999:0.67)) # high elevation
            
            rows = findall(isnan,mr_val)
            deleteat!(mr_x,rows); deleteat!(mr_y,rows); deleteat!(mr_val,rows)

            mr_z = findelev(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),mr_x,mr_y,10)

            mr_val[mr_val .== 0] .= 1.0
            lavd_val = fill(0.0,size(mr_val))

            for vx in eachindex(mr_val)
                if mr_z[vx] <= 900
                    lavd_val[vx] = temp_LA_L[Int.(mr_val[vx])]
                elseif 900 < mr_z[vx] <= 1500
                    lavd_val[vx] = temp_LA_M[Int.(mr_val[vx])]
                elseif mr_z[vx] > 1500
                    lavd_val[vx] = temp_LA_H[Int.(mr_val[vx])]
                end
            end

            chm_lavd_leafoff = findelev(copy(mr_x),copy(mr_y),copy(lavd_val),chm_x,chm_y,10,"cubic")

            cbh_leafoff = 0.0

        end

        if (st.phenology == "leafon") || (st.phenology == "both") # constant LA value across all canopy pixels

            mr_x, mr_y, mr_val, _ = read_griddata_window(fp.mrdf,limits_canopy,true,true)
            mr_val[mr_val .== 0] .= 1.0
            temp_LA = reverse(collect(1.2:(0.67-1.2)/9999:0.67))
            lavd_val = fill(0.0,size(mr_val))

            for vx in eachindex(mr_val)
                lavd_val[vx] = temp_LA[Int.(mr_val[vx])]
            end

            chm_lavd_leafon = findelev(copy(mr_x),copy(mr_y),copy(lavd_val),chm_x,chm_y,10,"cubic")

            if forest_flag == 1 # broadleaf
                # chm_lavd = fill(1.2,size(chm_x))
                cbh_leafon = 0.0
            elseif forest_flag == 2 # needleleaf
                # chm_lavd = fill(0.67,size(chm_x))
                cbh_leafon = 2.0
            end

        end

    elseif st.special_implementation == "oshd-alps"

        # load the forest type data
        ft_x, ft_y, for_typ = read_griddata_window(fp.ftdf,limits_canopy,true,true)
        for_typ = change_res(ft_x,ft_y,for_typ)[3]

        # get leafoff lavd everywhere
        cd_m = 6.98 .+(0.0612 .* chm_z)
        dbh_t = (0.974 .* (chm_z .* cd_m) .^ 0.748)
        chm_lavd_leafoff = (exp.((-8.31).+(2.61.*(log.(dbh_t.*10)))+(-0.07.*chm_z))) ./  (0.8 .* pi .* ((cd_m./2).^2) .* (chm_z.*0.8))
        chm_lavd_leafoff[chm_z .< 1] .= 0

        chm_base_leafoff = fill(0.0,size(chm_z))
    
        # now replace where for_typ == 2 with a low value for leaf-off trees
        chm_lavd_leafoff[for_typ .== 1] .= 0.05
        
    else

        # without special implementation settings for seasonality,
        # canopy base height is the same for leaf on and off conditions
        cbh_leafoff = copy(st.cbh)
        cbh_leafon  = copy(st.cbh)
        # eventually could allow flexibility as an input variable if needed

    end

    chm_z .+= chm_e # normalised to absolute elevation (chm point are now dsm points)

    # fetch leaf area density dataset if supplied 
    if !isempty(fp.lavdf)

        if st.phenology == "leafon"
            chm_lavd_leafon = read_griddata_window(fp.lavdf,limits_canopy,true,true)[3]
        elseif st.phenology == "leafoff"
            chm_lavd_leafoff = read_griddata_window(fp.lavdf,limits_canopy,true,true)[3]
        elseif st.phenology == "none"
            chm_lavd_egreen = read_griddata_window(fp.lavdf,limits_canopy,true,true)[3]
        end

    end

    # create the canopy base height vector based on above settings
    if !isempty(fp.cbhf)

        if st.phenology == "leaf_on"
            chm_base_leafon = read_griddata_window(fp.cbhf,limits_canopy,true,true)[3] .+ chm_e
        elseif st.phenology == "leafoff"
            chm_base_leafoff = read_griddata_window(fp.cbhf,limits_canopy,true,true)[3] .+ chm_e
        elseif st.phenology == "none"
            chm_base_egreen = read_griddata_window(fp.cbhf,limits_canopy,true,true)[3] .+ chm_e
        end

    elseif @isdefined(chm_base_leafoff) && (st.special_implementation == "oshd-alps")

        chm_base_leafoff .+= chm_e # canopy base height already defined above

    elseif st.cbh > 0 || st.special_implementation in ("swissrad","oshd")

        if st.forest_type == "evergreen"

            chm_base_egreen = fill(st.cbh,size(chm_z))
            chm_base_egreen .+= chm_e

        elseif st.forest_type in ("deciduous", "mixed")

            if st.phenology in ("leafon", "both")
                chm_base_leafon = fill(cbh_leafon,size(chm_z))
                chm_base_leafon .+= chm_e
            end
        
            if st.phenology in ("leafoff", "both")
                chm_base_leafoff = fill(cbh_leafoff,size(chm_z))
                chm_base_leafoff .+= chm_e
            end
        end

    end

    if !isempty(fp.cbhf) || st.cbh > 0 || (st.special_implementation in ("swissrad","oshd","oshd-alps"))
        cbh_flag = true
    else
        cbh_flag = false
    end


    # load the dbh
    if st.trunks
        limits_trees =  getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,st.trunks_peri)
        dbh_x, dbh_y, dbh_z, dbh_r, _ = loaddbh(fp.dbhf,limits_trees)
        if !isempty(dbh_x)
            dbh_e = findelev!(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),dbh_x,dbh_y,limits_trees,10.0,Vector{Float64}(undef,size(dbh_x))            )
            tsm_x, tsm_y, tsm_z  = calculate_trunks(dbh_x,dbh_y,dbh_z,dbh_r,30,0.1,dbh_e)
            trunks_flag = true
        else; trunks_flag = false; end
        rbins = collect(0:(forest_peri-0)/5:forest_peri)
        knum_t = Int.(round.(collect(15:(diff([10,1])/5)[1]:1)))
    else; trunks_flag = false;
    end

    ###############################################################################
    # > Image matrix  preparation

    # create an empty image matrix
    @unpack radius, mat2ev, g_rad, g_coorcrt = canrad
    g_rad[g_rad .> radius] .= NaN
    outside_img = isnan.(g_rad)
    g_coorcrt .= ((g_coorcrt .- radius) ./ radius) .* 90

    # make g_coorcrt a KDtree for easy look up
    kdtree = scipyspat.cKDTree(g_coorcrt)

    @unpack ring_radius, ring_tht, surf_area_p, surf_area_h, relevant_pix = canrad
    for rix = 1:1:size(ring_radius,1)-1
        relevant_pix[:,rix] = (ring_radius[rix] .< vec(g_rad) .< ring_radius[rix+1])
        surf_area_p[rix] = sin(ring_tht[rix+1]/360*2*pi)^2 - sin(ring_tht[rix]/360*2*pi)^2
        surf_area_h[rix] = 2*pi*(cos(ring_tht[rix]/360*2*pi) - cos(ring_tht[rix+1]/360*2*pi))/2/pi
    end

    if st.forest_type == "evergreen"
        mat2ev_egreen = copy(mat2ev)
    elseif st.forest_type in ("deciduous", "mixed")
        if st.phenology in ("leafon", "both")
            mat2ev_leafon = copy(mat2ev)
        end
        if st.phenology in ("leafoff", "both")
            mat2ev_leafoff = copy(mat2ev)
        end
    end
    
    st.calc_terrain && (mat2ev = copy(mat2ev))

    ###############################################################################
    # > Initialise variables for transmissivity and swr calculations

    if st.calc_trans
        loc_time     = collect(Dates.DateTime(st.t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(st.t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(st.t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(st.tstep):Dates.DateTime(st.t2,"dd.mm.yyyy HH:MM:SS"))
        solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = st.tstep, radius = canrad.radius, time_zone = st.time_zone)

        @unpack trans_for = solar

        if st.calc_swr > 0
            radiation = RADIATION(loc_time = loc_time)
            if st.calc_swr == 2
                swrdat   = readdlm(fp.swrf)
                radiation.swr_open = float.(swrdat[:,3])
            end
        end
    end

    ###############################################################################
    # > Create the output file structure

    if st.calc_trans
        dataset = createfiles(outdir,outstr,pts,st,fileio,loc_time_agg)
        pts_lat, pts_lon = get_latlon(pts_x,pts_y,st.epsg_code)
    else    
        dataset = createfiles(outdir,outstr,pts,st,fileio)
    end

    st.save_images && (images = create_exmat(outdir,outstr,pts,canrad.mat2ev,st))

    st.save_horizon && (hlm = create_exhlm(outdir,outstr,pts,ter2rad))

    ###############################################################################
    # start detailed progress reporting
    if st.step_progress
        progdir = joinpath(outdir,"ProgressLastPoint")
        progtextinit = "0. Pre-calc took "*cfmt.("%.$(2)f", time() - start)*" seconds"
        ispath(progdir) ? (rm(progdir,recursive=true); mkdir(progdir)) : mkdir(progdir)
        writedlm(joinpath(progdir,progtextinit),NaN)            
    end

    ###############################################################################
    # > Loop through the points

    @simd for crx in eachindex(pts_x)

        st.step_progress && (start = time())

        if !st.terrainmask_precalc || !flat_terrain_flag

            if st.terrain_highres
                # get the high-res local terrain
                pt_hrdtm_x, pt_hrdtm_y, pt_hrdtm_z = getsurfdat(copy(hrdtm_x),copy(hrdtm_y),copy(hrdtm_z),
                                                pts_x[crx],pts_y[crx],pts_e[crx],st.highres_peri);
                pt_hrdtm_x, pt_hrdtm_y = pcd2pol2cart!(ter2rad,pt_hrdtm_x, pt_hrdtm_y, pt_hrdtm_z,pts_x[crx],
                                        pts_y[crx],pts_e[crx],"terrain",rbins_hrdtm,pts_z[crx])

                st.save_horizon && copy!(hrdtm_mintht,ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
        
            end

            if st.terrain_lowres
                # get the low-res regional terrain
                pt_lrdtm_x, pt_lrdtm_y, pt_lrdtm_z = getsurfdat(copy(lrdtm_x),copy(lrdtm_y),copy(lrdtm_z),
                                                pts_x[crx],pts_y[crx],pts_e[crx],st.lowres_peri);
                pt_lrdtm_x, pt_lrdtm_y = pcd2pol2cart!(ter2rad,pt_lrdtm_x, pt_lrdtm_y, pt_lrdtm_z,pts_x[crx],
                                        pts_y[crx],pts_e_lrdtm[crx],"terrain",rbins_lrdtm,pts_z[crx]);

                if st.save_horizon
                    if st.terrain_highres && !isempty(hrdtm_x)
                        copy!(lrdtm_mintht,ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
                        hlm["tht"][:,crx] = Int8.(round.(minimum(hcat(hrdtm_mintht,lrdtm_mintht),dims=2)))
                    else
                        hlm["tht"][:,crx] = Int8.(round.(copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])))
                    end
                end
        
            end

            if st.hlm_precalc # use pre-calculated horizon line
                pt_hrdtm_x,pt_hrdtm_y = hlm2cart(ter2rad,hlm_tht[:,crx])
            elseif st.terrain_highres && st.terrain_lowres # combine highres and lowres horizon lines
                prepterdat!(append!(pt_hrdtm_x,pt_lrdtm_x),append!(pt_hrdtm_y,pt_lrdtm_y));
            elseif st.terrain_highres # use highres only
                prepterdat!(pt_hrdtm_x,pt_hrdtm_y)
            elseif st.terrain_lowres # use lowres only
                pt_hrdtm_x, pt_hrdtm_y = prepterdat(pt_lrdtm_x,pt_lrdtm_y)
            end

        end

        if buildings_flag

            pt_bhm_x, pt_bhm_y, pt_bhm_z = getsurfdat(copy(bhm_x),copy(bhm_y),copy(bhm_z),pts_x[crx],pts_y[crx],pts_e[crx],st.highres_peri);
            
            pt_bhm_x, pt_bhm_y =  pcd2pol2cart!(ter2rad,pt_bhm_x,pt_bhm_y,pt_bhm_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                "buildings",rbins_bhm,pts_z[crx],pts_slp[crx])
            if !terrainmask_precalc
                prepterdat!(append!(pt_hrdtm_x,pt_bhm_x),append!(pt_hrdtm_y,pt_bhm_y));
            else
                pt_hrdtm_x, pt_hrdtm_y = prepterdat(pt_bhm_x,pt_bhm_y)
            end

        end

        if pts_m[crx] .== 1 # if running the forest model for this point

            if st.forest_type == "evergreen"

                pt_chm_x, pt_chm_y, pt_chm_z, pt_chm_b, pt_lavd_egreen = getsurfdat_chm(copy(chm_x),copy(chm_y),copy(chm_z),
                    copy(chm_base_egreen),copy(chm_lavd_egreen),pts_x[crx],pts_y[crx],pts_e[crx],st.forest_peri)

                if st.leaf_type in ("needleleaf","both")

                    # pts from the CHM
                    pt_chm_x_pts, pt_chm_y_pts = pcd2pol2cart!(copy(pt_chm_x),copy(pt_chm_y),
                        copy(pt_chm_z),pts_x[crx],pts_y[crx],pts_e[crx],pts_z[crx],chm_cellsize) 

                end

                pt_chm_x, pt_chm_y, pt_chm_x_thick, pt_chm_y_thick = 
                    calcCHM_Ptrans!(chm2rad,pt_chm_x,pt_chm_y,pt_chm_z,pt_chm_b,pt_lavd_egreen,
                    pts_x[crx],pts_y[crx],pts_e[crx],pts_z[crx],chm_cellsize,rbins_chm,cbh_flag) # calculated points


            elseif st.forest_type in ("deciduous", "mixed")
            
                if st.phenology in ("leafon", "both")

                    pt_chm_x, pt_chm_y, pt_chm_z, pt_chm_b_leafon, pt_lavd_leafon =
                        getsurfdat_chm(copy(chm_x),copy(chm_y),copy(chm_z),
                        copy(chm_base_leafon),copy(chm_lavd_leafon),pts_x[crx],pts_y[crx],pts_e[crx],st.forest_peri)

                    pt_chm_x_leafon, pt_chm_y_leafon, pt_chm_x_thick_leafon, pt_chm_y_thick_leafon =
                        calcCHM_Ptrans!(chm2rad,pt_chm_x,pt_chm_y,pt_chm_z,pt_chm_b_leafon,pt_lavd_leafon,
                        pts_x[crx],pts_y[crx],pts_e[crx],pts_z[crx],chm_cellsize,rbins_chm,cbh_flag) # calculated points

                end
            
                if st.phenology in ("leafoff", "both")

                    pt_chm_x, pt_chm_y, pt_chm_z, pt_chm_b_leafoff, pt_lavd_leafoff =
                        getsurfdat_chm(copy(chm_x),copy(chm_y),copy(chm_z),
                        copy(chm_base_leafoff),copy(chm_lavd_leafoff),pts_x[crx],pts_y[crx],pts_e[crx],st.forest_peri)

                    pt_chm_x_leafoff, pt_chm_y_leafoff, pt_chm_x_thick_leafoff, pt_chm_y_thick_leafoff =
                        calcCHM_Ptrans!(chm2rad,pt_chm_x,pt_chm_y,pt_chm_z,pt_chm_b_leafoff,pt_lavd_leafoff,
                        pts_x[crx],pts_y[crx],pts_e[crx],pts_z[crx],chm_cellsize,rbins_chm,cbh_flag) # calculated points
                    
                end

                if st.leaf_type in ("needleleaf","both")

                    # pts from the CHM
                    pt_chm_x_pts, pt_chm_y_pts = pcd2pol2cart!(copy(pt_chm_x),copy(pt_chm_y),copy(pt_chm_z),
                        pts_x[crx],pts_y[crx],pts_e[crx],pts_z[crx],chm_cellsize) 

                end

            end

        end

        if trunks_flag
            pt_tsm_x, pt_tsm_y, pt_tsm_z = getsurfdat(copy(tsm_x),copy(tsm_y),copy(tsm_z),pts_x[crx],pts_y[crx],
                pts_e[crx],st.forest_peri)
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
                pt_tsm_x, pt_tsm_y, pt_tsm_z = pcd2pol2cart!(append!(pt_tsm_x,tsm_tmp[1]),append!(pt_tsm_y,tsm_tmp[2]),
                    append!(pt_tsm_z,tsm_tmp[3]),pts_x[crx],pts_y[crx],pts_e[crx],"trunks",pts_z[crx])

            else
                pcd2pol2cart!(pt_tsm_x,pt_tsm_y,pt_tsm_z,pts_x[crx],pts_y[crx],pts_e[crx],"trunks",pts_z[crx])
            end
        end

        if st.step_progress
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "1. Transferring to polar took "*cfmt.("%.$(2)f", time() - start)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        ######################
        ### Classify image

        st.step_progress && (start = time())

        st.terrainmask_precalc ? copy!(mat2ev,terrain_mask[:,:,crx]) : fill!(mat2ev,1)

        # create base + terrian image matrix
        fillmat!(canrad,kdtree,hcat(pt_hrdtm_x,pt_hrdtm_y),10,mat2ev);

        # occupy matrices
        if st.special_implementation in ("swissrad","oshd")
            
            copy!(mat2ev_leafoff,mat2ev)
            copy!(mat2ev_leafon,mat2ev)

            if pts_m[crx] .== 1
                fillmat!(canrad,kdtree,hcat(pt_chm_x_leafoff,pt_chm_y_leafoff),30,mat2ev_leafoff)
                fillmat!(canrad,kdtree,hcat(pt_chm_x_leafon,pt_chm_y_leafon),30,mat2ev_leafon)
                if forest_flag == 2 
                    # thick canopy treated as opaque in leafon or in evergreen forests
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_thick_leafon,pt_chm_y_thick_leafon),15,mat2ev_leafon); # distance canopy is opaque and treated with terrain
                    # include canopy surface points for more definition at the tops of trees
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_pts,pt_chm_y_pts),20,mat2ev_leafon); # include canopy surface points
                end

                mat2ev_leafoff[outside_img] .= 1
                mat2ev_leafon[outside_img] .= 1

                st.save_images && (images["SHI_leafon"][:,:,crx] = mat2ev_leafon)
                st.save_images && (images["SHI_leafoff"][:,:,crx] = mat2ev_leafoff)

            end

        elseif (st.special_implementation == "oshd-alps")

            copy!(mat2ev_leafoff,mat2ev)

            if pts_m[crx] .== 1
                fillmat!(canrad,kdtree,hcat(pt_chm_x_leafoff,pt_chm_y_leafoff),30,mat2ev_leafoff)
                fillmat!(canrad,kdtree,hcat(pt_chm_x_thick_leafoff,pt_chm_y_thick_leafoff),15,mat2ev_leafoff) # distant canopy is opaque and treated with terrain
                mat2ev_leafoff[outside_img] .= 1
                st.save_images && (images["SHI_leafoff"][:,:,crx] = mat2ev_leafoff)
            end

        else

            if st.forest_type == "evergreen"

                copy!(mat2ev_egreen,mat2ev)
                fillmat!(canrad,kdtree,hcat(pt_chm_x,pt_chm_y),30,mat2ev_egreen)
                fillmat!(canrad,kdtree,hcat(pt_chm_x_thick,pt_chm_y_thick),15,mat2ev_egreen) # distant canopy is opaque and treated with terrain
                fillmat!(canrad,kdtree,hcat(pt_chm_x_pts,pt_chm_y_pts),20,mat2ev_egreen) # canopy surface points included for definition
                mat2ev_egreen[outside_img] .= 1
                st.save_images && (images["SHI_evergreen"][:,:,crx] = mat2ev_egreen)

            elseif st.forest_type in ("deciduous","mixed")

                if st.phenology in ("leafon","both")

                    copy!(mat2ev_leafon,mat2ev)
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_leafon,pt_chm_y_leafon),30,mat2ev_leafon)
                    if st.leaf_type == "needleleaf"
                        fillmat!(canrad,kdtree,hcat(pt_chm_x_pts,pt_chm_y_pts),20,mat2ev_leafon); # include canopy surface points for more definition at the tops of trees
                        fillmat!(canrad,kdtree,hcat(pt_chm_x_thick_leafon,pt_chm_y_thick_leafon),15,mat2ev_leafon); # thick canopy treated as opaque in leafon
                    end

                    if trunks_flag
                        for zdx = 1:1:size(rbins,1)-1
                            tridx = findall(rbins[zdx] .<= pt_tsm_z .< rbins[zdx+1])
                            fillmat!(canrad,kdtree,hcat(pt_tsm_x[tridx],pt_tsm_y[tridx]),knum_t[zdx],mat2ev_leafon)
                        end
                    end
                    
                    mat2ev_leafon[outside_img] .= 1
                    st.save_images && (images["SHI_leafon"][:,:,crx] = mat2ev_leafon)

                end

                if st.phenology in ("leafoff","both")

                    copy!(mat2ev_leafoff,mat2ev)
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_leafoff,pt_chm_y_leafoff),30,mat2ev_leafoff)
                    mat2ev_leafoff[outside_img] .= 1
                    st.save_images && (images["SHI_leafoff"][:,:,crx] = mat2ev_leafoff)

                end

            end

        end

        mat2ev[outside_img] .= 1
        (st.calc_terrain && st.save_images) && (images["SHI_terrain"][:,:,crx] = mat2ev)

        if st.step_progress
            if crx != 1; try; rm(joinpath(progdir,progtext2)); catch; end; end
            global progtext2 = "2. Classifying image took "*cfmt.("%.$(2)f", time() - start)*" seconds"
            writedlm(joinpath(progdir,progtext2),NaN)
        end

        ########################
        ### Perform calculations

        st.step_progress && (start = time())

        if st.calc_trans
            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],st.time_zone)
            @unpack trans_for = solar
        end

        if st.forest_type == "evergreen"

            svf_p_egreen, svf_h_egreen = calc_svf(canrad,mat2ev_egreen)
            dataset["svf_planar_evergreen"][crx] = Int8(round(svf_p_egreen*100))
            dataset["svf_hemi_evergreen"][crx]   = Int8(round(svf_h_egreen*100))

            if st.calc_trans
                fill!(trans_for,0)
                calc_transmissivity!(canrad,solar,trans_for,float(mat2ev_egreen),sol_phi,sol_tht)
                dataset["tvt_evergreen"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100));

                if st.calc_swr > 0
                    swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_egreen,st.calc_swr)
                    dataset["swr_total_evergreen"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                    dataset["swr_direct_evergreen"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                end
            
            end

        elseif st.forest_type in ("deciduous","mixed")
        
            if st.phenology in ("leafon","both")

                svf_p_leafon, svf_h_leafon = calc_svf(canrad,mat2ev_leafon)
                dataset["svf_planar_leafon"][crx] = Int8(round(svf_p_leafon*100))
                dataset["svf_hemi_leafon"][crx]  = Int8(round(svf_h_leafon*100))

                if st.calc_trans
                    fill!(trans_for,0)
                    calc_transmissivity!(canrad,solar,trans_for,float(mat2ev_leafon),sol_phi,sol_tht)
                    dataset["tvt_leafon"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

                    if st.calc_swr > 0
                        swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_leafon,st.calc_swr)
                        dataset["swr_total_leafon"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                        dataset["swr_direct_leafon"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                    end
                end

            end
        
            if st.phenology in ("leafoff","both")

                svf_p_leafoff, svf_h_leafoff = calc_svf(canrad,mat2ev_leafoff)
                dataset["svf_planar_leafoff"][crx] = Int8(round(svf_p_leafoff*100))
                dataset["svf_hemi_leafoff"][crx]   = Int8(round(svf_h_leafoff*100))

                if st.calc_trans
                    fill!(trans_for,0)
                    calc_transmissivity!(canrad,solar,trans_for,float(mat2ev_leafoff),sol_phi,sol_tht)
                    dataset["tvt_leafoff"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100));

                    if st.calc_swr > 0
                        swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_leafoff,st.calc_swr)
                        dataset["swr_total_leafoff"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                        dataset["swr_direct_leafoff"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                    end
                end

            end

        end
        
        if st.calc_terrain

            svf_p_terrain, svf_h_terrain = calc_svf(canrad,mat2ev)
            dataset["svf_planar_terrain"][crx] = Int8(round(svf_p_terrain*100))
            dataset["svf_hemi_terrain"][crx]   = Int8(round(svf_h_terrain*100))

            if st.calc_trans 
                fill!(trans_for,0)
                calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
                dataset["tvt_terrain"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100));

                if st.calc_swr > 0
                    swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_terrain,st.calc_swr)
                    dataset["swr_total_terrain"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                    dataset["swr_direct_terrain"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                end
            end
        
        end

        if st.step_progress
            if crx != 1; try; rm(joinpath(progdir,progtext3)); catch; end; end
            global progtext3 = "3. Calculations and export "*cfmt.("%.$(2)f", time() - start)*" seconds"
            writedlm(joinpath(progdir,progtext3),NaN)
        end

        # save the loop progress
        percentdone = Int(floor((crx / size(pts,1)) * 100))
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*cfmt.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end # end crx

    close(dataset)
    st.save_images && close(images)
    st.save_horizon && close(hlm)

    (st.save_images && st.make_pngs) && make_SHIs(outdir)

    # Save settings and data information to text file if not running batch
    if !st.batch 
        write_metadata(outdir,outstr,st,fp)
    end

end
