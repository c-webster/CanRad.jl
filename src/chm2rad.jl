function chm2rad!(pts::Matrix{Float64},dat_in::Dict{String, String},par_in::Dict{String, Any},
    exdir::String,taskID="task")

    ################################################################################
    # Initialise

    # run compatability check then extract settings
    compatability_check!(par_in)

    eval(extract(dat_in))
    eval(extract(par_in))

    progress && (start = time())

    # separate the points to vectors
    pts_x, pts_y =[pts[:,x] for x in 1:size(pts,2)]

    # get model number (1 = canopy, 2 = terrain only)
    size(pts,2) > 2 ? pts_m = pts[:,3] : pts_m = ones(size(pts_x))

    canrad  = CANRAD()
    chm2rad = CHM2RAD()

    ################################################################################
    # > Get constants, organise the output and initiate progress reporting

    outdir, outstr  = organise_outf(taskID,exdir,batch)
    global outtext = "Processing "*cfmt.("%.$(0)f", 0)*"% ... "*string(0)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    if calc_trans
        loc_time     = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,forest_type,season,calc_terrain,loc_time_agg,time_zone)
        pts_lat, pts_lon = calc_latlon(pts_x,pts_y,coor_system)
        solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = tstep, radius = canrad.radius, time_zone = time_zone)
    else    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr)
    end

    save_images && (images = create_exmat(outdir,outstr,pts,canrad.mat2ev,forest_type,season,calc_terrain))

    ################################################################################
    # > Import surface data

    limits_canopy = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,forest_peri)

    chm_x, chm_y, chm_z, chm_cellsize = read_griddata_window(chmf,limits_canopy,true,true)

    if isempty(chm_x) # create a dataset of zeros if chm data is unavailable for tile

        warning("chm empty for requested area. running with no canopy")

        limits_chmfill = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,5)

        chm_cellsize = 2.0
        chm_x = vec((limits_chmfill[1]:chm_cellsize:limits_chmfill[2]-1)'  .* ones(Int(limits_chmfill[4]-limits_chmfill[3])))
        chm_y = vec(ones(Int(limits_chmfill[2]-limits_chmfill[1]))' .* (limits_chmfill[3]:chm_cellsize:limits_chmfill[4]-1))
        chm_z = zeros(size(chm_x))

    end

    rbins_chm = collect(4*chm_cellsize:sqrt(2).*chm_cellsize:forest_peri)

    ################################################################################
    # > Import prepare terrain data

    # set terrain = true if any terrain is to be plotted in the final image matrix
    if terrainmask_precalc
        terrain_mask = getterrainmask(canrad,terf,pts_x,pts_y)
    elseif hlm_precalc
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
        if special_implementation == "oshd"
            oshd_flag = true
            pt_dem_x, pt_dem_y = load_hlm_oshd(hlmf)
        else
            oshd_flag = false
            hlm_tht = load_hlm(ter2rad,hlmf,pts_x,pts_y)
        end
    end

    if terrain_highres || terrain_lowres || hlm_precalc 
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
    end

    if save_horizon
        hlm = create_exhlm(outdir,outstr,pts,ter2rad)
        dtm_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
        dem_mintht = copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
    end

    if terrain_highres

        @unpack limits_highres, pts_e = ter2rad
        getlimits!(limits_highres,pts_x,pts_y,highres_peri)
        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits_highres,true, true)
        rbins_dtm = collect(2*dtm_cellsize:sqrt(2).*dtm_cellsize:highres_peri)

        pts_e = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y,limits_highres,10.0,pts_e)

    end

    # Get the low-res terrain data
    if terrain_lowres && !hlm_precalc

        @unpack limits_lowres, pts_e_dem = ter2rad
        getlimits!(limits_lowres,pts_x,pts_y,lowres_peri)
        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits_lowres,true,true)
        rbins_dem = collect(2*dem_cellsize:sqrt(2).*dem_cellsize:lowres_peri)
        findelev!(copy(dem_x),copy(dem_y),copy(dem_z),pts_x,pts_y,limits_lowres,dem_cellsize*3,pts_e_dem)

    end

    # Get elevation of surface and evaluation points
    if terrain_highres
        chm_e     = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),chm_x,chm_y,limits_highres,10,Vector{Float64}(undef,size(chm_x)))
    elseif terrain_lowres # use the lower resolution dataset if a higres isn't available
        chm_e     = findelev!(copy(dem_x),copy(dem_y),copy(dem_z),chm_x,chm_y,limits_highres,dem_cellsize*3,Vector{Float64}(undef,size(chm_x)))
    else # if there's no terrain input, just assume the ground is flat
        pts_e     = zeros(size(pts_x))
        chm_e     = zeros(size(chm_x))
    end

    chm_z .+= chm_e # normalised to absolute elevation (chm point are now dsm points)

    # Calculate slope and aspect from dtm data
    if tilt # currently disabled in Preparatory_Functions.jl
        pts_slp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_s),pts_x,pts_y)
        # pts_asp = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_a),pts_x,pts_y)
    else
        pts_slp = zeros(size(pts_x,1))
    end

    # load the buildings
    if buildings
        bhm_x, bhm_y, bhm_z, bhm_cellsize = read_griddata_window(bhmf,limits_canopy,true,true)
        rbins_bhm = collect(2*bhm_cellsize:sqrt(2).*bhm_cellsize:forest_peri)
        if !isempty(bhm_x)
            build = true
    	else; build = false; end
    else; build = false; end

    ###############################################################################
    # > Set up the canopy structure information

    if (special_implementation == "swissrad") || (special_implementation == "oshd")

        limits_tile = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,0.0)

        # define forest type for the tile (broadleaf or needleleaf)
        # make sure higher elevation areas are larch, not broadleaf
        if median(pts_e) > 1500
            forest_flag = 2
        else # load the forest type data for Switzerland
            # get the dominant forest type for the tile (with check for incomplete
            #   tiles along the border)
            types_tile = read_griddata_window(fcdf,limits_tile,true,true)[3]
            if !isempty(types_tile[types_tile .> 0])
                forest_flag = mode(types_tile[types_tile .> 0])
            else # empty data tiles on the border get to be broadleaf
                forest_flag = 1
            end
        end

        if (season == "winter") || (season == "both") # LA value varies with mix ratio of deciduous/evergreen

            # load forest mix ratio and correct for lavd
            # get the range of LA values for the tile
            if forest_flag == 1 # broadleaf (broadleaf -> conifer)
                mr_x, mr_y, mr_val, _ = read_griddata_window(mrdf,limits_canopy,true,true)
            elseif forest_flag == 2 # needleleaf (larch -> conifer)
            # check incompatibilities between mix rate and forest type
                # if mix rate says evergreen (<50%), but copernicus says deciduous, force larch in mr_val 75%
                # but only above 1500m
                mr_x, mr_y, mr_val, _ = read_griddata_window(mrdf,limits_canopy,true,false)
                if median(pts_e) > 1500
                    _, _, ft_val,_ = read_griddata_window(ftdf,limits_canopy,true,false)
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

            mr_z = findelev(copy(dtm_x),copy(dtm_y),copy(dtm_z),mr_x,mr_y,10)

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

            chm_lavd_winter = findelev(copy(mr_x),copy(mr_y),copy(lavd_val),chm_x,chm_y,10,"cubic")

            cbh_w = 0.0

        end

        if (season == "summer") || (season == "both") # constant LA value across all canopy pixels

            mr_x, mr_y, mr_val, _ = read_griddata_window(mrdf,limits_canopy,true,true)
            mr_val[mr_val .== 0] .= 1.0
            temp_LA = reverse(collect(1.2:(0.67-1.2)/9999:0.67))
            lavd_val = fill(0.0,size(mr_val))

            for vx in eachindex(mr_val)
                lavd_val[vx] = temp_LA[Int.(mr_val[vx])]
            end

            chm_lavd_summer = findelev(copy(mr_x),copy(mr_y),copy(lavd_val),chm_x,chm_y,10,"cubic")

            if forest_flag == 1 # broadleaf
                # chm_lavd = fill(1.2,size(chm_x))
                cbh_s = 0.0
            elseif forest_flag == 2 # needleleaf
                # chm_lavd = fill(0.67,size(chm_x))
                cbh_s = 2.0
            end

        end

    end

    if @isdefined(lavdf)

        chm_lavd = read_griddata_window(lavdf,limits_canopy,true,true)[3]

    end

    # create the canopy base height vector based on above settings
    if @isdefined(cbhf)

        println("using cbh raster")
        chm_b = read_griddata_window(cbhf,limits_canopy,true,true)[3]

    else

        if forest_type == "evergreen"

            chm_b = fill(cbh,size(chm_z))
            chm_b .+= chm_e

        elseif (forest_type == "deciduous") || (forest_type == "mixed")

            if (season == "summer") || (season == "both")
                chm_b_s = fill(cbh,size(chm_z))
                chm_b_s .+= chm_e
            end
        
            if (season == "winter") || (season == "both")
                chm_b_w = fill(cbh,size(chm_z))
                chm_b_w .+= chm_e
            end
        end

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

    @unpack trans_for = solar

    if forest_type == "evergreen"
        mat2ev_e = copy(mat2ev)
    elseif (forest_type == "deciduous") || (forest_type == "mixed")
        if (season == "summer") || (season == "both")
            mat2ev_s = copy(mat2ev)
        end
        if (season == "winter") || (season == "both")
            mat2ev_w = copy(mat2ev)
        end
    end
    
    calc_terrain && (mat2ev = copy(mat2ev))

    ###############################################################################
    # > SWR calculations

    if calc_swr > 0
        radiation = RADIATION(loc_time = loc_time)
        if calc_swr == 2
            swrdat   = readdlm(swrf)
            radiation.swr_open = float.(swrdat[:,3])
        end
    end

    ###############################################################################
    # start detailed progress reporting
    if progress
        elapsed = time() - start
        progdir = joinpath(outdir,"ProgressLastPoint")
        progtextinit = "0. Pre-calc took "*cfmt.("%.$(2)f", elapsed)*" seconds"
        if ispath(progdir)
            rm(progdir,recursive=true); mkdir(progdir)
        else
            mkdir(progdir)
        end
        writedlm(joinpath(progdir,progtextinit),NaN)            
    end

    ###############################################################################
    # > Loop through the points

    @simd for crx = 1:size(pts_x,1)

        progress && (start = time())

        if !terrainmask_precalc

            if terrain_highres
                # get the high-res local terrain
                pt_dtm_x, pt_dtm_y, pt_dtm_z = getsurfdat(copy(dtm_x),copy(dtm_y),copy(dtm_z),
                                                pts_x[crx],pts_y[crx],pts_e[crx],highres_peri);
                pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,pt_dtm_x, pt_dtm_y, pt_dtm_z,pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)

                save_horizon && copy!(dtm_mintht,ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
        
            end

            if terrain_lowres
                # get the low-res regional terrain
                pt_dem_x, pt_dem_y, pt_dem_z = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e[crx],lowres_peri);
                pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,pt_dem_x, pt_dem_y, pt_dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],"terrain",rbins_dem,image_height);

                if save_horizon
                    if !isempty(dtm_x) && terrain_highres
                        copy!(dem_mintht,ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])
                        hlm["tht"][:,crx] = Int8.(round.(minimum(hcat(dtm_mintht,dem_mintht),dims=2)))
                    else
                        hlm["tht"][:,crx] = Int8.(round.(copy(ter2rad.mintht[ter2rad.dx1:ter2rad.dx2-1])))
                    end
                end
        
            end

            if hlm_precalc && !oshd_flag
                pt_dtm_x,pt_dtm_y = hlm2cart(ter2rad,hlm_tht[:,crx])
            elseif terrain_highres && (terrain_lowres || (hlm_precalc && oshd_flag))
                prepterdat!(append!(pt_dtm_x,pt_dem_x),append!(pt_dtm_y,pt_dem_y));
            elseif terrain_highres
                prepterdat!(pt_dtm_x,pt_dtm_y)
            elseif (terrain_lowres && !terrain_highres)
                pt_dtm_x, pt_dtm_y = prepterdat(pt_dem_x,pt_dem_y)
            end

        end

        if build

            pt_bhm_x, pt_bhm_y, pt_bhm_z = getsurfdat(copy(bhm_x),copy(bhm_y),copy(bhm_z),pts_x[crx],pts_y[crx],pts_e[crx],300);
            
            pt_bhm_x, pt_bhm_y =  pcd2pol2cart!(ter2rad,pt_bhm_x,pt_bhm_y,pt_bhm_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                                "buildings",rbins_bhm,image_height,pts_slp[crx])
            if !terrainmask_precalc
                prepterdat!(append!(pt_dtm_x,pt_bhm_x),append!(pt_dtm_y,pt_bhm_y));
            else
                pt_dtm_x, pt_dtm_y = prepterdat(pt_bhm_x,pt_bhm_y)
            end

        end

        if pts_m[crx] .== 1 # if running the forest model for this point

            if forest_type == "evergreen"

                pt_chm_x, pt_chm_y, pt_chm_z, pt_chm_b, pt_lavd = getsurfdat_chm(copy(chm_x),copy(chm_y),copy(chm_z),
                    copy(chm_b),copy(chm_lavd),pts_x[crx],pts_y[crx],pts_e[crx],forest_peri)

                if (tree_species == "needleleaf") || (tree_species == "both")

                    # pts from the CHM
                    pt_chm_x_pts, pt_chm_y_pts = pcd2pol2cart!(copy(pt_chm_x),copy(pt_chm_y),copy(pt_chm_z),pts_x[crx],pts_y[crx],
                        pts_e[crx],image_height,chm_cellsize) 

                end

                pt_chm_x, pt_chm_y, pt_chm_x_thick, pt_chm_y_thick = calcCHM_Ptrans!(chm2rad,pt_chm_x,pt_chm_y,pt_chm_z,pt_chm_b,pt_lavd,
                                pts_x[crx],pts_y[crx],pts_e[crx],image_height,chm_cellsize,rbins_chm,cbh) # calculated points


            elseif (forest_type == "deciduous") || (forest_type == "mixed")
            
                if (season == "summer") || (season == "both")

                    pt_chm_x, pt_chm_y, pt_chm_z, pt_chm_b_s, pt_lavd_s = getsurfdat_chm(copy(chm_x),copy(chm_y),copy(chm_z),
                        copy(chm_b_s),copy(chm_lavd_summer),pts_x[crx],pts_y[crx],pts_e[crx],forest_peri)

                    pt_chm_x_s, pt_chm_y_s, pt_chm_x_thick_s, pt_chm_y_thick_s = calcCHM_Ptrans!(chm2rad,pt_chm_x,pt_chm_y,pt_chm_z,pt_chm_b_s,pt_lavd_s,
                        pts_x[crx],pts_y[crx],pts_e[crx],image_height,chm_cellsize,rbins_chm,cbh_s) # calculated points

                end
            
                if (season == "winter") || (season == "both")

                    pt_chm_x, pt_chm_y, pt_chm_z, pt_chm_b_w, pt_lavd_w = getsurfdat_chm(copy(chm_x),copy(chm_y),copy(chm_z),
                        copy(chm_b_w),copy(chm_lavd_winter),pts_x[crx],pts_y[crx],pts_e[crx],forest_peri)

                    pt_chm_x_w, pt_chm_y_w, _, _ = calcCHM_Ptrans!(chm2rad,pt_chm_x,pt_chm_y,pt_chm_z,pt_chm_b_w,pt_lavd_w,
                        pts_x[crx],pts_y[crx],pts_e[crx],image_height,chm_cellsize,rbins_chm,cbh_w) # calculated points
                    
                end

                if (tree_species == "needleleaf" || tree_species == "both")

                    # pts from the CHM
                    pt_chm_x_pts, pt_chm_y_pts = pcd2pol2cart!(copy(pt_chm_x),copy(pt_chm_y),copy(pt_chm_z),pts_x[crx],pts_y[crx],
                        pts_e[crx],image_height,chm_cellsize) 

                end

            end

        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "1. Transferring to polar took "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        ######################
        ### Classify image

        progress && (start = time())

        terrainmask_precalc ? copy!(mat2ev,terrain_mask[:,:,crx]) : fill!(mat2ev,1)

        # create base + terrian image matrix
        fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),10,mat2ev);

        # occupy matrices
        if ((special_implementation == "swissrad") || (special_implementation == "oshd")) && (forest_flag == 2)
            
            copy!(mat2ev_w,mat2ev)
            copy!(mat2ev_s,mat2ev)
            fillmat!(canrad,kdtree,hcat(pt_chm_x_w,pt_chm_y_w),30,mat2ev_w)
            fillmat!(canrad,kdtree,hcat(pt_chm_x_s,pt_chm_y_s),30,mat2ev_s)
            if forest_flag == 2 
                # thick canopy treated as opaque in summer or in evergreen forests
                fillmat!(canrad,kdtree,hcat(pt_chm_x_thick_s,pt_chm_y_thick_s),15,mat2ev_s); # distance canopy is opaque and treated with terrain
                # include canopy surface points for more definition at the tops of trees
                fillmat!(canrad,kdtree,hcat(pt_chm_x_pts,pt_chm_y_pts),20,mat2ev_s); # include canopy surface points
            end

            mat2ev_w[outside_img] .= 1
            mat2ev_s[outside_img] .= 1


        else

            if forest_type == "evergreen"

                copy!(mat2ev_e,mat2ev)
                fillmat!(canrad,kdtree,hcat(pt_chm_x,pt_chm_y),30,mat2ev_e)
                fillmat!(canrad,kdtree,hcat(pt_chm_x_thick,pt_chm_y_thick),15,mat2ev_e) # distant canopy is opaque and treated with terrain
                fillmat!(canrad,kdtree,hcat(pt_chm_x_pts,pt_chm_y_pts),20,mat2ev_e) # canopy surface points included for definition
                mat2ev_e[outside_img] .= 1
                save_images && (images["SHI_evergreen"][:,:,crx] = mat2ev_e)

            elseif (forest_type == "deciduous") || (forest_type == "mixed")

                if (season == "summer") || (season == "both")

                    copy!(mat2ev_s,mat2ev)
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_s,pt_chm_y_s),30,mat2ev_s)
                    if tree_species == "needleleaf"
                        fillmat!(canrad,kdtree,hcat(pt_chm_x_pts,pt_chm_y_pts),20,mat2ev_s); # include canopy surface points for more definition at the tops of trees
                        fillmat!(canrad,kdtree,hcat(pt_chm_x_thick_s,pt_chm_y_thick_s),15,mat2ev_s); # thick canopy treated as opaque in summer
                    end
                    mat2ev_s[outside_img] .= 1
                    save_images && (images["SHI_summer"][:,:,crx] = mat2ev_s)

                end

                if (season == "winter") || (season == "both")

                    copy!(mat2ev_w,mat2ev)
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_w,pt_chm_y_w),30,mat2ev_w)
                    mat2ev_w[outside_img] .= 1
                    save_images && (images["SHI_winter"][:,:,crx] = mat2ev_w)

                end

            end

        end

        mat2ev[outside_img] .= 1
        (calc_terrain && save_images) && (images["SHI_terrain"][:,:,crx] = mat2ev)

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext2)); catch; end; end
            global progtext2 = "2. Classifying image took "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext2),NaN)
        end

        ########################
        ### Perform calculations

        progress && (start = time())

        if calc_trans
            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],time_zone)
            @unpack trans_for = solar
        end

        if forest_type == "evergreen"

            svf_p_e, svf_h_e = calc_svf(canrad,mat2ev_e)
            dataset["svf_planar_e"][crx] = Int8(round(svf_p_e*100))
            dataset["svf_hemi_e"][crx]   = Int8(round(svf_h_e*100))

            if calc_trans
                fill!(trans_for,0)
                calc_transmissivity!(canrad,solar,trans_for,float(mat2ev_e),sol_phi,sol_tht)
                dataset["for_trans_e"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100));

                if calc_swr > 0
                    swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_e,calc_swr)
                    dataset["swr_total_e"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                    dataset["swr_direct_e"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                end
            
            end

        elseif (forest_type == "deciduous") || (forest_type == "mixed")
        
            if (season == "summer") || (season == "both")

                svf_p_s, svf_h_s = calc_svf(canrad,mat2ev_s)
                dataset["svf_planar_s"][crx] = Int8(round(svf_p_s*100))
                dataset["svf_hemi_s"][crx]  = Int8(round(svf_h_s*100))

                if calc_trans
                    fill!(trans_for,0)
                    calc_transmissivity!(canrad,solar,trans_for,float(mat2ev_s),sol_phi,sol_tht)
                    dataset["for_trans_s"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100))

                    if calc_swr > 0
                        swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_s,calc_swr)
                        dataset["swr_total_s"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                        dataset["swr_direct_s"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                    end
                end

            end
        
            if (season == "winter") || (season == "both")

                svf_p_w, svf_h_w = calc_svf(canrad,mat2ev_w)
                dataset["svf_planar_w"][crx] = Int8(round(svf_p_w*100))
                dataset["svf_hemi_w"][crx]   = Int8(round(svf_h_w*100))

                if calc_trans
                    fill!(trans_for,0)
                    calc_transmissivity!(canrad,solar,trans_for,float(mat2ev_w),sol_phi,sol_tht)
                    dataset["for_trans_w"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100));

                    if calc_swr > 0
                        swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_w,calc_swr)
                        dataset["swr_total_w"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                        dataset["swr_direct_w"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                    end
                end

            end

        end
        
        if calc_terrain

            svf_p_t, svf_h_t = calc_svf(canrad,mat2ev)
            dataset["svf_planar_t"][crx] = Int8(round(svf_p_t*100))
            dataset["svf_hemi_t"][crx]   = Int8(round(svf_h_t*100))

            if calc_trans 
                fill!(trans_for,0)
                calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
                dataset["trans_t"][:,crx] = Int8.(round.((vec(aggregate_data(solar,trans_for)))*100));

                if calc_swr > 0
                    swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,svf_p_t,calc_swr)
                    dataset["swr_total_t"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                    dataset["swr_direct_t"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
                end
            end
        
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext3)); catch; end; end
            global progtext3 = "3. Calculations and export "*cfmt.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext3),NaN)
        end

        # save the progress
        percentdone = Int(floor((crx / size(pts,1)) * 100))
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*cfmt.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end # end crx

    close(dataset)
    save_images && close(images)
    save_horizon && close(hlm)

    save_images && make_SHIs(outdir,forest_type,season,calc_terrain)

    println("done with "*taskID)

end