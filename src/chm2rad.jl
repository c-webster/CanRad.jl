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

    outdir, outstr, crxstart, append_file, percentdone = organise_outf(taskID,exdir,batch,size(pts_x,1))
    crxstart = 1; append_file = false; percentdone = 0;     # force restart the tile
    global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crxstart-1)*" of "*string(size(pts,1))*".txt"
    writedlm(joinpath(outdir,outtext),NaN)

    if calc_trans
        loc_time     = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(2):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
        loc_time_agg = collect(Dates.DateTime(t1,"dd.mm.yyyy HH:MM:SS"):Dates.Minute(tstep):Dates.DateTime(t2,"dd.mm.yyyy HH:MM:SS"))
    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file,loc_time_agg,time_zone)
        pts_lat, pts_lon = calc_latlon(pts_x,pts_y,coor_system)
        solar = SOLAR(loc_time = loc_time, loc_time_agg = loc_time_agg, tstep = tstep, radius = canrad.radius, time_zone = time_zone)
    else    
        dataset = createfiles(outdir,outstr,pts,calc_trans,calc_swr,append_file)
    end

    if save_images 
        images = create_exmat(outdir,outstr,pts,canrad.mat2ev,append_file)
    end

    ################################################################################
    # > Import surface data

    limits_canopy = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,surf_peri)

    chm_x, chm_y, chm_z, chm_cellsize = read_griddata_window(chmf,limits_canopy,true,true)

    if isempty(chm_x) # create a dataset of zeros if chm data is unavailable for tile

        warning("chm empty for requested area. running with no canopy")

        limits_chmfill = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,5)

        chm_cellsize = 2.0
        chm_x = vec((limits_chmfill[1]:chm_cellsize:limits_chmfill[2]-1)'  .* ones(Int(limits_chmfill[4]-limits_chmfill[3])))
        chm_y = vec(ones(Int(limits_chmfill[2]-limits_chmfill[1]))' .* (limits_chmfill[3]:chm_cellsize:limits_chmfill[4]-1))
        chm_z = zeros(size(chm_x))

    end

    rbins_chm = collect(4*chm_cellsize:sqrt(2).*chm_cellsize:surf_peri)

    ################################################################################
    # > Import prepare terrain data

    # set terrain = true if any terrain is to be plotted in the final image matrix
    if terrainmask_precalc
        terrain_mask = getterrainmask(canrad,terf,pts_x,pts_y)
    elseif horizon_line
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
        if oshd_flag
            pt_dem_x, pt_dem_y = load_hlm_oshd(hlmf)
        elseif !oshd_flag
            hlm_tht = load_hlm(ter2rad,hlmf,pts_x,pts_y)
        end
    end

    if terrain_highres || terrain_lowres || horizon_line
        ter2rad = TER2RAD(pts_sz = size(pts_x,1))
    end

    if terrain_highres

        @unpack limits_highres, pts_e = ter2rad
        getlimits!(limits_highres,pts_x,pts_y,dtm_peri)
        dtm_x, dtm_y, dtm_z, dtm_cellsize = read_griddata_window(dtmf,limits_highres,true, true)
        rbins_dtm = collect(2*dtm_cellsize:sqrt(2).*dtm_cellsize:dtm_peri)

        pts_e = findelev!(copy(dtm_x),copy(dtm_y),copy(dtm_z),pts_x,pts_y,limits_highres,10.0,pts_e)

    end

    # Get the low-res terrain data
    if terrain_lowres && !horizon_line

        @unpack limits_lowres, pts_e_dem = ter2rad
        getlimits!(limits_lowres,pts_x,pts_y,terrain_peri)
        dem_x, dem_y, dem_z, dem_cellsize = read_griddata_window(demf,limits_lowres,true,true)
        rbins_dem = collect(2*dem_cellsize:sqrt(2).*dem_cellsize:terrain_peri)
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

    # get the tile horizon line
    if terrain_tile && (!horizon_line || !terrainmask_precalc)
        pt_dem_x, pt_dem_y, pt_dem_z = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),mean(pts_x),mean(pts_y),mean(pts_e_dem),terrain_peri);

        pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,pt_dem_x, pt_dem_y, pt_dem_z,mean(pts_x),mean(pts_y),mean(pts_e_dem),
            "terrain",rbins_dem,image_height,0.0)
    end

    # load the buildings
    if buildings
        bhm_x, bhm_y, bhm_z, bhm_cellsize = read_griddata_window(bhmf,limits_canopy,true,true)
        rbins_bhm = collect(2*bhm_cellsize:sqrt(2).*bhm_cellsize:surf_peri)
        if !isempty(bhm_x)
            build = true
    	else; build = false; end
    else; build = false; end

    ###############################################################################
    # > Set up the canopy structure information

    if @isdefined(mrdf)

        limits_tile = getlimits!(Vector{Float64}(undef,4),pts_x,pts_y,0.0)

        # define forest type for the tile (broadleaf or needleleaf)
        # make sure higher elevation areas are larch, not broadleaf
        if median(pts_e) > 1500
            for_type = 2
        else # load the forest type data for Switzerland
            # get the dominant forest type for the tile (with check for incomplete
            #   tiles along the border)
            types_tile = read_griddata_window(fcdf,limits_tile,true,true)[3]
            if !isempty(types_tile[types_tile .> 0])
                for_type = mode(types_tile[types_tile .> 0])
            else # empty data tiles on the border get to be broadleaf
                for_type = 1
            end
        end

        if season == "winter" # LA value varies with mix ratio of deciduous/evergreen

            # load forest mix ratio and correct for lavd
            # get the range of LA values for the tile
            if for_type == 1 # broadleaf (broadleaf -> conifer)
                mr_x, mr_y, mr_val, _ = read_griddata_window(mrdf,limits_canopy,true,true)
            elseif for_type == 2 # needleleaf (larch -> conifer)
            # check incompatibilities between mix rate and forest type
                # if mix rate says evergreen (<50%), but copernicus says deciduous, force larch in mr_val 75%
                # but only above 1500m
                mr_x, mr_y, mr_val, _ = read_griddata_window(mrdf,limits_canopy,true,true)
                if median(pts_e) > 1500
                    _, _, ft_val,_ = read_griddata_window(ftdf,limits_canopy,true,true)
                    tmp_dx = findall(mr_val .< 5000 .& ft_val .== 1)
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
            temp_LA_M = reverse(collect(0.02:(0.67-0.05)/9999:0.67)) # mid elevation
            temp_LA_H = reverse(collect(0.2:(0.67-0.05)/9999:0.67)) # high elevation

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

            chm_lavd = findelev(copy(mr_x),copy(mr_y),copy(lavd_val),chm_x,chm_y,10,"cubic")

            cbh = 0.0

        elseif season == "summer" # constant LA value across all canopy pixels

            mr_x, mr_y, mr_val, _ = read_griddata_window(mrdf,limits_canopy,true,true)
            mr_val[mr_val .== 0] .= 1.0
            temp_LA = reverse(collect(1.2:(0.67-1.2)/9999:0.67))
            lavd_val = fill(0.0,size(mr_val))

            for vx in eachindex(mr_val)
                lavd_val[vx] = temp_LA[Int.(mr_val[vx])]
            end

            chm_lavd = findelev(copy(mr_x),copy(mr_y),copy(lavd_val),chm_x,chm_y,10,"cubic")

            if for_type == 1 # broadleaf
                # chm_lavd = fill(1.2,size(chm_x))
                cbh = 0.0
            elseif for_type == 2 # needleleaf
                # chm_lavd = fill(0.67,size(chm_x))
                cbh = 2.0
            end

        end

    elseif @isdefined(lavdf)

        ch_x, ch_y, chm_lavd, _ = read_griddata(lavdf,true,true)
        clipdat!(copy(ch_x),copy(ch_y),chm_lavd,limits_canopy,Vector{Bool}(undef,size(chm_lavd,1)))
        cbh = 2.0

        for_type = 2 # currently assumes evergreen needleleaf but needs improvement of implementation

    end

    # create the canopy base height vector based on above settings
    chm_b = fill(cbh,size(chm_z))
    chm_b .+= chm_e

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
        progtextinit = "0. Pre-calc took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
        if ispath(progdir)
            rm(progdir,recursive=true); mkdir(progdir)
        else
            mkdir(progdir)
        end
        writedlm(joinpath(progdir,progtextinit),NaN)            
    end

    ###############################################################################
    # > Loop through the points

    @simd for crx = crxstart:size(pts_x,1)

        if progress; start = time(); end

        if !terrainmask_precalc

            if terrain_highres
                # get the high-res local terrain
                pt_dtm_x, pt_dtm_y, pt_dtm_z = getsurfdat(copy(dtm_x),copy(dtm_y),copy(dtm_z),
                                                pts_x[crx],pts_y[crx],pts_e[crx],dtm_peri);
                pt_dtm_x, pt_dtm_y = pcd2pol2cart!(ter2rad,pt_dtm_x, pt_dtm_y, pt_dtm_z,pts_x[crx],pts_y[crx],pts_e[crx],"terrain",rbins_dtm,image_height)
        
            end

            if terrain_lowres && !terrain_tile
                # get the low-res regional terrain
                pt_dem_x, pt_dem_y, pt_dem_z = getsurfdat(copy(dem_x),copy(dem_y),copy(dem_z),pts_x[crx],pts_y[crx],pts_e[crx],terrain_peri);
                pt_dem_x, pt_dem_y = pcd2pol2cart!(ter2rad,pt_dem_x, pt_dem_y, pt_dem_z,pts_x[crx],pts_y[crx],pts_e_dem[crx],"terrain",rbins_dem,image_height);
        
            end

            if horizon_line && !oshd_flag
                pt_dtm_x,pt_dtm_y = hlm2cart(ter2rad,hlm_tht[:,crx])
            elseif terrain_highres && (terrain_lowres || (horizon_line && oshd_flag))
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
           
            pt_chm_x, pt_chm_y, pt_chm_z, pt_chm_b, pt_lavd = getsurfdat_chm(copy(chm_x),copy(chm_y),copy(chm_z),
                                copy(chm_b),copy(chm_lavd),pts_x[crx],pts_y[crx],pts_e[crx],surf_peri)
            
            if pt_corr

                pt_chm_x_pts, pt_chm_y_pts = pcd2pol2cart!(copy(pt_chm_x),copy(pt_chm_y),copy(pt_chm_z),pts_x[crx],pts_y[crx],
                                pts_e[crx],image_height,chm_cellsize) # pts from the CHM

                pt_chm_x, pt_chm_y, pt_chm_x_thick, pt_chm_y_thick = calcCHM_Ptrans!(chm2rad,pt_chm_x,pt_chm_y,pt_chm_z,pt_chm_b,pt_lavd,
                                pts_x[crx],pts_y[crx],pts_e[crx],image_height,chm_cellsize,rbins_chm,cbh) # calculated points

            else
                #  100% opaque canopy:
                pt_chm_x, pt_chm_y = pcd2pol2cart!(ter2rad,pt_chm_x, pt_chm_y, pt_chm_z,pts_x[crx],pts_y[crx],pts_e[crx],
                                            "terrain",rbins_chm,image_height,pts_slp[crx])
                if !terrainmask_precalc # merge if using terrain
                    prepterdat!(append!(pt_dtm_x,pt_chm_x),append!(pt_dtm_y,pt_chm_y));
                else
                    pt_dtm_x, pt_dtm_y = prepterdat(pt_chm_x,pt_chm_y)
                end
            end
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext1)); catch; end; end
            global progtext1 = "1. Transferring to polar took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext1),NaN)
        end

        ######################
        ### Classify image
        if progress; start = time(); end

        if terrainmask_precalc
            mat2ev .= terrain_mask[:,:,crx]
        else
            fill!(mat2ev,1)
        end

        # occupy matrix
        if pt_corr
            if pts_m[crx] .== 1
                fillmat!(canrad,kdtree,hcat(pt_chm_x,pt_chm_y),30,mat2ev)
                if season == "summer" || for_type == 2 
                    # thick canopy treated as opaque in summer or in evergreen forests
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_thick,pt_chm_y_thick),15,mat2ev); # distance canopy is opaque and treated with terrain
                    # include canopy surface points for more definition at the tops of trees
                    fillmat!(canrad,kdtree,hcat(pt_chm_x_pts,pt_chm_y_pts),20,mat2ev); # include canopy surface points
                end
            end
            if !terrainmask_precalc
                fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),10,mat2ev);
            end
        else # treat all canopy as opaque (like terrain) -> all points come in as terrain points from above section
            fillmat!(canrad,kdtree,hcat(pt_dtm_x,pt_dtm_y),10,mat2ev); # use this line if plotting opaque canpoy
        end      

        mat2ev[outside_img] .= 1

        if save_images
            images["SHI"][:,:,crx] = mat2ev;
        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext2)); catch; end; end
            global progtext2 = "2. Classifying image took "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext2),NaN)
        end

        ########################
        ### Perform calculations

        if progress; start = time(); end

        ##### Calculate Vf and save
        Vf_p, Vf_h = calculateVf(canrad,mat2ev)

        dataset["Vf_planar"][crx] = Int(round(Vf_p*100))
        dataset["Vf_hemi"][crx]   = Int(round(Vf_h*100))

        ##### CalcSWR
        if calc_trans
            sol_tht, sol_phi, sol_sinelev  = calc_solar_track(solar,pts_lat[crx],pts_lon[crx],time_zone)

            @unpack trans_for = solar
            calc_transmissivity!(canrad,solar,trans_for,float(mat2ev),sol_phi,sol_tht)
            
            dataset["Forest_Transmissivity"][:,crx] = Int.(round.((vec(aggregate_data(solar,trans_for))).*100));

            if calc_swr > 0
                swrtot, swrdir = calculateSWR(radiation,trans_for,sol_sinelev,Vf_p,calc_swr)
                dataset["SWR_total"][:,crx]  = Int.(round.(vec(aggregate_data(solar,swrtot))))
                dataset["SWR_direct"][:,crx] = Int.(round.(vec(aggregate_data(solar,swrdir))))
            end

        end

        if progress
            elapsed = time() - start
            if crx != 1; try; rm(joinpath(progdir,progtext3)); catch; end; end
            global progtext3 = "3. Calculations and export "*sprintf1.("%.$(2)f", elapsed)*" seconds"
            writedlm(joinpath(progdir,progtext3),NaN)
        end

        # save the progress
        percentdone = Int(floor((crx / size(pts,1)) * 100))
        rm(joinpath(outdir,outtext))
        global outtext = "Processing "*sprintf1.("%.$(0)f", percentdone)*"% ... "*string(crx)*" of "*string(size(pts,1))*".txt"
        writedlm(joinpath(outdir,outtext),NaN)

    end # end crx

    close(dataset)
    if save_images; close(images); end

    println("done with "*taskID)

end