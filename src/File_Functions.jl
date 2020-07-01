extension(url::String) = match(r"\.[A-Za-z0-9]+$", url).match

function readlas(infile::String,ground=nothing)
    if extension(infile) == ".laz"
        header, dsmdat = LazIO.load(infile)
    elseif extension(infile) == ".las"
        header, dsmdat = FileIO.load(infile)
    else
        error("Unknown DSM file extension")
    end

    dsm_x = fill(NaN,size(dsmdat))
    dsm_y = fill(NaN,size(dsmdat))
    dsm_z = fill(NaN,size(dsmdat))
    dsm_c = fill(NaN,size(dsmdat))
    for dx in eachindex(dsmdat)
        dsm_c[dx] = trunc(Int,float(classification(dsmdat[dx])))
        if dsm_c[dx] == 2
            continue
        else
            dsm_x[dx] = xcoord(dsmdat[dx],header)
            dsm_y[dx] = ycoord(dsmdat[dx],header)
            dsm_z[dx] = zcoord(dsmdat[dx],header)
        end
    end
    rows = findall(isnan,dsm_x)
    deleteat!(dsm_x,rows)
    deleteat!(dsm_y,rows)
    deleteat!(dsm_z,rows)
    return dsm_x, dsm_y, dsm_z
end

function importdtm(dtmf::String,tilt::Bool)
    if tilt
        file = matopen(dtmf); dtmdat = read(file,"dtm"); close(file)
        dtm_x = dtmdat["x"]
        dtm_y = dtmdat["y"]
        dtm_z = dtmdat["z"]
        dtm_s = dtmdat["s"]
        dtm_a = dtmdat["a"]
        dtm_cellsize = dtmdat["cellsize"]
        return dtm_x, dtm_y, dtm_z, dtm_s, dtm_a, dtm_cellsize
    else
        if extension(dtmf) == ".mat"
            file = matopen(dtmf); dtmdat = read(file,"dtm"); close(file)
            dtm_x = vec(dtmdat["x"])
            dtm_y = vec(dtmdat["y"])
            dtm_z = vec(dtmdat["z"])
            dtm_cellsize = dtmdat["cellsize"]
            rows = findall(isnan,dtm_z)
            deleteat!(dtm_x,rows)
            deleteat!(dtm_y,rows)
            deleteat!(dtm_z,rows)
        elseif extension(dtmf) == ".asc" || extension(dtmf) == ".txt"
            dtm_x, dtm_y, dtm_z, dtm_cellsize = read_ascii(dtmf)
        end
    return dtm_x, dtm_y, dtm_z, dtm_cellsize
    end
end

function read_ascii(demf::String)

    f = open(demf)
    ncols     = parse(Int64,split(readline(f))[2])
    nrows     = parse(Int64,split(readline(f))[2])
    xllcorner = parse(Float64,split(readline(f))[2])
    yllcorner = parse(Float64,split(readline(f))[2])
    cellsize  = parse(Float64,split(readline(f))[2])
    nodatval  = parse(Float64,split(readline(f))[2])
    close(f)

    demdat = readdlm(demf,skipstart=6)
    replace!(demdat, -9999=>NaN)

    # GC.gc()
    xdem = collect(xllcorner:cellsize:(xllcorner+cellsize*(ncols-1))) .+ cellsize/2;
    ydem = collect(yllcorner:cellsize:(yllcorner+cellsize*(nrows-1))) .+ cellsize/2;
    tgrid = Matlab.meshgrid(xdem,ydem)

    dem_x = vec(tgrid[1]);
    dem_y = vec(tgrid[2]);
    dem_z = vec(reverse(demdat,dims=1));

    rows = findall(isnan,dem_z)
    deleteat!(dem_x,rows)
    deleteat!(dem_y,rows)
    deleteat!(dem_z,rows)

    return dem_x, dem_y, dem_z, cellsize
end

function loaddbh(fname::String,bounds::Array{Int64,2})
    dat, _ = readdlm(fname,'\t',header=true)

    dat_x, dat_y, dat_z, rmdx = clipdat(dat[:,1],dat[:,2],dat[:,3],bounds,0)

    dat_r = deleteat!(dat[:,4],rmdx) # diameter at breast height
    dat_r /= 2
    dat_r /= 100

    dat_h = deleteat!(dat[:,5],rmdx) # height

    return dat_x, dat_y, dat_z, dat_r
end

function loadltc(fname::String,bounds::Array{Int64,2})
    ltc, _ = readdlm(fname,'\t',header=true)
    replace!(ltc, -9999=>NaN)
    _, _, _, rmdx = clipdat(ltc[:,1],ltc[:,2],ltc[:,3],bounds,0)
    ltc = ltc[setdiff(1:end, rmdx), :]
    return ltc
end


function createfiles(outdir::String,outstr::String,pts::Array{Float64,2},loc_time::Array{DateTime,1},t1::String,t2::String,int::Int64,calc_swr::Bool,append::Bool)

    outfile  = outdir*"/Output_"*outstr*".nc"

    if append
        dataset  = netcdf.Dataset(outfile,"r+")

        Vf_weighted  = dataset.variables["Vf_weighted"]
        Vf_flat      = dataset.variables["Vf_flat"]

        if calc_swr
            swr_tot = dataset.variables["SWR_tot"]
            swr_dir = dataset.variables["SWR_dir"]
            for_tau = dataset.variables["Forest_Transmissivity"]
            return swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset
        else
            return Vf_weighted, Vf_flat, dataset
        end

    else

        writedlm(outdir*"/Time_"*outstr*".txt",Dates.format.(loc_time, "yyyy.mm.dd HH:MM:SS"))

        dataset = netcdf.Dataset(outfile,"w",format="NETCDF4_CLASSIC")

        locxy = dataset.createDimension("loc_XY",size(pts,1))
        locdt = dataset.createDimension("loc_DT",size(loc_time,1))
        ptsn  = dataset.createDimension("ptsn",2)

        Vf_weighted  = dataset.createVariable("Vf_weighted",np.float32,("loc_XY"),zlib="True",
                                            least_significant_digit=3)
        Vf_flat      = dataset.createVariable("Vf_flat",np.float32,("loc_XY"),zlib="True",
                                            least_significant_digit=3)
        Coors        = dataset.createVariable("Coordinates",np.float32,("loc_XY","ptsn"),zlib="TRUE",
                                            least_significant_digit=1)

        for cx in eachindex(pts[:,1])
         Coors[cx] = np.array(pts[cx,:])
        end

        if calc_swr
            swr_tot = dataset.createVariable("SWR_tot",np.float32,("loc_XY","loc_DT"),zlib="True",
                                                least_significant_digit=3)

            swr_dir = dataset.createVariable("SWR_dir",np.float32,("loc_XY","loc_DT"),zlib="True",
                                                least_significant_digit=3)

            for_tau = dataset.createVariable("Forest_Transmissivity",np.float32,("loc_XY","loc_DT"),zlib="True",
                                                least_significant_digit=3)

            return swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset

        else
            return Vf_weighted, Vf_flat, dataset
        end

    end

end


function create_exmat(outdir::String,outstr::String,pts::Array{Float64,2},g_img::Array{Int64,2},append::Bool)

    outfile  = outdir*"/SHIs_"*outstr*".nc"

    if append
        images = netcdf.Dataset(outfile,"r+")
        SHIs   = images.variables["SHI"]
    else

        images   = netcdf.Dataset(outfile,"w",format="NETCDF4_CLASSIC")

        dims     = size(g_img)

        loc1  = images.createDimension("loc1",size(g_img)[1])
        loc2  = images.createDimension("loc2",size(g_img)[2])
        locxy = images.createDimension("locxy",size(pts,1))
        ptsn  = images.createDimension("ptsn",2)

        SHIs  = images.createVariable("SHI",np.int8,("locxy","loc1","loc2"),zlib="TRUE",
                                        least_significant_digit=1)

        Coors = images.createVariable("Coordinates",np.float32,("locxy","ptsn"),zlib="TRUE",
                                        least_significant_digit=1)

        for cx in eachindex(pts[:,1])
         Coors[cx] = np.array(pts[cx,:])
        end
    end

    return SHIs, images

end
