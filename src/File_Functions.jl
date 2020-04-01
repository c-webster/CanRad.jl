extension(url::String) = match(r"\.[A-Za-z0-9]+$", url).match

function readlas(infile,ground=nothing)
    if extension(infile) == ".laz"
        header, dsmdat = LazIO.load(infile)
    elseif extension(infile) == ".las"
        header, dsmdat = FileIO.load(infile)
    else
        error("Unknown DSM file extension")
    end

    dsm = zeros(size(dsmdat,1),4) .* NaN
    dsm_c = zeros(size(dsmdat,1)) .* NaN
    for dx in eachindex(dsmdat)
        dsm_c[dx] = trunc(Int,float(classification(dsmdat[dx])))
        if dsm_c[dx] == 2
            continue
        else
            dsm[dx,1] = xcoord(dsmdat[dx],header)
            dsm[dx,2] = ycoord(dsmdat[dx],header)
            dsm[dx,3] = zcoord(dsmdat[dx],header)
        end
    end
    rows = findall(isnan,dsm[:,1])
    dsm = dsm[setdiff(1:end, rows), :]
    return dsm
end

function importdtm(dtmf,tilt)
    if tilt
        file = matopen(dtmf); dtmdat = read(file,"dtm"); close(file)
        dtm1 = dtmdat["x"]
        dtm2 = dtmdat["y"]
        dtm3 = dtmdat["z"]
        dtm4 = dtmdat["s"]
        dtm5 = dtmdat["a"]
        dtm = [dtm1 dtm2 dtm3 dtm4 dtm5]
        dtm_cellsize = dtmdat["cellsize"]
    else
        if extension(dtmf) == ".mat"
            file = matopen(dtmf); dtmdat = read(file,"dtm"); close(file)
            dtm1 = dtmdat["x"]
            dtm2 = dtmdat["y"]
            dtm3 = dtmdat["z"]
            dtm = [dtm1 dtm2 dtm3]
            dtm_cellsize = dtmdat["cellsize"]
        elseif extension(dtmf) == ".asc"
            dtm, cellsize = read_ascii(dtmf)
        end
    end
    return dtm, dtm_cellsize
end

function read_ascii(demf)

    f = open(demf)
    ncols     = parse(Float64,split(readline(f))[2])
    nrows     = parse(Float64,split(readline(f))[2])
    xllcorner = parse(Float64,split(readline(f))[2])
    yllcorner = parse(Float64,split(readline(f))[2])
    cellsize  = parse(Float64,split(readline(f))[2])
    nodatval  = parse(Float64,split(readline(f))[2])
    close(f)

    demdat = readdlm(demf,skipstart=6)
    replace!(demdat, -9999=>NaN)

    xdem  = collect(xllcorner:cellsize:(xllcorner+cellsize*(ncols-1))) .+ cellsize/2
    ydem  = collect(yllcorner:cellsize:(yllcorner+cellsize*(nrows-1))) .+ cellsize/2
    tgrid = Matlab.meshgrid(xdem,ydem)

    dem      = zeros(size(demdat,1) * size(demdat,2),3) .* NaN
    dem[:,1] = vec(tgrid[1])
    dem[:,2] = vec(tgrid[2])
    dem[:,3] = vec(reverse(demdat,dims=1))

    return dem, cellsize
end

function createfiles(outdir,pts,loc_time,t1,t2,int,calc_swr)
    # writedlm(outdir*"/Coords_"*string(Int(pts[1,1]))*"_"*string(Int(pts[1,2]))*".txt",pts)

    file = matopen(outdir*"/Time_"*string(Int(pts[1,1]))*"_"*string(Int(pts[1,2]))*".mat", "w")
    write(file,"coords",pts)
    write(file,"T1",t1)
    write(file,"int",int)
    write(file,"T2",t2)
    close(file)

    outfile  = outdir*"/Output_"*string(Int(pts[1,1]))*"_"*string(Int(pts[1,2]))*".nc"

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


function create_exmat(outdir,pts,g_img)

    outfile  = outdir*"/SHIs_"*string(Int(pts[1,1]))*"_"*string(Int(pts[1,2]))*".nc"

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

    return SHIs, images

end
