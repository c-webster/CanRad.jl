function readlas(infile,lastoolspath)
    if infile[end-3:end] == ".laz"
        header, dsmdat = LazIO.load(infile)
    elseif infile[end-3:end] == ".las"
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
        file = matopen(dtmf); dtmdat = read(file,"dtm"); close(file)
        dtm1 = dtmdat["x"]
        dtm2 = dtmdat["y"]
        dtm3 = dtmdat["z"]
        dtm = [dtm1 dtm2 dtm3]
        dtm_cellsize = dtmdat["cellsize"]
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

function createfiles(outdir,pts,loc_time,t1,t2,int)
    writedlm(outdir*"/Coords_"*string(Int(pts[1,1]))*"_"*string(Int(pts[1,2]))*".txt",pts)

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

    swr_tot = dataset.createVariable("SWR_tot",np.float32,("loc_XY","loc_DT"),zlib="True",
                                        least_significant_digit=3)

    swr_dir = dataset.createVariable("SWR_dir",np.float32,("loc_XY","loc_DT"),zlib="True",
                                        least_significant_digit=3)

    for_tau = dataset.createVariable("Forest_Transmissivity",np.float32,("loc_XY","loc_DT"),zlib="True",
                                        least_significant_digit=3)

    Vf_weighted  = dataset.createVariable("Vf_weighted",np.float32,("loc_XY"),zlib="True",
                                        least_significant_digit=3)
    Vf_flat      = dataset.createVariable("Vf_flat",np.float32,("loc_XY"),zlib="True",
                                        least_significant_digit=3)

    return swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset
end
