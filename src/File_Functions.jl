function loaddbh(fname::String,limits::Array{Int64,2},peri::Int64)
    dat, _ = readdlm(fname,'\t',header=true)

    dat_x, dat_y, dat_z, rmdx = clipdat(dat[:,1],dat[:,2],dat[:,3],limits,peri)

    dat_r = deleteat!(dat[:,4],rmdx) # diameter at breast height
    dat_r /= 2
    dat_r /= 100

    # dat_h = deleteat!(dat[:,5],rmdx) # height

    if size(dat,2) == 5
        dat_tc = deleteat!(dat[:,5],rmdx)
    else
        dat_tc = fill(NaN,size(dat_x))
    end

    return dat_x, dat_y, dat_z, dat_r, dat_tc
end

function loadltc_txt(fname::String,limits::Array{Int64,2},peri::Int64)
    ltc, _ = readdlm(fname,'\t',header=true)
    replace!(ltc, -9999=>NaN)
    _, _, _, rmdx = clipdat(ltc[:,1],ltc[:,2],ltc[:,3],limits,peri)
    ltc = ltc[setdiff(1:end, findall(rmdx.==1)), :]
    return ltc
end



function loadltc_laz(fname::String,limits::Array{Int64,2},peri::Int64,
            dbh_x::Array{Float64,1},dbh_y::Array{Float64,1},dbh_e::Array{Float64,1},
            lastc::Array{Float64,1})

        # calculate random branch angle for each tree
        ang = rand(Uniform(60,100),size(lastc,1)) # Norway Spruce

        header, ltcdat = LazIO.load(fname)
        ltc = fill(NaN,size(ltcdat,1),8)
        ltc_c = fill(NaN,size(ltcdat,1),1)
        for dx in eachindex(ltcdat)
            ltc_c[dx] = trunc(Int,float(classification(ltcdat[dx])))
            if ltc_c[dx] == 2
                continue
            else
                ltc[dx,1] = xcoord(ltcdat[dx],header) #lasx
                ltc[dx,2] = ycoord(ltcdat[dx],header) #lasy
                ltc[dx,3] = zcoord(ltcdat[dx],header) #lasx
                ltc[dx,7] = Int64.(LasIO.intensity(ltcdat[dx])) # treeclass
                tree_dx = (lastc .== ltc[dx,7])
                if sum(tree_dx) > 0
                    ltc[dx,4] = dbh_x[tree_dx][1] #treex
                    ltc[dx,5] = dbh_y[tree_dx][1] #treey
                    ltc[dx,6] = dist(ltc[dx,1],ltc[dx,2],ltc[dx,4],ltc[dx,5])
                    ltc[dx,8] = ang[tree_dx][1]
                end
            end
        end
        ltc = ltc[setdiff(1:end, findall(isnan.(ltc[:,4]).==1)), :]
        _, _, _, rmdx = clipdat(ltc[:,1],ltc[:,2],ltc[:,3],limits,peri)
        ltc = ltc[setdiff(1:end, findall(rmdx.==1)), :]
    return ltc
end


function createfiles(outdir::String,outstr::String,pts::Array{Float64,2},calc_trans::Bool,calc_swr::Int64,
            append_file::Bool,loc_time=nothing)

    outfile  = outdir*"/Output_"*outstr*".nc"

    if append_file
        dataset  = netcdf.Dataset(outfile,"r+")

        Vf_weighted  = dataset.variables["Vf_weighted"]
        Vf_flat      = dataset.variables["Vf_flat"]

        if calc_trans
            # swr_tot = dataset.variables["SWR_tot"]
            # swr_dir = dataset.variables["SWR_dir"]
            for_tau = dataset.variables["Forest_Transmissivity"]
            return for_tau, Vf_weighted, Vf_flat, dataset
        else
            return Vf_weighted, Vf_flat, dataset
        end

    else

        if calc_trans
            writedlm(outdir*"/Time_"*outstr*".txt",Dates.format.(loc_time, "yyyy.mm.dd HH:MM:SS"))
        end

        dataset = netcdf.Dataset(outfile,"w",format="NETCDF4_CLASSIC")

        locxy = dataset.createDimension("loc_XY",size(pts,1))
        ptsn  = dataset.createDimension("ptsn",2)

        Vf_weighted  = dataset.createVariable("Vf_weighted",np.float32,("loc_XY"),zlib="True",
                                            least_significant_digit=3)
        Vf_flat      = dataset.createVariable("Vf_flat",np.float32,("loc_XY"),zlib="True",
                                            least_significant_digit=3)
        Coors        = dataset.createVariable("Coordinates",np.float32,("loc_XY","ptsn"),zlib="TRUE",
                                            least_significant_digit=1)

        if size(pts,1) == 1
            Coors[1] = np.array(pts[1,1:2])
        else
            for cx in eachindex(pts[:,1])
             Coors[cx] = np.array(pts[cx,:])
            end
        end

        if calc_trans

            locdt = dataset.createDimension("loc_DT",size(loc_time,1))

            for_tau = dataset.createVariable("Forest_Transmissivity",np.float32,("loc_XY","loc_DT"),zlib="True",
                                                least_significant_digit=3)

            if calc_swr > 0

                swr_tot = dataset.createVariable("SWR_tot",np.float32,("loc_XY","loc_DT"),zlib="True",
                                                    least_significant_digit=3)

                swr_dir = dataset.createVariable("SWR_dir",np.float32,("loc_XY","loc_DT"),zlib="True",
                                                    least_significant_digit=3)

                return swr_tot, swr_dir, for_tau, Vf_weighted, Vf_flat, dataset

            else
                return for_tau, Vf_weighted, Vf_flat, dataset
            end
        else
            return Vf_weighted, Vf_flat, dataset
        end

    end

end


function create_exmat(outdir::String,outstr::String,pts::Array{Float64,2},g_img::Array{Int64,2},append_file::Bool)

    outfile  = outdir*"/SHIs_"*outstr*".nc"

    if append_file
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

        if size(pts,1) == 1
            Coors = np.array(pts[1,1:2])
        else
            for cx in eachindex(pts[:,1])
             Coors[cx] = np.array(pts[cx,:])
            end
        end

    end

    return SHIs, images

end
