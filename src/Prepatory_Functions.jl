function extract(d::Dict)
    expr = quote end
    for (k, v) in d
       push!(expr.args, :($(Symbol(k)) = $v))
    end
    # eval(expr)
    return expr
end

function clipdat(pc_x::Array{Float64,1},pc_y::Array{Float64,1},pc_z::Array{Float64,1},pts::Array{Int64,2},peri::Int64)
    rmidx = (pc_x.<(minimum(pts[:,1])-peri)) .| (pc_x.>(maximum(pts[:,1])+peri)) .|
                (pc_y.<(minimum(pts[:,2])-peri)) .| (pc_y.>(maximum(pts[:,2])+peri))
    deleteat!(pc_x,rmidx)
    deleteat!(pc_y,rmidx)
    deleteat!(pc_z,rmidx)
    return pc_x, pc_y, pc_z, rmidx
end

function create_tiles(basefolder::String,ptsf::String)

    tsize = readdlm(basefolder*"/TileSize.txt")

    segs = readdir(basefolder*"/Segments")

    if !ispath(basefolder*"/Tiles")
        try
            mkdir(basefolder*"/Tiles"); catch
        end
    end

    pts_all = Int.(readdlm(ptsf))

    tdx = 0

    ptsfname     = String[]
    inputsegname = String[]

    for sdx in eachindex(segs)

        limits = readdlm(basefolder*"/Segments"*"/"*segs[sdx]*"/"*segs[sdx]*"_analysisarea.txt")

        limx = Int.(collect(limits[1]:tsize[1]:limits[2]))
        limy = Int.(collect(limits[3]:tsize[1]:limits[4]))

        for x in eachindex(limx[1:end-1]), y in eachindex(limy[1:end-1])
            # create the tasks
            idx = (limx[x] .<= pts_all[:,1] .< limx[x+1]) .&
                    (limy[y] .<= pts_all[:,2] .< limy[y+1])

            if sum(idx) > 0
                tdx += 1
                push!(ptsfname,sprintf1.("%03.$(0)f", tdx)*".txt")
                writedlm(basefolder*"/Tiles"*"/"*ptsfname[end],Int.(pts_all[idx,:]))

                push!(inputsegname,segs[sdx])
            end
        end

    end

    exdir = basefolder*"/Output/"
    if !ispath(exdir)
        try
            mkdir(exdir); catch
        end
    end

    return ptsfname, inputsegname, exdir

end

function check_output(exdir,pts,batch)

    if batch
        outstr = string(Int(floor(pts[1,1])))*"_"*string(Int(floor(pts[1,2])))
        outdir = exdir*"/"*outstr
    else
        outstr = String(split(exdir,"/")[end])
        outdir = exdir
    end

    if ispath(outdir)
        try
            crxstart = parse(Int,(split(reduce(1,readdir(outdir)[findall(occursin.("Processing",readdir(outdir)))])))[4])
            if crxstart == size(pts,1)
                out = true
            else
                out = false
            end
        catch
            out = false
        end
    else
        out = false
    end
    return out
end
