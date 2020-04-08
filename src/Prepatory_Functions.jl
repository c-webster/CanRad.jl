function extract(d::Dict)
    expr = quote end
    for (k, v) in d
       push!(expr.args, :($(Symbol(k)) = $v))
    end
    # eval(expr)
    return expr
end

function clipdat(pc::Array{Float64,2},pts::Array{Int64,2},peri::Int64)
    kpidx = Int.((minimum(pts[:,1])-peri).<pc[:,1].<(maximum(pts[:,1])+peri)) .*
                Int.((minimum(pts[:,2])-peri).<pc[:,2].<(maximum(pts[:,2])+peri))
    pc = pc[kpidx.==1,:]
    return pc
end

function create_tiles(basefolder::String,ptsf::String)

    tsize = readdlm(basefolder*"/TileSize.txt")

    segs = readdir(basefolder*"/Segments")

    if !ispath(basefolder*"/Tiles")
        mkdir(basefolder*"/Tiles")
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
        mkdir(exdir)
    end

    return ptsfname, inputsegname, exdir

end
