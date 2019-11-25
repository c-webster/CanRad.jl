function extract(d)
    expr = quote end
    for (k, v) in d
       push!(expr.args, :($(Symbol(k)) = $v))
    end
    eval(expr)
    return
end

function clipdat(pc,pts,peri)
    idxx = (minimum(pts[:,1])-peri).<pc[:,1].<(maximum(pts[:,1])+peri)
    idxy = (minimum(pts[:,2])-peri).<pc[:,2].<(maximum(pts[:,2])+peri)
    kpidx = Int.(idxx) .* Int.(idxy)
    pcn = pc[kpidx.==1,:]
    return pcn
end
