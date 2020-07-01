function findelev(tmpc_x::Array{Float64,1},tmpc_y::Array{Float64,1},tmpc_z::Array{Float64,1},x::Any,y::Any,peri::Int64=10,method::String="linear")
    pc_x, pc_y, pc_z, _ = clipdat(tmpc_x,tmpc_y,tmpc_z,Int.(floor.([x y])),peri)
    v = pyinterp.griddata(hcat(pc_x,pc_y), pc_z, (x, y), method="linear")
    return round.(v,digits=2)
end


function trunkpoints(x1::Float64,y1::Float64,h::Float64,r::Float64,bh::Float64,npt::Any,hint::Any,e::Float64)
    # lower trunk, below breast height, cylinder shape
    th = collect(0:(pi/npt[1]):(2*pi))
    xunit = (r * cos.(th) .+ x1)
    yunit = (r * sin.(th) .+ y1)
    zinc = collect(0.1:hint[1]:bh)

    xnew1 = repeat(xunit,length(zinc))
    ynew1 = repeat(yunit,length(zinc))
    znew1 = repeat(collect(0.1:hint[1]:bh),inner=length(yunit))

    # upper drunk, above breast height, cone shaped
    inttop = Int((floor(h) - bh)/hint[1])
    rchange = collect(range(r,stop=0.001,length=inttop)) ##### CHECK

    xnew2 = vec(rchange' .* cos.(repeat(th',inttop)') .+ x1)
    ynew2 = vec(rchange' .* sin.(repeat(th',inttop)') .+ y1)
    znew2 = repeat(collect((bh+hint[1]):hint[1]:floor(h)),inner=size(th,1))

    enew =  fill(e,((size(znew1,1)+size(znew2,1)),))

    xnew = vcat(xnew1,xnew2)
    ynew = vcat(ynew1,ynew2)
    znew = vcat(znew1,znew2)

    return xnew, ynew, znew, enew
end

function preallo_trunks(dat_x::Array{Float64,1},dat_y::Array{Float64,1},dat_z::Array{Float64,1},dat_r::Array{Float64,1},npt::Any,hint::Any)
    # if size(dat,2) == 1
        # dims = fill(NaN,size(dat,2))
        # dat = dat'
    # else
        dims = fill(NaN,size(dat_x,1))
    # end
    @inbounds @simd for tix in eachindex(dims)

        if size(npt,1) > 1; nix = tix
        else; nix = 1; end

        th = collect(0:pi/npt[nix]:2*pi)
        yunit = (dat_r[tix] * sin.(th) .+ dat_y[tix])
        znew1 = repeat(collect(0.1:hint[nix]:1.5),inner=length(yunit))
        znew2 = repeat(collect((1.5+hint[nix]):hint[nix]:floor(dat_z[tix])),inner=size(th,1))
        global dims[tix,1] = (size(znew1,1) + size(znew2,1))
    end
    return dims
end

function calculate_trunks(dbh_x::Array{Float64,1},dbh_y::Array{Float64,1},dbh_z::Array{Float64,1},dbh_r::Array{Float64,1},npt::Any,hint::Any,e::Any)
    # bh  = 1.5 # breast height
    dims =  preallo_trunks(dbh_x,dbh_y,dbh_z,dbh_r,npt,hint)

    global tsm_x = fill(NaN,Int(sum(dims))); global tsm_y = fill(NaN,Int(sum(dims)))
    global tsm_z = fill(NaN,Int(sum(dims))); global tsm_e = fill(NaN,Int(sum(dims)))

    @inbounds @simd for tix = 1:1:size(dbh_x,1)
        global tsm_x, tsm_y, tsm_z, tsm_e

        if length(npt) > 1
            xnew, ynew, znew, enew = trunkpoints(dbh_x[tix],dbh_y[tix],dbh_z[tix],dbh_r[tix],1.5,Int.(npt[tix]),float.(hint[tix]),e[tix])
        else
            xnew, ynew, znew, enew = trunkpoints(dbh_x[tix],dbh_y[tix],dbh_z[tix],dbh_r[tix],1.5,npt,hint,e[tix])
        end
        if tix == 1
            tsm_x = copy(xnew); tsm_y = copy(ynew)
            tsm_z = copy(znew); tsm_e = copy(enew)
        else
            append!(tsm_x,xnew); append!(tsm_y,ynew)
            append!(tsm_z,znew); append!(tsm_e,enew)
        end
    end
    return tsm_x, tsm_y, tsm_z.+tsm_e
end


function make_branches(ltc::Array{Float64,2})
    spacing = 0.1
    ltc_ad = zeros(size(ltc[:,6],1)) .* NaN

    ltc_ad[ltc[:,8].>90]  = ltc[ltc[:,8].>90,6]./cosd.(90 .- ltc[ltc[:,8].>90,8])
    ltc_ad[ltc[:,8].<90]  = ltc[ltc[:,8].<90,6]./cosd.(ltc[ltc[:,8].<90,8].-90)
    ltc_ad[ltc[:,8].==90] = ltc[ltc[:,8].==90,6]

    ltc_zd = sqrt.(ltc_ad .^2 .- ltc[:,6] .^ 2)

    # calculate height at which branch intersects trunk
    ltc_tz = ltc[:,3] + ltc_zd

    nwmatx = fill!(zeros(Int.(ceil(maximum(ltc_ad./spacing))),size(ltc,1)),NaN)
    nwmaty = fill!(zeros(Int.(ceil(maximum(ltc_ad./spacing))),size(ltc,1)),NaN)
    nwmatz = fill!(zeros(Int.(ceil(maximum(ltc_ad./spacing))),size(ltc,1)),NaN)

    npts = []
    # for bidx = 1:1:size(ltc,1)
    @inbounds for bidx in eachindex(ltc[:,1])
        try
            npts = Int(floor(ltc_ad[bidx,1] / spacing))
            nwmatx[1:npts,bidx] = ltc[bidx,1] .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc[bidx,4] - ltc[bidx,1])
            nwmaty[1:npts,bidx] = ltc[bidx,2] .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc[bidx,5] - ltc[bidx,2])
            nwmatz[1:npts,bidx] = ltc[bidx,3] .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc_tz[bidx,1] - ltc[bidx,3])
        catch
            npts = Int(ceil(ltc_ad[bidx,1] / spacing))
            nwmatx[1:npts,bidx] = ltc[bidx,1] .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc[bidx,4] - ltc[bidx,1])
            nwmaty[1:npts,bidx] = ltc[bidx,2] .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc[bidx,5] - ltc[bidx,2])
            nwmatz[1:npts,bidx] = ltc[bidx,3] .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc_tz[bidx,1] - ltc[bidx,3])
        end
    end

    # bsm = fill(NaN,(size(nwmatx,1)*size(nwmatx,2)),4)
    bsm_x = vec(nwmatx)
    bsm_y = vec(nwmaty)
    bsm_z = vec(nwmatz)
    rows = findall(isnan,bsm_z)
    deleteat!(bsm_x,rows)
    deleteat!(bsm_y,rows)
    deleteat!(bsm_z,rows)
    return bsm_x, bsm_y, bsm_z
end


function create_mat(radius::Int64)
    tgrid = Matlab.meshgrid(collect(1:radius*2),collect(1:radius*2))

    grad  = sqrt.((tgrid[1] .- radius .- 1).^2 .+ (tgrid[2] .- radius .- 1).^2)
    grad[grad .> radius] .= NaN
    gtht  = (grad./radius).*90

    gphi  = reverse(map(atan,(tgrid[2].-radius),(tgrid[1].-radius)),dims=1)

    grid  = fill(1,(size(gphi)))
    # grid[grad .> radius] .= NaN
    gcoors = [vec(gphi) vec(gtht)]
    gcoorcart = float([vec(tgrid[1]) vec(tgrid[2])])

    return grad, gcoors, gcoorcart, grid
end

function dist(dat_x::Array{Float64,1},dat_y::Array{Float64,1},xcoor::Float64,ycoor::Float64)
            return @fastmath (sqrt.(((dat_x .- xcoor).^2) + ((dat_y .- ycoor).^2)))
end

function getsurfdat(dsm_x::Array{Float64,1},dsm_y::Array{Float64,1},dsm_z::Array{Float64,1},
                        bsm_x::Array{Float64,1},bsm_y::Array{Float64,1},bsm_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,peri::Int64)

    dsm_d = dist(dsm_x,dsm_y,xcoor,ycoor).<peri
    bsm_d = dist(bsm_x,bsm_y,xcoor,ycoor).<peri*0.5

    return append!(dsm_x[dsm_d],bsm_x[bsm_d]), append!(dsm_y[dsm_d],bsm_y[bsm_d]),
                append!(dsm_z[dsm_d],bsm_z[bsm_d])
end

function getsurfdat(dsm_x::Array{Float64,1},dsm_y::Array{Float64,1},dsm_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,peri::Int64)
    dsm_d = dist(dsm_x,dsm_y,xcoor,ycoor).<peri
    return dsm_x[dsm_d], dsm_y[dsm_d], dsm_z[dsm_d]
end

function cart2sph(in_x::Array{Float64,1},in_y::Array{Float64,1},in_z::Array{Float64,1},peri::Int64)
    out_phi    = atan.(in_y,in_x) # az
    out_tht    = ((pi/2) .- (atan.(in_z,hypot.(in_x,in_y)))) * (180/pi) # elev/zenith
    out_tht[out_tht.> 90] .= 90
    out_rad    = hypot.(hypot.(in_x,in_y),in_z) # r

    return filterbyradius(out_phi,out_tht,out_rad,peri)
end

function cart2pol(in_x::Array{Float64,1},in_y::Array{Float64,1})
    out_phi    = atan.(in_y,in_x) # az
    out_tht    = ((pi/2) .- (atan.(in_z,hypot.(in_x,in_y)))) * (180/pi) # elev/zenith
    out_tht[out_tht.> 90] .= 90

    return out_phi, out_tht
end


function fillterrain(rphi::Array{Float64,1},rtht::Array{Float64,1},slp::Float64)

    min_rphi = minimum(filter(!isnan,rtht))
    int = size(collect(min_rphi:0.5:(90+slp)),1)
    if int==1; int=2; end
    temp1 = repeat(collect(range(0,stop=1,length=int))',outer=[size(rtht,1),1])
    temp2 = repeat((ones(size(rtht,1),1) .* (90+slp)) - rtht,outer=[1,int])

    tht = vec(repeat(rtht,outer=[1,int]) + temp1 .* temp2)
    phi = vec(repeat(rphi,outer=[1,int]))

    return phi, tht # azimuth, zenith
end


# dist(dat::Array{Float64,2},xcoor::Float64,ycoor::Float64) = @fastmath (sqrt.(((dat[:,1] .- xcoor).^2) +
#                         ((dat[:,2] .- ycoor).^2)))


function normalise(pcd_x::Array{Float64,1},pcd_y::Array{Float64,1},pcd_z::Array{Float64,1},
                    xcoor::Float64,ycoor::Float64,ecoor::Float64,ch=0.0::Float64)
    pcd_x .-= xcoor
    pcd_y .-= ycoor
    pcd_z .-= ecoor
    pcd_z .-= ch
    return pcd_x, pcd_y, pcd_z
end

function pcd2pol2cart(pcd_x::Array{Float64,1},pcd_y::Array{Float64,1},pcd_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64,dat::String,
                        ch=0.0::Float64,slp=0.0::Float64,cellsize=0.0::Float64)

    if dat=="terrain"
        pcd_x, pcd_y, pcd_z = getsurfdat(pcd_x,pcd_y,pcd_z,xcoor,ycoor,peri)
    end
    pcd_x, pcd_y, pcd_z = normalise(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,ch)

    pol_phi, pol_tht, pol_rad = cart2sph(pcd_x,pcd_y,pcd_z,peri)

    # convert phi/tht to cartesian
    if dat=="terrain"
        pol_phi, pol_tht = calc_horizon_lines(cellsize,peri,pol_phi,pol_tht,pol_rad,slp)
    end

    # convert to cartesian
    pol_phi, pol_tht = pol2cart(pol_phi,pol_tht)
    # pol_phi, pol_tht = [pol_tht .* cos.(pol_phi) pol_tht .* sin.(pol_phi)]

    if dat=="terrain"
        return pol_phi, pol_tht # returns cartesian coordinates
    else
        return prepcrtdat(round.(pol_phi,digits=3),round.(pol_tht,digits=3),round.(pol_rad,digits=3))
    end
end

function pol2cart(pol_phi::Array{Float64,1},pol_tht::Array{Float64,1})
    return pol_tht .* cos.(pol_phi), pol_tht .* sin.(pol_phi)
end

function filterbyradius(phi::Array{Float64,1},tht::Array{Float64,1},rad::Array{Float64,1},peri::Int64)
    rmdx = rad .> peri
    return deleteat!(phi,rmdx), deleteat!(tht,rmdx), deleteat!(rad,rmdx)
end

function prepcrtdat(mat_in_x::Array{Float64,1},mat_in_y::Array{Float64,1},mat_in_r::Array{Float64,1})

    idx    = (1:1:size(mat_in_x,1));

    # remove duplicates dat
    sdx = sortperm(mat_in_r);
    mat_in_r = mat_in_r[sdx];
    idx2 = vec(idx[sdx,:]);

    rmdx   = (nonunique(convert(DataFrames.DataFrame,hcat(mat_in_x,mat_in_y))));
    deleteat!(mat_in_x,rmdx)
    deleteat!(mat_in_y,rmdx)
    deleteat!(mat_in_r,rmdx)
    deleteat!(idx2,rmdx)

    rmdx2 = (idx2 .> size(mat_in_x,1))
    return deleteat!(mat_in_x,rmdx2), deleteat!(mat_in_y,rmdx2), deleteat!(mat_in_r,rmdx2)

end

function prepterdat(matcrt_x::Array{Float64,1},matcrt_y::Array{Float64,1})

    matcrt_x = round.(matcrt_x,digits = 1)
    matcrt_y = round.(matcrt_y,digits = 1)

    rmdx   = (nonunique(convert(DataFrames.DataFrame,hcat(matcrt_x,matcrt_y))));

    return deleteat!(matcrt_x,rmdx), deleteat!(matcrt_y,rmdx)

end

function getimagecentre(slp::Float64,asp::Float64)
    if asp == 0 || asp == 360
        X = 0
        Y = slp
    elseif asp == 90
        X = slp
        Y = 0
    elseif asp == 180
        X = 0
        Y = -slp
    elseif asp == 270
        X = -slp
        Y = 0
    elseif asp > 0 && asp < 90
        X = cosd(90-asp) * slp
        Y = cosd(asp) * slp
    elseif asp > 90 && asp < 180
        X = cosd(asp-90) * slp
        Y = -(cosd(180-asp) * slp)
    elseif asp > 180 && asp < 270
        X = -(cosd(270-asp) * slp)
        Y = -(cosd(asp-180) * slp)
    elseif asp > 270 && asp < 360
        X = -(cosd(asp-270) * slp)
        Y = cosd(360-asp) * slp
    end
    return X,Y
end

function findpairs(kdtree::Any,datcrt::Array{Float64,2},tol::Float64,kdtreedims::Int64,knum::Int64)

    distances, indices = scipyspat.cKDTree.query(kdtree,datcrt, k=knum, n_jobs=-1)
    lia = fill(0,kdtreedims,1)
    lia[indices[distances .<= tol],:] .= 1

    return lia
end

function fillmat(kdtree::Any,datcrt::Array{Float64,2},tol::Float64,kdtreedims::Int64,knum::Int64,radius::Int64,mat2ev::Array{Int64,2})
    imdx = reshape(findpairs(kdtree,datcrt,tol,kdtreedims,knum),(radius*2,radius*2))
    mat2ev[imdx.==1] .= 0
    return mat2ev
end

function findmincol(row)
    return findmin(row)[2]
end

function frbins(pcd::Array{Float64,1},rbin1::Float64,rbin2::Float64)
    return ((pcd .>= rbin1) .& (pcd .< rbin2))
end

function calcmintht(mintht::Array{Float64,2},rbins::Array{Float64,1},phi_bins::Array{Float64,1},idx::Array{Int64,1},
                        pcdpol_phi::Array{Float64,1},pcdpol_tht::Array{Float64,1},pcdpol_rad::Array{Float64,1},
                        fix1::Array{Int64,1},temp::Array{Float64,2},tdx)

    @inbounds @simd for rbix = length(rbins)-1:-1:1
            fix1 = idx[frbins(pcdpol_rad,rbins[rbix],rbins[rbix+1])]
            if size(fix1) > (1,)
                tdx = (minimum(pcdpol_phi[fix1[sortperm(pcdpol_phi[fix1])]]) .<= phi_bins
                                    .<= maximum(pcdpol_phi[fix1[sortperm(pcdpol_phi[fix1])]]))
                temp2 = copy(temp)
                temp2[tdx] = LinearInterpolation(pcdpol_phi[fix1[sortperm(pcdpol_phi[fix1])]],
                            pcdpol_tht[fix1[sortperm(pcdpol_phi[fix1])]])(phi_bins[tdx])
                mintht = minimum(hcat(mintht,temp2),dims=2)
            end
    end
    return mintht

end

function calc_horizon_lines(cellsize::Float64,peri::Int64,pcdpol_phi::Array{Float64,1},
                            pcdpol_tht::Array{Float64,1},pcdpol_rad::Array{Float64,1},slp::Float64)

    rbins = collect(2*cellsize:sqrt(2).*cellsize:peri)
    phi_bins = collect(-pi+((pi/360)*3):pi/720:pi-((pi/360)*3))

    idx = collect(1:1:size(pcdpol_phi,1))
    fix1 = Array{Int64,1}(undef,10000)
    tdx = Array{Bool,size(phi_bins,1)}
    mintht = float.(fill(90,size(phi_bins,1),1));

    mintht = calcmintht(mintht,rbins,phi_bins,idx,pcdpol_phi,pcdpol_tht,pcdpol_rad,fix1,copy(mintht),tdx);

    rphi = collect(-pi:pi/1080:pi)
    rtht = Array{Float64,1}(undef,size(rphi,1))

    rtht = LinearInterpolation([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
                        vec([mintht[end];mintht;mintht[1]]))(rphi)

    pol_phi, pol_tht = fillterrain(rphi,rtht,slp)

    return pol_phi, pol_tht
end




#
# function calcThickness(pcd,xcoor,ycoor,ecoor,peri1,peri2,cellsize)
#
#     pcdpol = pcd2pol(pcd,xcoor,ycoor,ecoor,peri2)
#
#     rbins = collect(peri1:sqrt(2).*cellsize:peri2)
#
#     phi_bins = collect(-pi+((pi/360)*3):pi/720:pi-((pi/360)*3)) # azimuth bins
#
#     global mintht = ones(size(phi_bins,1)) .* 90
#     global minrad = ones(size(phi_bins,1)) .* peri2
#     global thick  = zeros(size(temptht,1),90)
#
#     for rbix = length(rbins)-1:-1:1
#         global mintht, minrad
#             fix1 = findall((pcdpol[:,3] .>= rbins[rbix]) .& (pcdpol[:,3] .< rbins[rbix+1]))
#             if size(fix1) > (1,)
#                 fix2 = sortperm(pcdpol[fix1,1])
#                 tdx = findall(minimum(pcdpol[fix1[fix2],1]) .<= phi_bins
#                                     .<= maximum(pcdpol[fix1[fix2],1]))
#
#                 temptht = ones(size(phi_bins,1)) .* 90
#                 temptht[tdx] =  pyinterp.interp1d(pcdpol[fix1[fix2],1],
#                             pcdpol[fix1[fix2],2])(phi_bins[tdx])
#
#                 flag = vec(mintht .> temptht)
#                 minrad[flag.==1,:] .= mean([rbins[rbix],rbins[rbix+1]])
#                 mintht = minimum(hcat(mintht,temptht),dims=2)
#
#                 tempthick = zeros(size(temptht,1),90)
#                 temptht = Int.(round.(temptht))
#                 for tx = 1:1:size(temptht,1)
#                     if temptht[tx] == 90
#                         continue
#                     else
#                         tempthick[tx,1:temptht[tx]] .= cellsize
#                     end
#                 end
#
#                 thick = thick .+ tempthick
#
#             end
#     end
#
#     rphi = collect(-pi:pi/360:pi)
#     rtht = pyinterp.interp1d([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
#                         vec([mintht[end];mintht;mintht[1]]),fill_value="extrapolate")(rphi)
#     rrad = pyinterp.interp1d([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
#                         vec([minrad[end];minrad;minrad[1]]),fill_value="extrapolate")(rphi)
#
#     return rtht,rrad,thick
#
# end
