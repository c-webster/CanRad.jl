function findelev(tmpc_x::Array{Float64,1},tmpc_y::Array{Float64,1},tmpc_z::Array{Float64,1},x::Any,y::Any,
                        peri::Int64=10,method::String="linear")
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

function preallo_trunks(dat_x::Array{Float64,1},dat_y::Array{Float64,1},dat_z::Array{Float64,1},
            dat_r::Array{Float64,1},npt::Any,hint::Any)
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

function calculate_trunks(dbh_x::Array{Float64,1},dbh_y::Array{Float64,1},dbh_z::Array{Float64,1},
                dbh_r::Array{Float64,1},npt::Any,hint::Any,e::Any)
    # bh  = 1.5 # breast height
    dims =  preallo_trunks(dbh_x,dbh_y,dbh_z,dbh_r,npt,hint)

    global tsm_x = fill(NaN,Int(sum(dims))); global tsm_y = fill(NaN,Int(sum(dims)))
    global tsm_z = fill(NaN,Int(sum(dims))); global tsm_e = fill(NaN,Int(sum(dims)))

    @inbounds @simd for tix = 1:1:size(dbh_x,1)
        global tsm_x, tsm_y, tsm_z, tsm_e

        if dbh_z[tix] > 2
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

function dist(pcd_x::Array{Float64,1},pcd_y::Array{Float64,1},xcoor::Float64,ycoor::Float64)
            return hypot.(pcd_x.-xcoor,pcd_y.-ycoor)
            # return @fastmath (sqrt.(((dat_x .- xcoor).^2) + ((dat_y .- ycoor).^2)))
end

function dist(dat_x::Float64,dat_y::Float64,xcoor::Float64,ycoor::Float64)
            return hypot.(dat_x.-xcoor,dat_y.-ycoor)
            # return @fastmath (sqrt.(((dat_x .- xcoor).^2) + ((dat_y .- ycoor).^2)))
end

function dist3d(dat_x::Array{Float64,1},dat_y::Array{Float64,1},dat_z::Array{Float64,1},xcoor::Float64,ycoor::Float64,zcoor::Float64)
            return hypot.(hypot.(dat_x.-xcoor,dat_y.-ycoor),dat_z.-zcoor)
end

function dist3d(dat_x::Float64,dat_y::Float64,dat_z::Float64,xcoor::Float64,ycoor::Float64,zcoor::Float64)
            return hypot(hypot.(dat_x.-xcoor,dat_y.-ycoor),dat_z.-zcoor)
end

# for lidar points and branches:
function getsurfdat(dsm_x::Array{Float64,1},dsm_y::Array{Float64,1},dsm_z::Array{Float64,1},
                        bsm_x::Array{Float64,1},bsm_y::Array{Float64,1},bsm_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64)

    dsm_d = dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor).<peri
    bsm_d = dist3d(bsm_x,bsm_y,bsm_z,xcoor,ycoor,ecoor).<peri*0.5

    return append!(dsm_x[dsm_d],bsm_x[bsm_d]), append!(dsm_y[dsm_d],bsm_y[bsm_d]),
                append!(dsm_z[dsm_d],bsm_z[bsm_d])
end

# for lidar/terrain points
function getsurfdat(dsm_x::Array{Float64,1},dsm_y::Array{Float64,1},dsm_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64)
    dsm_d = dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor).<peri
    return dsm_x[dsm_d], dsm_y[dsm_d], dsm_z[dsm_d]
end

# for lidar/terrain points but with a minimum and maximum perimeter
function getsurfdat(dsm_x::Array{Float64,1},dsm_y::Array{Float64,1},dsm_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri1::Int64,peri2::Int64)
    dsm_d = peri1 .< dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor) .< peri2
    return dsm_x[dsm_d], dsm_y[dsm_d], dsm_z[dsm_d]
end


function getsurfdat_lavd(dsm_x::Array{Float64,1},dsm_y::Array{Float64,1},dsm_z::Array{Float64,1},lavd::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64)
    dsm_d = dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor) .< peri
    return lavd[dsm_d]
end

function cart2sph(in_x::Array{Float64,1},in_y::Array{Float64,1},in_z::Array{Float64,1},peri::Int64)
    out_phi    = atan.(in_y,in_x) # az
    out_tht    = ((pi/2) .- (atan.(in_z,hypot.(in_x,in_y)))) * (180/pi) # elev/zenith
    out_tht[out_tht.> 90] .= 90
    out_rad    = hypot.(hypot.(in_x,in_y),in_z) # r

    # return filterbyradius(out_phi,out_tht,out_rad,peri)
    return out_phi, out_tht, out_rad
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

function normalise(pcd_x::Array{Float64,1},pcd_y::Array{Float64,1},pcd_z::Array{Float64,1},
                    xcoor::Float64,ycoor::Float64,ecoor::Float64,ch=0.0::Float64)
    pcd_x .-= xcoor
    pcd_y .-= ycoor
    pcd_z .-= ecoor
    pcd_z .-= ch
    return pcd_x, pcd_y, pcd_z
end

function pcd2pol(pcd_x::Array{Float64,1},pcd_y::Array{Float64,1},pcd_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64,dat::String,
                        ch=0.0::Float64,slp=0.0::Float64,cellsize=0.0::Float64)

    if dat=="terrain"
        pcd_x, pcd_y, pcd_z = getsurfdat(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,peri)
    elseif dat=="chm"
        pcd_x, pcd_y, pcd_z = getsurfdat(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,10,peri)
    elseif dat=="canopy"
        pcd_x, pcd_y, pcd_z = getsurfdat(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,peri)
    end
    pcd_x, pcd_y, pcd_z = normalise(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,ch)

    return cart2sph(pcd_x,pcd_y,pcd_z,peri)

end

function pcd2pol2cart(pcd_x::Array{Float64,1},pcd_y::Array{Float64,1},pcd_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64,dat::String,
                        ch=0.0::Float64,slp=0.0::Float64,cellsize=0.0::Float64)

    pol_phi, pol_tht, pol_rad = pcd2pol(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,peri,dat,ch,slp,cellsize)

    # convert phi/tht to cartesian
    if dat=="terrain"
        pol_phi, pol_tht = calc_horizon_lines(cellsize,peri,pol_phi,pol_tht,pol_rad,slp)
    elseif dat=="chm"
        rows = findall(vec(pol_rad).<(10 .* (sqrt(2)*cellsize))) # deletes the points right above the camera
        deleteat!(pol_phi,rows); deleteat!(pol_tht,rows); deleteat!(pol_rad,rows)
    end

    # convert to cartesian
    pol_phi, pol_tht = pol2cart(pol_phi,pol_tht)

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
    phi_bins = collect(-pi+((pi/360)*3):pi/90:pi-((pi/360)*3))

    idx = collect(1:1:size(pcdpol_phi,1))
    fix1 = Array{Int64,1}(undef,10000)
    tdx = Array{Bool,size(phi_bins,1)}
    mintht = (fill(90.0,size(phi_bins,1),1));

    mintht = calcmintht(mintht,rbins,phi_bins,idx,pcdpol_phi,pcdpol_tht,pcdpol_rad,fix1,copy(mintht),tdx);

    rphi = collect(-pi:pi/1080:pi)
    rtht = Array{Float64,1}(undef,size(rphi,1))

    rtht = LinearInterpolation([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
                        vec([mintht[end];mintht;mintht[1]]))(rphi)

    pol_phi, pol_tht = fillterrain(rphi,rtht,slp)

    return pol_phi, pol_tht
end

function calcCHM_Ptrans(pcd_x::Array{Float64,1},pcd_y::Array{Float64,1},pcd_z::Array{Float64,1},pcd_b::Array{Float64,1},lavd::Array{Float64,1},
                        # stm_x::Array{Float64,1},stm_y::Array{Float64,1},stm_z::Array{Float64,1},
                        xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64,ch::Float64,cellsize::Float64,dat::String="canopy")

    lavd = getsurfdat_lavd(copy(pcd_x),copy(pcd_y),copy(pcd_z),lavd,xcoor,ycoor,ecoor,peri)

    can_phi, can_tht, can_rad = pcd2pol(copy(pcd_x),copy(pcd_y),copy(pcd_z),xcoor,ycoor,ecoor,peri,dat,ch,0.0,cellsize)

    bse_phi, bse_tht, bse_rad = pcd2pol(copy(pcd_x),copy(pcd_y),copy(pcd_b),xcoor,ycoor,ecoor,peri,dat,ch,0.0,cellsize)

    # stm_phi, stm_tht, stm_rad = pcd2pol(copy(stm_x),copy(stm_y),copy(stm_z),xcoor,ycoor,ecoor,peri,dat,ch,0.0,cellsize)

    # allocate variables for the radial loops in calcThickness
    # phi_bins = collect(-pi+((pi/360)*3):pi/360:pi-((pi/360)*3))
  
#     phi_bins = collect(-pi:pi/360:pi)
    phi_bins = [float.(pi);collect(-pi:pi/360:pi);float.(-pi)] # the size of this vector will also have an effect on canopy density in the resulting image

    sum_lavd_thick, sum_thick, rdist = calcThickness(collect(4*cellsize:sqrt(2).*cellsize:peri),phi_bins,
                                            can_phi,can_tht,can_rad,lavd,bse_phi,bse_tht,bse_rad,
                                            Array{Int64,1}(undef,10000),Array{Int64,1}(undef,10000),
                                            collect(1:1:size(can_phi,1)),collect(1:1:size(bse_phi,1)),
                                            Array{Bool,size(phi_bins,1)},Array{Bool,size(phi_bins,1)},
                                            fill(90.0,size(phi_bins,1),1),fill(0.0,size(phi_bins,1),1),fill(0.0,(size(phi_bins,1),90)),
                                            fill(0.0,size(phi_bins,1),90),cellsize,fill(float(peri),(size(phi_bins,1),90)))

    # remove all points where canopy thickness is only one Rbin
    pt_chm_x, pt_chm_y, rdist = calcPtrans(sum_lavd_thick, phi_bins, sum_thick, cellsize, vec(rdist))

    pt_chm_x_thick, pt_chm_y_thick = calcPtrans_dist(fill(0.0,size(phi_bins,1),90), phi_bins, sum_thick)

    return pt_chm_x, pt_chm_y, rdist, pt_chm_x_thick, pt_chm_y_thick
end

function calcThickness(rbins::Array{Float64,1},phi_bins::Array{Float64,1},
                        can_phi::Array{Float64,1},can_tht::Array{Float64,1},can_rad::Array{Float64,1},lavd::Array{Float64,1},
                        bse_phi::Array{Float64,1},bse_tht::Array{Float64,1},bse_rad::Array{Float64,1},
                        fix_can::Array{Int64,1},fix_bse::Array{Int64,1},
                        idx_can::Array{Int64,1},idx_bse::Array{Int64,1},
                        tdx_can,tdx_bse,
                        temp_mintht::Array{Float64,2},thick::Array{Float64,2},sum_thick::Array{Float64,2},
                        sum_lavd_thick::Array{Float64,2},cellsize::Float64,rdist::Array{Float64,2})

    for rbix = length(rbins)-1:-1:1

        # global mintht, minrad, sum_lavd_thick, tempthicks, temp, sum_thick, trunk_temptht, sum_trunk_thick
        fix_can = idx_can[frbins(can_rad,rbins[rbix],rbins[rbix+1])]
        if size(fix_can) > (1,)
            tdx_can = (minimum(can_phi[fix_can[sortperm(can_phi[fix_can])]]) .<= phi_bins
                                .<= maximum(can_phi[fix_can[sortperm(can_phi[fix_can])]]))
            chm_temptht = copy(temp_mintht) # copy(mintht)
            chm_temptht[tdx_can] = LinearInterpolation(can_phi[fix_can[sortperm(can_phi[fix_can])]],
                        can_tht[fix_can[sortperm(can_phi[fix_can])]])(phi_bins[tdx_can])
            if sum(isnan.(chm_temptht)) > 0; chm_temptht[isnan.(chm_temptht)] .= 90.0; end

            temp_lavd = copy(thick) # copy(temp2)
            temp_lavd[tdx_can] = LinearInterpolation(can_phi[fix_can[sortperm(can_phi[fix_can])]],
                        lavd[fix_can[sortperm(can_phi[fix_can])]])(phi_bins[tdx_can])
            tempthick = fill(0.0,(size(chm_temptht,1),90));

            chm_temptht = Int.(round.(chm_temptht));
            len = size(chm_temptht,1)

            #get the canopy base line
            fix_bse = idx_bse[frbins(bse_rad,rbins[rbix],rbins[rbix+1])]
            bse_temptht = copy(temp_mintht) # copy(mintht)
            if size(fix_bse) > (1,)
                tdx_bse = (minimum(bse_phi[fix_bse[sortperm(bse_phi[fix_bse])]]) .<= phi_bins
                                    .<= maximum(bse_phi[fix_bse[sortperm(bse_phi[fix_bse])]]))
                bse_temptht[tdx_bse] = LinearInterpolation(bse_phi[fix_bse[sortperm(bse_phi[fix_bse])]],
                            bse_tht[fix_bse[sortperm(bse_phi[fix_bse])]])(phi_bins[tdx_bse])
                if sum(isnan.(bse_temptht)) > 0; bse_temptht[isnan.(bse_temptht)] .= 90.0; end
            end
            bse_temptht = Int.(round.(bse_temptht));

            for tx = 1:1:len
                # smooth the horizon line
#                 if tx !== 1; lo = chm_temptht[tx-1]; else; lo = chm_temptht[end]; end
#                 if tx !== len; hi = chm_temptht[tx+1]; else; hi = chm_temptht[1]; end
#                 chm_temptht[tx] = Int(round(mean([lo,chm_temptht[tx],hi])))

                # fill the chm stuff
                if 90 - chm_temptht[tx] == 0
                    tempthick[tx,:] .= 0
                else
                    if (bse_temptht[tx])-(chm_temptht[tx]) > 0 # zenith angle
                        if chm_temptht[tx] == 0; chm_temptht[tx] = 1; end
                        tempthick[tx,91-bse_temptht[tx]:(91-chm_temptht[tx])] .= sqrt(2)*cellsize # elevation angle
                        rdist[tx,91-bse_temptht[tx]:(91-chm_temptht[tx])] .= rbins[rbix]
                    else
                        tempthick[tx,1:(90-chm_temptht[tx])] .= sqrt(2)*cellsize
                        rdist[tx,1:(90-chm_temptht[tx])] .= rbins[rbix]
                    end
                end

            end

            sum_lavd_thick = sum_lavd_thick .+ (tempthick .* temp_lavd)

            sum_thick = sum_thick .+ tempthick

        end

    end

    return sum_lavd_thick, sum_thick, vec(rdist)
end

function calcPtrans(sum_lavd_thick::Array{Float64,2},phi_bins::Array{Float64,1},sum_thick::Array{Float64,2},cellsize::Float64,rdist::Array{Float64,1})
    Ptrans = exp.(-0.5 * sum_lavd_thick)
    # phi_bins = collect(-pi:pi/356:pi); #push!(phi_bins,pi) # tidy up phi_bins so it starts and ends at pi

    # solve pt based on average pt and random distribution
#     mean_Ptrans = median(Ptrans[(Ptrans .> 0) .& (Ptrans .< 1)])
    rand_Ptrans = rand(Uniform(0,1),size(Ptrans))

    Ptrans[Ptrans .>= rand_Ptrans] .= 1
    Ptrans[Ptrans .< rand_Ptrans] .= 0
#     Ptrans[sum_thick .<= (2 .* (sqrt(2)*cellsize))] .= 1

    pt_chm_x, pt_chm_y, rdist = getPhiTht(phi_bins,Ptrans,rdist)

    return pt_chm_x, pt_chm_y, rdist
end

function calcPtrans_dist(canopy_thick::Array{Float64,2},phi_bins::Array{Float64,1},sum_thick::Array{Float64,2})

    # zero transmissivity where the canopy is > 15m thick.
    canopy_thick[sum_thick .> 20] .= 0
    canopy_thick[sum_thick .<= 20] .= 1
    phi = repeat(phi_bins,90);
    tht = repeat(90.0:-1:1,inner=size(phi_bins,1))
    rows = findall(vec(canopy_thick).==1)
    deleteat!(phi,rows)
    deleteat!(tht,rows)

    return pol2cart(phi,float.(tht))
end

function getPhiTht(phi_bins::Array{Float64,1},Ptrans::Array{Float64,2},rdist::Array{Float64,1})
    phi = repeat(phi_bins,90);
    tht = repeat(90.0:-1:1,inner=size(phi_bins,1))

    for thtdx = 1:size(phi_bins,1):size(tht,1)-size(phi_bins,1)
        tht[thtdx:thtdx+size(phi_bins,1)-1] = rand((tht[thtdx]-1:tht[thtdx]),size(thtdx:thtdx+size(phi_bins,1)-1))
    end

    temp = sortperm(phi)
    for phidx = 1:size(Ptrans,1):size(phi,1)-size(Ptrans,1)
        phi[temp[phidx:phidx+89]] = rand((phi[temp[phidx]]-0.01:phi[temp[phidx]+1]),size(phidx:phidx+89))
    end

    # remove regions of sky where transmissivity = 0
    rows = findall(vec(Ptrans).==1)
    deleteat!(phi,rows)
    deleteat!(tht,rows)
    deleteat!(rdist,rows)

    pt_chm_x, pt_chm_y = pol2cart(phi,float.(tht))

    return pt_chm_x, pt_chm_y, rdist
end

# function createVariables(pol_phi::Array{Float64,1},phi_bins::Array{Float64,1},peri::Int64)
#     idx = collect(1:1:size(pol_phi,1))
#     fix1 = Array{Int64,1}(undef,10000)
#     tdx = Array{Bool,size(phi_bins,1)}
#     mintht = fill(90.0,size(phi_bins,1),1)
#     minrad = float.(fill(peri,size(phi_bins,1),1))
#     thick  = fill(0.0,size(phi_bins,1),1)
#     sum_lavd_thick = fill(0.0,size(phi_bins,1),90)
#     tempthick = fill(0.0,(size(phi_bins,1),90))
#     rdist = fill(float(peri),(size(phi_bins,1),90))
#     return idx, fix1, tdx, mintht, minrad, thick, sum_lavd_thick, tempthick, rdist
# end



# function trunk_locs(dbh_x::Array{Float64,1},dbh_y::Array{Float64,1},dbh_r::Array{Float64,1})
#
#     t_x =  Vector{Float64}(); t_y =  Vector{Float64}();
#
#     for tx = 1:1:size(dbh_r,1)
#
#         xunit = (dbh_r[tx] * cos.((0:(pi/10):(2*pi))) .+ dbh_x[tx])
#         yunit = (dbh_r[tx] * sin.((0:(pi/10):(2*pi))) .+ dbh_y[tx])
#
#         points = vec(SVector.((minimum(xunit):0.05:maximum(xunit))',(minimum(yunit):0.05:maximum(yunit))))
#
#         polygon = SVector.(xunit,yunit)
#         inside = [inpolygon(p, polygon; in=true, on=true, out=false) for p in points]
#
#         rows = findall(inside.==0)
#         append!(t_x,first.(points)[rows])
#         append!(t_y,last.(points)[rows])
#
#     end
#
#     return t_x, t_y
#
# end
