function findelev(tmpc::Array{Float64,2},x::Any,y::Any)
    tmpc_clip = clipdat(tmpc,Int.(floor.([x y])),10)
    v = pyinterp.griddata(tmpc_clip[:,1:2], tmpc_clip[:,3], (x, y), method="linear")
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

function preallo_trunks(dat::Array{Float64,2},npt::Any,hint::Any)
    if size(dat,2) == 1
        dims = fill(NaN,size(dat,2))
        dat = dat'
    else
        dims = fill(NaN,size(dat,1))
    end
    for tix in eachindex(dims)
        x1 = dat[tix,1]
        y1 = dat[tix,2]
        h  = dat[tix,3]
        r  = dat[tix,4]

        if size(npt,1) > 1
            nix = tix
        else
            nix = 1
        end

        th = collect(0:pi/npt[nix]:2*pi)
        yunit = (r * sin.(th) .+ y1)
        znew1 = repeat(collect(0.1:hint[nix]:1.5),inner=length(yunit))
        znew2 = repeat(collect((1.5+hint[nix]):hint[nix]:floor(h)),inner=size(th,1))
        global dims[tix,1] = (size(znew1,1) + size(znew2,1))
    end
    return dims
end

function calculate_trunks(dbh1::Array{Float64,2},npt::Any,hint::Any,e::Any)
    bh  = 1.5 # breast height
    dims =  preallo_trunks(dbh1,npt,hint)

    if size(dbh1,2) == 1
        endix = 1
    else
        endix = size(dbh1,1)
    end

    @inbounds for tix = 1:1:endix
        if tix == 1
            global outtsm = fill(NaN,Int(sum(dims)),4)
        end

        if size(dbh1,2) == 1
            x1 = dbh1[1]; y1 = dbh1[2]
            h  = dbh1[3]; r  = dbh1[4]
        else
            x1 = dbh1[tix,1]; y1 = dbh1[tix,2]
            h  = dbh1[tix,3]; r  = dbh1[tix,4]
        end

        if length(npt) > 1
            xnew, ynew, znew, enew = trunkpoints(x1,y1,h,r,bh,Int.(npt[tix]),Float.(hint[tix]),e[tix])
        else
            xnew, ynew, znew, enew = trunkpoints(x1,y1,h,r,bh,npt,hint,e[tix])
        end

        dim = Int(size(xnew,1))

        if tix == 1
            global sdim = 1
            global edim = sdim+dim
            outtsm[1:edim-1,:] = [xnew ynew znew enew]
        else
            outtsm[edim:edim+dim-1,:] = [xnew ynew znew enew]
            edim += dim
        end

    end

    return outtsm
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

    bsm = fill(NaN,(size(nwmatx,1)*size(nwmatx,2)),4)
    bsm[:,1] = vec(nwmatx)
    bsm[:,2] = vec(nwmaty)
    bsm[:,3] = vec(nwmatz)
    rows = findall(isnan,bsm[:,1])
    bsm = bsm[setdiff(1:end, rows), :]

    return bsm
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

function getsurfdat(dsm::Array{Float64,2},bsm::Array{Float64,2},xcoor::Float64,ycoor::Float64,peri::Int64)
    vcat(dsm[dist(dsm,xcoor,ycoor).<peri,:],bsm[dist(bsm,xcoor,ycoor).<peri*0.5,:])
end

function cart2sph(in::Array{Float64,2})
    out = Array{Float64,2}(undef,size(in,1),3)
    out[:,1]    = atan.(in[:,2],in[:,1]) # az
    out[:,2]    = atan.(in[:,3],hypot.(in[:,1],in[:,2])) # elev
    out[:,3]    = hypot.(hypot.(in[:,1],in[:,2]),in[:,3]) # r
    return out
end

function dem2pol(dem::Array{Float64,2},pts_x::Float64,pts_y::Float64,ch::Float64,peri::Int64,dem_cellsize::Float64,slp::Float64)

    if size(pts_x,1) > 1
        loc_x = mean([maximum(filter(!isnan,pts[:,1])),minimum(filter(!isnan,pts[:,1]))])
        loc_y = mean([maximum(filter(!isnan,pts[:,2])),minimum(filter(!isnan,pts[:,2]))])
    else
        loc_x = pts_x
        loc_y = pts_y
    end

    dem = dem[dist(dem,loc_x,loc_y).< peri,:]
    dem[:,3] = dem[:,3] .- (pyinterp.griddata(dem[:,1:2], dem[:,3], (loc_x, loc_y), method="linear")) .- ch;
    dem[:,1] = dem[:,1] .- loc_x;
    dem[:,2] = dem[:,2] .- loc_y;

    dem = cart2sph(dem)

    dem = correct_sph(dem,peri);

    # calculate horizon lines
    dem = calc_horizon_lines(dem_cellsize*2,peri,dem,slp)

    @fastmath dem[:,1:2] = [dem[:,2] .* cos.(dem[:,1]) dem[:,2] .* sin.(dem[:,1])]

    return dem
end

function fillterrain(rphi::Array{Float64,1},rtht::Array{Float64,1},slp::Float64)

    min_rphi = minimum(filter(!isnan,rtht))
    int = size(collect(min_rphi:0.5:(90+slp)),1)
    temp1 = repeat(collect(range(0,stop=1,length=int))',outer=[size(rtht,1),1])
    temp2 = repeat((ones(size(rtht,1),1) .* (90+slp)) - rtht,outer=[1,int])

    tht = vec(repeat(rtht,outer=[1,int]) + temp1 .* temp2)
    phi = vec(repeat(rphi,outer=[1,int]))

    return hcat(phi,tht) # azimuth, zenith
end

dist(dat::Array{Float64,2},xcoor::Float64,ycoor::Float64) = @fastmath (sqrt.(((dat[:,1] .- xcoor).^2) +
                        ((dat[:,2] .- ycoor).^2)))


function normalise(pcd::Array{Float64,2},xcoor::Float64,ycoor::Float64,ecoor::Float64,ch=0.0::Float64)
    pcd[:,1] .-= xcoor
    pcd[:,2] .-= ycoor
    pcd[:,3] = pcd[:,3] .- ecoor .- ch
    return pcd
end

function correct_sph(pol::Array{Float64,2},peri::Int64)
    @fastmath @views pol[:,2] = ((pi/2) .- pol[:,2]) * (180/pi) # zenith (theta)
    pol[pol[:,2] .> 90,2] .= 90
    pol = pol[pol[:,3] .< peri,:]
    # return pol
end

function pcd2pol2cart(pcd::Array{Float64,2},xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64,dat::String,ch=0.0::Float64,slp=0.0::Float64,cellsize=0.0::Float64)

    if dat=="terrain"
        pcd = pcd[dist(pcd,xcoor,ycoor).< peri,:]
    end
    pcd = normalise(pcd,xcoor,ycoor,ecoor,ch)

    # pol = cart2sph(pcd);
    pol = correct_sph(cart2sph(pcd),peri);

    # convert phi/tht to cartesian
    if dat=="terrain"
        pol = calc_horizon_lines(cellsize,peri,pol,slp)
    end

    @fastmath @views pol[:,1:2] = [pol[:,2] .* cos.(pol[:,1]) pol[:,2] .* sin.(pol[:,1])]

    return pol
end


function prepcrtdat(mat_in::Array{Float64,2},res::Int64)

    mat_in = round.(mat_in,digits=res)
    idx    = collect(1:1:size(mat_in,1))

    # remove duplicates dat
    sdx = Array{Int8,1}(undef,size(mat_in[:,3],1))
    sdx = sortperm(mat_in[:,3])
    mat_in[:,3] = mat_in[sdx,3]
    idx2 = Array{Int8,1}(undef,size(mat_in,1))
    idx2 = vec(idx[sdx,:])

    kpdx   = Array{Int8,1}(undef,size(mat_in,1));
    kpdx   = findall(.!nonunique(convert(DataFrames.DataFrame,mat_in[:,1:2])));
    mat_in = mat_in[kpdx,:]
    idx2   = vec(idx2[kpdx,:])

    mat_in = mat_in[findall(idx2 .<= size(mat_in,1)),:]

    return mat_in
end

function prepterdat(matcrt::Array{Float64,2})

    matcrt = round.(matcrt,digits = 1)
    matcrt = matcrt[findall(.!nonunique(convert(DataFrames.DataFrame,matcrt))),:]

    return matcrt
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
                            pcdpol::Array{Float64,2},fix1::Array{Int64,1},temp::Array{Float64,2},tdx)

    for rbix = length(rbins)-1:-1:1
            fix1 = idx[frbins(pcdpol[:,3],rbins[rbix],rbins[rbix+1])]
            if size(fix1) > (1,)
                tdx = (minimum(pcdpol[fix1[sortperm(pcdpol[fix1,1])],1]) .<= phi_bins
                                    .<= maximum(pcdpol[fix1[sortperm(pcdpol[fix1,1])],1]))
                temp2 = copy(temp)
                temp2[tdx] = LinearInterpolation(pcdpol[fix1[sortperm(pcdpol[fix1,1])],1],
                            pcdpol[fix1[sortperm(pcdpol[fix1,1])],2])(phi_bins[tdx])
                mintht = minimum(hcat(mintht,temp2),dims=2)
            end
    end
    return mintht

end


function calc_horizon_lines(cellsize::Float64,peri::Int64,pcdpol::Array{Float64,2},slp::Float64)

    rbins = collect(2*cellsize:sqrt(2).*cellsize:peri)
    phi_bins = collect(-pi+((pi/360)*3):pi/720:pi-((pi/360)*3))

    idx = collect(1:1:size(pcdpol,1))
    fix1 = Array{Int64,1}(undef,10000)
    tdx = Array{Bool,size(phi_bins,1)}
    mintht = float.(fill(90,size(phi_bins,1),1));
    temp = copy(mintht);

    mintht = calcmintht(mintht,rbins,phi_bins,idx,pcdpol,fix1,temp,tdx);

    rphi = collect(-pi:pi/1080:pi)
    rtht = Array{Float64,1}(undef,size(rphi,1))

    rtht = LinearInterpolation([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
                        vec([mintht[end];mintht;mintht[1]]))(rphi)

    pol = fillterrain(rphi,rtht,slp)

    return pol
end

function calcThickness(pcd,xcoor,ycoor,ecoor,peri1,peri2,cellsize)

    pcdpol = pcd2pol(pcd,xcoor,ycoor,ecoor,peri2)

    rbins = collect(peri1:sqrt(2).*cellsize:peri2)

    phi_bins = collect(-pi+((pi/360)*3):pi/720:pi-((pi/360)*3)) # azimuth bins

    global mintht = ones(size(phi_bins,1)) .* 90
    global minrad = ones(size(phi_bins,1)) .* peri2
    global thick  = zeros(size(temptht,1),90)

    for rbix = length(rbins)-1:-1:1
        global mintht, minrad
            fix1 = findall((pcdpol[:,3] .>= rbins[rbix]) .& (pcdpol[:,3] .< rbins[rbix+1]))
            if size(fix1) > (1,)
                fix2 = sortperm(pcdpol[fix1,1])
                tdx = findall(minimum(pcdpol[fix1[fix2],1]) .<= phi_bins
                                    .<= maximum(pcdpol[fix1[fix2],1]))

                temptht = ones(size(phi_bins,1)) .* 90
                temptht[tdx] =  pyinterp.interp1d(pcdpol[fix1[fix2],1],
                            pcdpol[fix1[fix2],2])(phi_bins[tdx])

                flag = vec(mintht .> temptht)
                minrad[flag.==1,:] .= mean([rbins[rbix],rbins[rbix+1]])
                mintht = minimum(hcat(mintht,temptht),dims=2)

                tempthick = zeros(size(temptht,1),90)
                temptht = Int.(round.(temptht))
                for tx = 1:1:size(temptht,1)
                    if temptht[tx] == 90
                        continue
                    else
                        tempthick[tx,1:temptht[tx]] .= cellsize
                    end
                end

                thick = thick .+ tempthick

            end
    end

    rphi = collect(-pi:pi/360:pi)
    rtht = pyinterp.interp1d([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
                        vec([mintht[end];mintht;mintht[1]]),fill_value="extrapolate")(rphi)
    rrad = pyinterp.interp1d([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
                        vec([minrad[end];minrad;minrad[1]]),fill_value="extrapolate")(rphi)

    return rtht,rrad,thick

end
