function findelev(tmpc,x,y)
    tmpc_clip = clipdat(tmpc,[x y],10)
    v = pyinterp.griddata(tmpc_clip[:,1:2], tmpc_clip[:,3], (x, y), method="linear")
    return v
end

function trunkpoints(x1,y1,h,r,bh,npt,hint,e)
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

function preallo_trunks(dat,npt,hint)
    if size(dat,2) == 1
        dims = zeros(size(dat,2)) .* NaN
        dat = dat'
    else
        dims = zeros(size(dat,1)) .* NaN
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

function calculate_trunks(dbh1,npt,hint,e)
    bh  = 1.5 # breast height
    dims =  preallo_trunks(dbh1,npt,hint)

    if size(dbh1,2) == 1
        endix = 1
    else
        endix = size(dbh1,1)
    end

    @inbounds for tix = 1:1:endix
        if tix == 1
            global outtsm = zeros(Int(sum(dims)),4) .* NaN
        end

        if size(dbh1,2) == 1
            x1 = dbh1[1]; y1 = dbh1[2]
            h  = dbh1[3]; r  = dbh1[4]
        else
            x1 = dbh1[tix,1]; y1 = dbh1[tix,2]
            h  = dbh1[tix,3]; r  = dbh1[tix,4]
        end

        if length(npt) > 1
            xnew, ynew, znew, enew = trunkpoints(x1,y1,h,r,bh,npt[tix],hint[tix],e[tix])
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

function make_branches(ltc)
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

    bsm = zeros(size(nwmatx,1)*size(nwmatx,2),4) .* NaN
    bsm[:,1] = vec(nwmatx)
    bsm[:,2] = vec(nwmaty)
    bsm[:,3] = vec(nwmatz)
    rows = findall(isnan,bsm[:,1])
    bsm = bsm[setdiff(1:end, rows), :]

    return bsm
end

function create_mat(radius)
    tgrid = Matlab.meshgrid(collect(1:radius*2),collect(1:radius*2))

    grad  = sqrt.((tgrid[1] .- radius .- 1).^2 .+ (tgrid[2] .- radius .- 1).^2)
    grad[grad .> radius] .= NaN
    gtht  = (grad./radius).*90

    gphi  = reverse(map(atan,(tgrid[2].-radius),(tgrid[1].-radius)),dims=1)

    grid  = ones(size(gphi))
    grid[grad .> radius] .= NaN
    gcoors = [vec(gphi) vec(gtht)]
    gcoorcart = float([vec(tgrid[1]) vec(tgrid[2])])

    return grad, gcoors, gcoorcart, grid
end

function cart2sph(x,y,z)
    hypotxy = hypot.(x,y)
    r       = hypot.(hypotxy,z)
    elev    = atan.(z,hypotxy)
    az      = atan.(y,x)
    return (az,elev,r)
end

function dem2pol(dem,pts_x,pts_y,ch,peri,dem_cellsize,slp)

    if size(pts_x,1) > 1
        loc_x = mean([maximum(filter(!isnan,pts[:,1])),minimum(filter(!isnan,pts[:,1]))])
        loc_y = mean([maximum(filter(!isnan,pts[:,2])),minimum(filter(!isnan,pts[:,2]))])
    else
        loc_x = pts_x
        loc_y = pts_y
    end

    xpc = dem[:,1] .- loc_x
    ypc = dem[:,2] .- loc_y
    zpc = dem[:,3] .- (pyinterp.griddata(dem[:,1:2], dem[:,3], (loc_x, loc_y), method="linear")) .- ch

    hdist = sqrt.(((dem[:,1] .- loc_x).^2) + ((dem[:,2] .- loc_y).^2))

    xpc[hdist .> peri] .= NaN

    rows = findall(isnan,xpc)
    xpc = xpc[setdiff(1:end, rows), :]
    ypc = ypc[setdiff(1:end, rows), :]
    zpc = zpc[setdiff(1:end, rows), :]

    temppol = cart2sph(xpc,ypc,zpc)

    dpol = zeros(size(xpc,1),3) .* NaN
    dpol[:,1] = temppol[1] # azimuth (phi)
    dpol[:,2] = (pi/2 .- temppol[2]) .*(180/pi) # zenith (theta)
    # dpol[dpol[:,2] .> 90,2] .= 90
    dpol[:,3] = temppol[3] # rad
    dpol[dpol[:,3] .>  peri,3] .= NaN

    # remove nans
    kpidx = findall(.!isnan.(dpol[:,3]))
    dpol   = (dpol[kpidx,:])
    kpidx = findall(.!isnan.(dpol[:,2]))
    dpol   = (dpol[kpidx,:])

    # calculate horizon lines
    dempol = calc_horizon_lines(dem_cellsize,peri,dpol,slp)

    return dempol
end

function fillterrain(rphi,rtht,slp)

    min_rphi = minimum(filter(!isnan,rtht))
    int = size(collect(min_rphi:0.5:(90+slp)),1)
    temp1 = repeat(collect(range(0,stop=1,length=int))',outer=[size(rtht,1),1])
    temp2 = repeat((ones(size(rtht,1),1) .* (90+slp)) - rtht,outer=[1,int])

    tht = vec(repeat(rtht,outer=[1,int]) + temp1 .* temp2)
    phi = vec(repeat(rphi,outer=[1,int]))

    return hcat(phi,tht) # azimuth, zenith
end

function pcd2pol(pcd,xcoor,ycoor,ecoor,peri,ch=0)

    pcd = clipdat(pcd,[xcoor ycoor],peri)
    plc = zeros(size(pcd,1),3) .* NaN

    plc[:,1] = pcd[:,1] .- xcoor
    plc[:,2] = pcd[:,2] .- ycoor
    if size(pcd,2) == 4
        plc[:,3] = pcd[:,3] .+ pcd[:,4] .- ecoor .- ch
    else
        plc[:,3] = pcd[:,3] .- ecoor .- ch
    end

    @fastmath kpidx = (sqrt.(((pcd[:,1] .- xcoor).^2) +
                        ((pcd[:,2] .- ycoor).^2))) .< peri
    plcn = copy(plc[kpidx.==1,:])

    temppol = cart2sph(plcn[:,1],plcn[:,2],plcn[:,3])

    pol = zeros(size(plcn,1),3) .* NaN
    pol[:,1] = temppol[1] # azimuth (phi)
    pol[:,2] = ((pi/2) .- temppol[2]) * (180/pi) # zenith (theta)
    # pol[pol[:,2] .> 90,2] .= NaN
    pol[:,3] = temppol[3] # rad
    pol[pol[:,3] .>  peri,3] .= NaN

    kpidx = findall(.!isnan.(pol[:,3]))
    pol   = (pol[kpidx,:])
    kpidx = findall(.!isnan.(pol[:,2]))
    pol   = (pol[kpidx,:])

    return pol
end

function pol2cart(phi,tht)
    @fastmath datcrt = [tht .* cos.(phi) tht .* sin.(phi)]
    return datcrt
end

function prepcrtdat(matcrt,matrad,res)

    matcrt = round.(matcrt,digits=res)
    idx    = collect(1:1:size(matcrt,1))

    # remove duplicates dat
    sdx = sortperm(matrad[:,1])
    matrad = matrad[sdx,:]
    idx2   = vec(idx[sdx,:])

    kpdx   = findall(.!nonunique(convert(DataFrames.DataFrame,matcrt)))
    matcrt = matcrt[kpdx,:]
    matrad = matrad[kpdx,:]
    idx2   = vec(idx2[kpdx,:])

    dat1dx = findall(idx2 .<= size(matcrt,1))
    datcrt = matcrt[dat1dx,:]
    datrad = matrad[dat1dx,:]

    return datcrt, datrad
end

function prepterdat(datcrt1,datcrt2)

    matcrt = vcat(datcrt2,datcrt1)
    matcrt = round.(matcrt,digits = 1)
    idx = collect(1:1:size(matcrt,1))

    kpdx   = findall(.!nonunique(convert(DataFrames.DataFrame,matcrt)))
    matcrt = matcrt[kpdx,:]

    return matcrt
end

function getimagecentre(slp,asp)
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

function findpairs(kdtree,b,t,dim,knum)

    distances, indices = scipyspat.cKDTree.query(kdtree,b, k=knum, n_jobs=-1)
    d = distances .<= t

    dx1  = Int.((float(indices[findall(float(d) .== 1)])))

    lia = zeros(dim,1)
    lia[dx1,:] .= 1

    return lia
end

function fillmat(kdtree,datcrt,tol,kdtreedims,radius,mat2ev)

        idx = findpairs(kdtree,datcrt,tol,kdtreedims,10)
        imdx = float(reshape(idx,(radius*2,radius*2)))
        mat2ev[imdx.==1] .= 0

        return mat2ev
end

function findmincol(row)
    return findmin(row)[2]
end

function calc_horizon_lines(cellsize,peri,pcdpol,slp)

    rbins = collect(2*cellsize:sqrt(2).*cellsize:peri)
    phi_bins = collect(-pi+((pi/360)*3):pi/720:pi-((pi/360)*3))

    global mintht = ones(size(phi_bins,1)) .* 90

    for rbix = length(rbins)-1:-1:1
        global mintht
            fix1 = findall((pcdpol[:,3] .>= rbins[rbix]) .& (pcdpol[:,3] .< rbins[rbix+1]))
            if size(fix1) > (1,)
                fix2 = sortperm(pcdpol[fix1,1])
                tdx = findall(minimum(pcdpol[fix1[fix2],1]) .<= phi_bins
                                    .<= maximum(pcdpol[fix1[fix2],1]))

                temp = ones(size(phi_bins,1)) .* 90
                temp[tdx] =  pyinterp.interp1d(pcdpol[fix1[fix2],1],
                            pcdpol[fix1[fix2],2],fill_value="extrapolate")(phi_bins[tdx])

                mintht = minimum(hcat(mintht,temp),dims=2)
            end
    end

    rphi = collect(-pi:pi/1080:pi)
    rtht = pyinterp.interp1d([phi_bins[end]-2*pi; phi_bins; phi_bins[1]+2*pi],
                        vec([mintht[end];mintht;mintht[1]]),fill_value="extrapolate")(rphi)

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
