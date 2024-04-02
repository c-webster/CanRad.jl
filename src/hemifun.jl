### functions for calculating the synthetic hemispheric images


function clipdat(pc_x::Vector{Float64},pc_y::Vector{Float64},pc_z::Vector{Float64},limits,peri=0::Int64)
    rmidx = (pc_x.<(limits[1]-peri)) .| (pc_x.>(limits[2]+peri)) .|
                (pc_y.<(limits[3]-peri)) .| (pc_y.>(limits[4]+peri))
    deleteat!(pc_x,rmidx)
    deleteat!(pc_y,rmidx)
    deleteat!(pc_z,rmidx)
    return pc_x, pc_y, pc_z, rmidx
end

function clipdat!(inpcx::Vector{Float64},inpcy::Vector{Float64},inpcz::Vector{Float64},
    limits::Vector{Float64},rmidx::Vector{Bool},peri=0::Number)

    rmidx .= (inpcx.<(limits[1]-peri)) .| (inpcx.>(limits[2]+peri)) .|
        (inpcy.<(limits[3]-peri)) .| (inpcy.>(limits[4]+peri))

    deleteat!(inpcx,rmidx)
    deleteat!(inpcy,rmidx)
    deleteat!(inpcz,rmidx)

end

function findelev!(inpcx::Vector{Float64},inpcy::Vector{Float64},inpcz::Vector{Float64},x,y,
    limits::Vector{Float64},peri::Number,elev::Vector{Float64},interp_method::String="linear")

    getlimits!(limits,x,y,peri)

    clipdat!(inpcx,inpcy,inpcz,limits,Vector{Bool}(undef,size(inpcx,1))),peri;
    elev .= pyinterp.griddata(hcat(inpcx,inpcy), inpcz, (x, y), method=interp_method)

end

function findelev(inpcx::Vector{Float64},inpcy::Vector{Float64},inpcz::Vector{Float64},x,y,
    peri=20.0::Number,interp_method::String="linear")

    limits = getlimits!(Vector{Float64}(undef,4),x,y,peri)

    clipdat!(inpcx,inpcy,inpcz,limits,Vector{Bool}(undef,size(inpcx,1))),peri;
    return pyinterp.griddata(hcat(inpcx,inpcy), inpcz, (x, y), method=interp_method)

end

function getlimits!(limits,pts_x,pts_y,peri)

    limits[1] = (floor(minimum(pts_x)-peri))
    limits[2] = (ceil(maximum(pts_x)+peri))
    limits[3] = (floor(minimum(pts_y)-peri))
    limits[4] = (ceil(maximum(pts_y)+peri))
    
    return limits

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

function preallo_trunks(dat_x::Vector{Float64},dat_y::Vector{Float64},dat_z::Vector{Float64},
            dat_r::Vector{Float64},npt::Any,hint::Any)
    
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

function calculate_trunks(dbh_x::Vector{Float64},dbh_y::Vector{Float64},dbh_z::Vector{Float64},
                dbh_r::Vector{Float64},npt::Any,hint::Any,e::Any)
    
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


function make_branches(ltc_x::Vector{Float64},ltc_y::Vector{Float64},ltc_z::Vector{Float64},
    ltc_tx::Vector{Float64},ltc_ty::Vector{Float64},ltc_hd::Vector{Float64},ltc_ang::Vector{Float64},
    ltc_cls::Vector{Int32},spacing=0.1::Float64,season=nothing::String)

    # calculate angular (eucledian) distane of trunk from trunk to tip
    ltc_ad = zeros(size(ltc_hd,1)) .* NaN

    ltc_ad[ltc_ang.>90] = ltc_hd[ltc_ang.>90]./cosd.(90 .- ltc_ang[ltc_ang.>90])
    ltc_ad[ltc_ang.<90]  = ltc_hd[ltc_ang.<90]./cosd.(ltc_ang[ltc_ang.<90].-90)
    ltc_ad[ltc_ang.==90] = ltc_hd[ltc_ang.==90]

    # calculate vertical 2D distance between top and bottom of branch
    ltc_zd = sqrt.(ltc_ad .^2 .- ltc_hd .^ 2)

    # calculate height at which branch intersects trunk
    ltc_tz = ltc_z + ltc_zd

    bsm_x = Vector{Float64}(); bsm_y = Vector{Float64}(); bsm_z = Vector{Float64}();
    season == "complete" ? bsm_cls = Vector{Int32}() : bsm_cls = 0
    # bsm_cls = 1 if deciduous, 0 if evergreen

    for bidx in eachindex(ltc_x)
        npts = Int32(floor(ltc_ad[bidx,1] / spacing))-1
        if npts > 2 && ltc_hd[bidx] .< 4 # limit branches to < 4m 
            # the '+spacing' is so the dsm points aren't included in the bsm dataset
            append!(bsm_x,ltc_x[bidx]+spacing .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc_tx[bidx] - ltc_x[bidx]))
            append!(bsm_y,ltc_y[bidx]+spacing .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc_ty[bidx] - ltc_y[bidx]))
            append!(bsm_z,ltc_z[bidx]+spacing .+ collect(range(0,stop=1,length=npts)) .*
                                        (ltc_tz[bidx] - ltc_z[bidx]))
            if season == "complete"
                append!(bsm_cls,Int.(repeat([ltc_cls[bidx]],npts)))
            end
        end
    end

    # make the points a bit more random
    bsm_x .+= rand(Uniform(-branch_spacing,branch_spacing),size(bsm_x,1))
    bsm_y .+= rand(Uniform(-branch_spacing,branch_spacing),size(bsm_x,1))
    bsm_z .+= rand(Uniform(-branch_spacing,branch_spacing),size(bsm_x,1))

    return bsm_x, bsm_y, bsm_z, Bool.(bsm_cls)

end

function dist(pcd_x::Vector{Float64},pcd_y::Vector{Float64},xcoor::Float64,ycoor::Float64)
    return hypot.(pcd_x.-xcoor,pcd_y.-ycoor)
end

function dist(dat_x::Float64,dat_y::Float64,xcoor::Float64,ycoor::Float64)
    return hypot.(dat_x.-xcoor,dat_y.-ycoor)
end

function dist3d(dat_x::Vector{Float64},dat_y::Vector{Float64},dat_z::Vector{Float64},xcoor::Float64,ycoor::Float64,zcoor::Float64)
    return hypot.(hypot.(dat_x.-xcoor,dat_y.-ycoor),dat_z.-zcoor)
end

function dist3d(dat_x::Float64,dat_y::Float64,dat_z::Float64,xcoor::Float64,ycoor::Float64,zcoor::Float64)
    return hypot(hypot.(dat_x.-xcoor,dat_y.-ycoor),dat_z.-zcoor)
end

function getsurfdat(dsm_x::Vector{Float64},dsm_y::Vector{Float64},dsm_z::Vector{Float64},
    bsm_x::Vector{Float64},bsm_y::Vector{Float64},bsm_z::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64)

    # for lidar points and branches:
    dsm_d = dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor).<peri
    bsm_d = dist3d(bsm_x,bsm_y,bsm_z,xcoor,ycoor,ecoor).<peri*0.5

    return append!(deleteat!(dsm_x,.!dsm_d),deleteat!(bsm_x,.!bsm_d)), 
            append!(deleteat!(dsm_y,.!dsm_d),deleteat!(bsm_y,.!bsm_d)), 
            append!(deleteat!(dsm_z,.!dsm_d),deleteat!(bsm_z,.!bsm_d))

end

function getsurfdat(dsm_x::Vector{Float64},dsm_y::Vector{Float64},dsm_z::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Number)

    # for lidar/terrain points
    dsm_d = dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor).<peri
    return deleteat!(dsm_x,.!dsm_d), deleteat!(dsm_y,.!dsm_d), deleteat!(dsm_z,.!dsm_d)
end

function getsurfdat(dsm_x::Vector{Float64},dsm_y::Vector{Float64},dsm_z::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,peri1::Int64,peri2::Int64)

    # for lidar/terrain points but with a minimum and maximum perimeter
    dsm_d = peri1 .< dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor) .< peri2
    return deleteat!(dsm_x,.!dsm_d), deleteat!(dsm_y,.!dsm_d), deleteat!(dsm_z,.!dsm_d)
end

function getsurfdat_chm(dsm_x::Vector{Float64},dsm_y::Vector{Float64},dsm_z::Vector{Float64},dsm_b::Vector{Float64},
    dsm_lavd::Vector{Float64},xcoor::Float64,ycoor::Float64,ecoor::Float64,peri::Int64)

    dsm_d = dist3d(dsm_x,dsm_y,dsm_z,xcoor,ycoor,ecoor) .< peri

    return deleteat!(dsm_x,.!dsm_d),deleteat!(dsm_y,.!dsm_d), deleteat!(dsm_z,.!dsm_d),
            deleteat!(dsm_b,.!dsm_d), deleteat!(dsm_lavd,.!dsm_d)

end

function normalise!(pcd_x::Vector{Float64},pcd_y::Vector{Float64},pcd_z::Vector{Float64},
                    xcoor::Float64,ycoor::Float64,ecoor::Float64,image_height=0.0::Float64)
    pcd_x .-= xcoor
    pcd_y .-= ycoor
    pcd_z .-= ecoor
    pcd_z .-= image_height

end

function cart2sph!(in_x::Vector{Float64},in_y::Vector{Float64},in_z::Vector{Float64},
    pcd_phi::Vector{Float64},pcd_tht::Vector{Float64})
    
    pcd_phi    .= atan.(in_y,in_x) # az
    pcd_tht    .= ((pi/2) .- (atan.(in_z,hypot.(in_x,in_y)))) * (180/pi) # elev/zenith
    pcd_tht[pcd_tht.> 90] .= 90
    in_z    .= hypot.(hypot.(in_x,in_y),in_z) # r

end

function pcd2pol2cart!(pcd_x::Vector{Float64},pcd_y::Vector{Float64},pcd_z::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,
    image_height::Float64,cellsize=0.0::Float64)     

    # for chm points
    normalise!(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,image_height)

    pcd_phi = copy(pcd_x)
    pcd_tht = copy(pcd_y)  
    cart2sph!(pcd_x,pcd_y,pcd_z,pcd_phi,pcd_tht) # pcd_z is pcd_rad in radial distance (r)

    rows = (pcd_z .< (10 .* (sqrt(2)*cellsize)))
    deleteat!(pcd_phi,rows); deleteat!(pcd_tht,rows); deleteat!(pcd_z,rows)

    return pol2cart(pcd_phi,pcd_tht)

end

function pcd2pol2cart!(pcd_x::Vector{Float64},pcd_y::Vector{Float64},pcd_z::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,dat_type::String,image_height::Float64)     

    # for lidar points
    normalise!(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,image_height)

    pcd_phi = Vector{Float64}(undef,size(pcd_x,1))
    pcd_tht = Vector{Float64}(undef,size(pcd_x,1))  
    cart2sph!(pcd_x,pcd_y,pcd_z,pcd_phi,pcd_tht) # pcd_z is pcd_rad in radial distance (r)

    pol2cart!(pcd_phi,pcd_tht,pcd_x,pcd_y)

    # if dat_type=="trunks"
    # return pcd_x, pcd_y
    # else
    return pcd_x, pcd_y, pcd_z
    # end

end

function pcd2pol2cart!(ter2rad::TER2RAD,pcd_x::Vector{Float64},pcd_y::Vector{Float64},pcd_z::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,dat_type::String,rbins::Vector{Float64},
    image_height::Float64,slp=0.0::Float64)     

    # for terrain points
    normalise!(pcd_x,pcd_y,pcd_z,xcoor,ycoor,ecoor,image_height)

    pcd_phi = Vector{Float64}(undef,size(pcd_x,1))
    pcd_tht = Vector{Float64}(undef,size(pcd_x,1))  
    cart2sph!(pcd_x,pcd_y,pcd_z,pcd_phi,pcd_tht) # pcd_z is pcd_rad in radial distance (r)

    pcd_phi, pcd_tht = calc_horizon_lines(ter2rad,pcd_phi,pcd_tht,pcd_z,rbins,slp,dat_type)

    return pol2cart(pcd_phi,pcd_tht)

end

function pol2cart!(pol_phi::Vector{Float64},pol_tht::Vector{Float64},pcd_x::Vector{Float64},pcd_y::Vector{Float64})
    pcd_x .= pol_tht .* cos.(pol_phi)
    pcd_y .= pol_tht .* sin.(pol_phi)
end

function pol2cart(pcd_phi::Vector{Float64},pcd_tht::Vector{Float64})
    return pcd_tht .* cos.(pcd_phi), pcd_tht .* sin.(pcd_phi)
end

function hlm2cart(ter2rad::TER2RAD,rtht::Vector{Float64})

    @unpack rphi = ter2rad

    pol_phi, pol_tht = fillterrain(rphi,rtht,0.0)
    pt_dtm_x, pt_dtm_y = pol2cart(pol_phi,pol_tht)
    return prepterdat(pt_dtm_x,pt_dtm_y)

end

function filterbyradius(phi::Vector{Float64},tht::Vector{Float64},rad::Vector{Float64},peri::Int64)
    rmdx = rad .> peri
    return deleteat!(phi,rmdx), deleteat!(tht,rmdx), deleteat!(rad,rmdx)
end


function prepterdat!(matcrt_x::Vector{Float64},matcrt_y::Vector{Float64})

    matcrt_x .= round.(matcrt_x,digits = 1)
    matcrt_y .= round.(matcrt_y,digits = 1)

    rmdx = nonunique(DataFrame(hcat(matcrt_x,matcrt_y),:auto));

    deleteat!(matcrt_x,rmdx)
    deleteat!(matcrt_y,rmdx)

end

function prepterdat(matcrt_x::Vector{Float64},matcrt_y::Vector{Float64})

    matcrt_x .= round.(matcrt_x,digits = 1)
    matcrt_y .= round.(matcrt_y,digits = 1)

    rmdx = nonunique(DataFrame(hcat(matcrt_x,matcrt_y),:auto));

    return deleteat!(matcrt_x,rmdx), deleteat!(matcrt_y,rmdx)

end

function prepsurfdat!(matcrt_x::Vector{Float64},matcrt_y::Vector{Float64},matcrt_z::Vector{Float64})

    matcrt_x .= round.(matcrt_x,digits = 1)
    matcrt_y .= round.(matcrt_y,digits = 1)
    matcrt_z .= round.(matcrt_z,digits = 1)

    rmdx = nonunique(DataFrame(hcat(matcrt_x,matcrt_y),:auto));

    deleteat!(matcrt_x,rmdx)
    deleteat!(matcrt_y,rmdx)
    deleteat!(matcrt_z,rmdx)

end

function findpairs(kdtree::Any,datcrt::Matrix{Float64},knum::Number,lia::BitVector)

    lia[scipyspat.cKDTree.query(kdtree,datcrt, k=knum)[2],:] .= 0
    return lia

end

function fillmat!(canrad::CANRAD,kdtree::PyObject,datcrt::Matrix{Float64},
    knum::Number,mat2ev::Matrix{Int64})

    @unpack diameter, lia = canrad
    fill!(lia,1)
    mat2ev .*= (reshape(findpairs(kdtree,datcrt,knum,lia),(diameter,diameter)))

end

function findmincol(row)
    return findmin(row)[2]
end

function frbins(pcd::Vector{Float64},rbin1::Float64,rbin2::Float64)
    return ((pcd .>= rbin1) .& (pcd .< rbin2))
end

function remove_duplicates!(indat::Vector{Float64})
    # remove indices in fix1 that create duplicates
    # method taken from dedeplicate_knot!() in Interpolations.jl [line 121:gridded.jl]
    # used here to prevent warning in Interpolations caused phi values along 
    # same azimuth included in zenith ring (happens along N,E,S,W trajectories)

    last_val = first(indat)
    for i in 2:length(indat)
        if indat[i] == last_val
            indat[i] = nextfloat(indat[i-1])
        else
            last_val = indat[i]
        end
    end

end

function calcmintht!(ter2rad::TER2RAD,mintht::Vector{Float64},pcd_phi::Vector{Float64},pcd_tht::Vector{Float64},pcd_rad::Vector{Float64},
    rbins::Vector{Float64})

    @unpack phi_bins, fix1, tdx, tempmintht = ter2rad

    fill!(fix1,0)
    fill!(tdx,0)
    fill!(tempmintht,90)

    idx = 1:size(pcd_phi,1)

    @inbounds @simd for rbix = length(rbins)-1:-1:1
        
            fix1 = idx[frbins(pcd_rad,rbins[rbix],rbins[rbix+1])]

            if size(fix1) > (1,)

                phi_srtdx = sortperm(pcd_phi[fix1])
                phi_vals = pcd_phi[fix1][phi_srtdx]
                !allunique(phi_vals) && remove_duplicates!(phi_vals)

                tdx = minimum(phi_vals) .<= phi_bins .<= maximum(phi_vals)
                fill!(tempmintht,90)
                tempmintht[tdx] = linear_interpolation(phi_vals,pcd_tht[fix1[phi_srtdx]])(phi_bins[tdx])
                mintht .= vec(minimum(hcat(mintht,tempmintht),dims=2))

            end

    end

end

function fillterrain(rphi::Vector{Float64},rtht::Vector{Float64},slp=0.0::Float64)

    min_rphi = minimum(filter(!isnan,rtht))
    int = size(collect(min_rphi:0.5:(90+slp)),1)
    int==1 && (int=2) # below steps don't work if interval is 1
    temp1 = repeat(collect(range(0,stop=1,length=int))',outer=[size(rtht,1),1])
    temp2 = repeat((ones(size(rtht,1),1) .* (90+slp)) - rtht,outer=[1,int])

    # tht = vec(repeat(rtht,outer=[1,int]) + temp1 .* temp2)
    # phi = vec(repeat(rphi,outer=[1,int]))

    return vec(repeat(rphi,outer=[1,int])), vec(repeat(rtht,outer=[1,int]) + temp1 .* temp2) # azimuth, zenith

end

function calc_horizon_lines(ter2rad::TER2RAD,pcd_phi::Vector{Float64},pcd_tht::Vector{Float64},pcd_rad::Vector{Float64},
    rbins::Vector{Float64},slp::Float64,dat_type=nothing::String)

    pcd_phitemp = copy(pcd_phi)

    # create two loops around horizon line to avoid artefacts at -pi/pi
    pcd_phitemp[pcd_phi .<= 0]  .+= 2*pi
    pcd_phitemp[pcd_phi .> 0] .-= 2*pi

    append!(pcd_phi,pcd_phitemp)
    pcd_tht = repeat(pcd_tht,outer=2)
    pcd_rad = repeat(pcd_rad,outer=2)

    @unpack mintht = ter2rad
    fill!(mintht,90)
    calcmintht!(ter2rad,mintht,pcd_phi,pcd_tht,pcd_rad,rbins);

    if dat_type=="buildings"
        # calculate moving average across building tops to smooth surface
        temp = vcat(fill(90.0,4),mintht,fill(90.0,4))
        dx = mintht .> 75
        temp[temp .> 75] .= NaN
        temp = sma(temp,9)
        temp[dx] .= mintht[dx]
        mintht = temp
        # rphi = collect(-pi:pi/720:pi)
    end

    if sum(isnan.(mintht)) .> 0
        mintht .= Impute.interp(replace(mintht,NaN=>missing))
    end

    @unpack rphi, rtht, dx1, dx2, phi_bins = ter2rad
    # increase sampling along horizonline to create opaque terrain
    rtht = linear_interpolation(phi_bins[dx1:dx2],vec(mintht[dx1:dx2]))(rphi)

    return fillterrain(rphi,rtht,slp)

end

function normalise_chmdat!(pcd_x::Vector{Float64},pcd_y::Vector{Float64},pcd_z::Vector{Float64},pcd_b::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,image_height=0.0::Float64)
    
    pcd_x .-= xcoor
    pcd_y .-= ycoor
    pcd_z .-= ecoor
    pcd_z .-= image_height
    pcd_b .-= ecoor
    pcd_b .-= image_height

end

function calcCHM_Ptrans!(chm2rad::CHM2RAD,pcd_x::Vector{Float64},pcd_y::Vector{Float64},pcd_z::Vector{Float64},
    pcd_b::Vector{Float64},lavd::Vector{Float64},
    xcoor::Float64,ycoor::Float64,ecoor::Float64,
    image_height::Float64,cellsize::Float64,rbins::Vector{Float64},cbh::Float64)

    normalise_chmdat!(pcd_x,pcd_y,pcd_z,pcd_b,xcoor,ycoor,ecoor,image_height)

    chm_phi = Vector{Float64}(undef,size(pcd_x,1))
    chm_tht = Vector{Float64}(undef,size(pcd_x,1))
    cart2sph!(pcd_x,pcd_y,pcd_z,chm_phi,chm_tht)

    bse_phi = Vector{Float64}(undef,size(pcd_x,1))
    bse_tht = Vector{Float64}(undef,size(pcd_x,1))
    cart2sph!(pcd_x,pcd_y,pcd_b,bse_phi,bse_tht)

    # note pcd_z and pcd_b are no longer height but radial distance from camera)

    pcd_phitemp = copy(chm_phi)
    pcd_phitemp[chm_phi .<= 0] .+= 2*pi
    pcd_phitemp[chm_phi .> 0]  .-= 2*pi

    append!(chm_phi,pcd_phitemp)
    chm_tht = repeat(chm_tht,outer=2)
    pcd_z   = repeat(pcd_z,outer=2) # radial distance
    lavd    = repeat(lavd,outer=2)

    @unpack sum_lavd_thick, sum_thick = chm2rad
    fill!(sum_lavd_thick,0.0)
    fill!(sum_thick,0.0)
    calcThickness!(chm2rad,rbins,sum_lavd_thick, sum_thick,
        chm_phi,chm_tht,pcd_z,lavd,bse_phi,bse_tht,pcd_b,cellsize,cbh)

    pt_chm_x, pt_chm_y             = calcPtrans(chm2rad,sum_lavd_thick,cellsize)
    pt_chm_x_thick, pt_chm_y_thick = calcPtrans_dist(chm2rad,sum_thick)

    # return calcPtrans(canrad,sum_lavd_thick,cellsize), calcPtrans_dist(sum_thick)
    return pt_chm_x, pt_chm_y, pt_chm_x_thick, pt_chm_y_thick

end

function calcThickness!(chm2rad::CHM2RAD,rbins::Vector{Float64},sum_lavd_thick::Matrix{Float64},sum_thick::Matrix{Float64},
                        chm_phi::Vector{Float64},chm_tht::Vector{Float64},chm_rad::Vector{Float64},lavd::Vector{Float64},
                        bse_phi::Vector{Float64},bse_tht::Vector{Float64},bse_rad::Vector{Float64},
                        cellsize::Float64,cbh::Float64)

    @unpack phi_bins_long, phi_bins_short, dx1_pbl, dx2_pbl = chm2rad
    @unpack chm_temptht_long, chm_temptht_short, bse_temptht, tdx_chm, tdx_bse = chm2rad
    @unpack temp_lavd, temp_thick  = chm2rad

    idx_chm = 1:size(chm_phi,1)
    idx_bse = 1:size(bse_phi,1)

    for rbix = length(rbins)-1:-1:1

        fix_chm = idx_chm[frbins(chm_rad,rbins[rbix],rbins[rbix+1])]
        if size(fix_chm) > (1,)

            phi_srtdx = sortperm(chm_phi[fix_chm])
            phi_vals  = chm_phi[fix_chm][phi_srtdx]

            !allunique(phi_vals) && remove_duplicates!(phi_vals)

            tdx_chm .= (minimum(phi_vals) .<= phi_bins_long .<= maximum(phi_vals))

            fill!(chm_temptht_long,90.0)
            chm_temptht_long[tdx_chm] = linear_interpolation(phi_vals,chm_tht[fix_chm[phi_srtdx]])(phi_bins_long[tdx_chm])

            (sum(isnan.(chm_temptht_long)) > 0) && (chm_temptht_long[isnan.(chm_temptht_long)] .= 90.0)

            fill!(temp_lavd,0.0) 
            temp_lavd[tdx_chm] = linear_interpolation(phi_vals,lavd[fix_chm[phi_srtdx]])(phi_bins_long[tdx_chm])

            fill!(chm_temptht_short,90)
            chm_temptht_short .= Int.(round.(chm_temptht_long[dx1_pbl:dx2_pbl]))
            chm_temptht_short[chm_temptht_short .== 0] .= 1 # minimum elevation angle is 1

            if cbh .> 0
                # get the canopy base line
                fix_bse = idx_bse[frbins(bse_rad,rbins[rbix],rbins[rbix+1])]
                if size(fix_bse) > (1,)
                    phi_srtdx_bse = sortperm(bse_phi[fix_bse])
                    phi_vals_bse  = bse_phi[fix_bse][phi_srtdx_bse]
                    !allunique(phi_vals_bse) && remove_duplicates!(phi_vals_bse)
                    tdx_bse .= (minimum(phi_vals_bse) .<= phi_bins_short .<= maximum(phi_vals_bse))
                    fill!(bse_temptht,90.0)
                    bse_temptht[tdx_bse] = linear_interpolation(phi_vals_bse,
                                bse_tht[fix_bse[phi_srtdx_bse]])(phi_bins_short[tdx_bse])
                    sum(isnan.(bse_temptht)) > 0 && (bse_temptht[isnan.(bse_temptht)] .= 90.0)
                end
            end

            # add the thickness to each cell
            # would be great to work out how to do this without the loop
            fill!(temp_thick,0.0)
            for tx = 1:size(chm_temptht_short,1)
                if (bse_temptht[tx])-(chm_temptht_short[tx]) > 0 # zenith angle
                    temp_thick[tx,round(Int,91-bse_temptht[tx]):(91-chm_temptht_short[tx])] .= sqrt(2)*cellsize # elevation angle
                else
                    temp_thick[tx,1:(90-Int.(chm_temptht_short[tx]))] .= sqrt(2)*cellsize
                end
            end

            sum_lavd_thick .+= (temp_thick .* temp_lavd[dx1_pbl:dx2_pbl])
            sum_thick .+= temp_thick

        end

    end

end

function calcPtrans(chm2rad::CHM2RAD,sum_lavd_thick::Matrix{Float64},cellsize::Float64)

    @unpack Ptrans = chm2rad

    Ptrans .= exp.(-0.5 * sum_lavd_thick)

    # solve pt based on average pt and random distribution
    rand_Ptrans = rand(Uniform(0,1),size(Ptrans))

    Ptrans[Ptrans .>= rand_Ptrans] .= 1
    Ptrans[Ptrans .< rand_Ptrans] .= 0

    # smooths horizon line in dense forests by removing all points where canopy thickness is only one Rbin
    # Ptrans[sum_thick .<= ((sqrt(2)*cellsize))] .= 1

    return getPhiTht(chm2rad,Ptrans)

end


function calcPtrans_dist(chm2rad::CHM2RAD,sum_thick::Matrix{Float64})

    @unpack phi_bins_short, canopy_thick  = chm2rad

    fill!(canopy_thick,0)
    phi_dist = repeat(phi_bins_short,90);
    tht_dist = repeat(90.0:-1:1,inner=size(phi_bins_short,1))

    # zero transmissivity where the canopy is > X m thick.
    canopy_thick[sum_thick .> 30] .= 0
    canopy_thick[sum_thick .<= 30] .= 1

    rows = (vec(canopy_thick).==1)
    deleteat!(phi_dist,rows)
    deleteat!(tht_dist,rows)

    return pol2cart(phi_dist,float.(tht_dist))

end

function getPhiTht(chm2rad::CHM2RAD,Ptrans::Matrix{Float64}) 
    # randomises the points within the phi/tht bins so looks more like a point cloud

    @unpack phi_bins_short = chm2rad
    
    phi = repeat(phi_bins_short,90) 
    phi .+= rand(Uniform(-pi/360,pi/360),size(phi))

    tht = repeat(90.0:-1:1,inner=size(phi_bins_short,1))
    tht .+= rand(Uniform(0,1),size(tht))
    
    # remove regions of sky where transmissivity = 1
    rows = (vec(Ptrans).==1)
    deleteat!(phi,rows)
    deleteat!(tht,rows)

    # pt_chm_x, pt_chm_y = pol2cart(phi,float.(tht))
    return pol2cart(phi,float.(tht))

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
