function calc_horizon_canopy(pcd,xcoor,ycoor,ecoor,peri1,peri2,cellsize)

    pcdpol = Array{Float64,1}(undef,3, size(pcd,2))
    pcdpol = pcd2pol(pcd,xcoor,ycoor,ecoor,peri2)

    rbins = Array{Float64,1}(undef, Int(round((peri2-peri1)/(sqrt(2).*cellsize))))
    rbins = collect(peri1:sqrt(2).*cellsize:peri2)

    tbins = Array{Float64,1}(undef, 1429)
    tbins = collect(-pi+((pi/360)*3):pi/720:pi-((pi/360)*3))

    global minphi = ones(size(tbins,1)) .* 90
    global minrad = ones(size(tbins,1)) .* peri2
    global thick  = zeros(size(tempphi,1),90)

    for rbix = length(rbins)-1:-1:1
        global minphi, minrad
            fix1 = findall((pcdpol[:,3] .>= rbins[rbix]) .& (pcdpol[:,3] .< rbins[rbix+1]))
            if size(fix1) > (1,)
                fix2 = sortperm(pcdpol[fix1,1])
                tdx = findall(minimum(pcdpol[fix1[fix2],1]) .<= tbins
                                    .<= maximum(pcdpol[fix1[fix2],1]))

                tempphi = ones(size(tbins,1)) .* 90
                tempphi[tdx] =  pyinterp.interp1d(pcdpol[fix1[fix2],1],
                            pcdpol[fix1[fix2],2])(tbins[tdx])

                flag = vec(minphi .> tempphi)
                minrad[flag.==1,:] .= mean([rbins[rbix],rbins[rbix+1]])
                minphi = minimum(hcat(minphi,tempphi),dims=2)

                tempthick = zeros(size(tempphi,1),90)
                tempphi = Int.(round.(tempphi))
                for tx = 1:1:size(tempphi,1)
                    if tempphi[tx] == 90
                        continue
                    else
                        tempthick[tx,1:tempphi[tx]] .= cellsize
                    end
                end

                thick = thick .+ tempthick

            end
    end

    rtht = collect(-pi:pi/360:pi)
    rphi = pyinterp.interp1d([tbins[end]-2*pi; tbins; tbins[1]+2*pi],
                        vec([minphi[end];minphi;minphi[1]]),fill_value="extrapolate")(rtht)
    rrad = pyinterp.interp1d([tbins[end]-2*pi; tbins; tbins[1]+2*pi],
                        vec([minrad[end];minrad;minrad[1]]),fill_value="extrapolate")(rtht)

    return rphi,rrad,thick
end
