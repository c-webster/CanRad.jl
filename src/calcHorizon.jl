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
