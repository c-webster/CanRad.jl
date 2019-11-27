function DTM_create(set_in)

    eval(extract(set_in))

    if !ispath(outdir)
        mkdir(outdir)
    end

    if clip
        new[1,1:2] = xlimits[1]-buffer,ylimits[1]-buffer
        new[2,1:2] = xlimits[2]+buffer,ylimits[1]-buffer
        new[3,1:2] = xlimits[2]+buffer,ylimits[2]+buffer
        new[4,1:2] = xlimits[1]-buffer,ylimits[2]+buffer
        new[5,1:2] = new[1,1:2]


    end

end


function DTM_segment(set_in)

    eval(extract(set_in))

    # import the dtm


    # loop through the segments, clipping and saving the dtm to the analysis area dimensions
    # @rput 
    #
    # @rimport raster as raster
    #
    # R"""
    # chm <- raster()
    #
    # """







end


function DTM_resample(set_in)




end
