#----------------------------------------------------------------------
# svx_triangle is a track type for plotting SV endpoint pairs as dots on a triangle plot
# expects:
#   one or more track$settings$items from showTrackSamplesDialog()
#   dataFn(track, sourceId, samples) must return data.table(samples, center = integer64, size = integer64, cex, color, ...)
# data may arise from one or more sources, derived from Adjust Sample Sources
# each source may have one or more samples, as a result of multi-sample find actions in the relevant svX pipelines
#----------------------------------------------------------------------

# constructor
new_svx_triangleTrack <- function(trackId, expandReactive) {
    list(
        click  = TRUE,
        hover  = FALSE,
        brush  = FALSE,
        items  = TRUE,
        expand = expandReactive,
        expand2 = TRUE,
        navigation = TRUE
    )
}

# track build function
build.svx_triangle_track <- function(track, reference, coord, layout, trackBuffer, loadFn){
    req(coord, coord$chromosome)

    # get the data to plot
    selectedSources <- getSourcesFromTrackSamples(track$settings$items())

    # set the plot layout and axes
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    Max_SV_Size <- getBrowserTrackSetting(track, "Filters", "Max_SV_Size", 0)
    Color_By <- getBrowserTrackSetting(track, "Points", "Color_By", "edgeType")
    if(Max_SV_Size == 0) Max_SV_Size <- as.numeric(coord$width)
    ylim <- Max_SV_Size * c(-0.05, 1.05)

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svx_triangle", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = "SV Size", #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 
        svx_trianglechromLines(reference, coord)
        jxns <- svx_getTrackJunctions(
            track, selectedSources, loadFn, coord, 
            "centerAll", chromOnly = TRUE
        )[sample(.N)]
        jxns <- if(Color_By == "sample") dt_colorBySelectedSample(jxns, selectedSources) 
                else svx_colorCompositeGenomes(jxns, reference)
        points(
            jxns$center, 
            jxns$size,
            pch = 19,
            cex = jxns$cex,
            col = jxns$color
        )
        if(nrow(jxns) == 0) trackNoData(coord, ylim, "no matching junctions in window")
        svx_junctionsLegend(track, coord, ylim, selectedSources, 
                            hasIntergenome = "interGenome" %in% names(jxns))
        trackBuffer[[track$id]] <- jxns
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

# special diagonal chromsome demarcation lines
svx_trianglechromLines <- function(reference, coord, color = CONSTANTS$plotlyColors$grey, lwd = 0.5){
    if(!isProperChromosome(coord$chromosome)) {
        x <- getChromosomeSizes(reference$genome, reference$metadata)$size
        for(xx in c(0, cumsum(as.double(x)))){
            minX <- as.numeric(coord$start)
            maxX <- as.numeric(coord$end)            
            lines(c(xx, maxX), c(0,  2 * maxX + 2 * -xx), col = color, lwd = lwd)
            lines(c(minX, xx), c(-2 * minX + 2 *  xx, 0), col = color, lwd = lwd)
        }
    }
}
