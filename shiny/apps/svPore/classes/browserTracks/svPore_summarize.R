# svPore summary of all filtered junctions
svx_summarizeJunctions_app <- function(jxn, track, layout){
    req(jxn)
    jxn[, alnBaseQual := as.integer(alnBaseQual)]
    jxn[, alnSize_kb  := as.integer(alnSize / 1000)]

    plotH <- function(dt, lwd = 2){
        plot(
            dt, 
            typ = "h", 
            ylim = c(0, max(dt[[2]]) * 1.05), 
            xlab = names(dt)[1],
            ylab = "N",
            lwd = lwd
        )
    }    

    padding <- padding(track, layout)
    height <- getBrowserTrackSetting(track, "Summarize", "Plot_Height", 3) # height in inches
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "svPore_triangle summary", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        layout(matrix(1:6, nrow = 2, byrow = TRUE)) 
        par(mar = c(4.6, 4.1, 0.6, 0.1), cex = 1)

        # mapQ            = max(mapQ), # these values aggregate over all junction in the cluster, even the fuzzy matched ones
        # gapCompressedIdentity = max(gapCompressedIdentity),

        plotH(jxn[, .N, keyby = .(nSamples)])
        plotH(jxn[, .N, keyby = .(nInstances)]) # nInstances <= 50
        plotH(jxn[between(insertSize, -50, 50), .N, keyby = .(insertSize)])
        plotH(jxn[, .N, keyby = .(alnBaseQual)])
        plotH(jxn[, .N, keyby = .(alnSize_kb)])
        barplot(N ~ edgeType, jxn[, .N, keyby = .(edgeType)])
    })
    stopSpinner(session)  

    # return the track's magick image and associated metadata
    list(
        # xlim  = xlim, # will be rescaled to c(0, 1) across all stacked zoom images
        ylim  = NA,
        mai   = mai,
        image = image
    )
}
