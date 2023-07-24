# plot histograms represent the summary of all filtered junctions

svx_summarizeJunctions <- function(jxns, track, layout, customPlots = NULL, ncol = 3){
    req(jxns)

    distributions <- c(list(
        jxns[, .N, keyby = .(nSamples)],
        jxns[, .N, keyby = .(nInstances)], # nInstances <= 50
        jxns[between(insertSize, -50, 50), .N, keyby = .(insertSize)]
    ), customPlots)
    nPlots <- length(distributions) + 2

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
    height <- getBrowserTrackSetting(track, "Summarize", "Summary_Plot_Height", 4) # height in inches
    mai <- NULL
    image <- mdiTrackImage(layout, height, message = "creating distributions", function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        layout(matrix(1:(ceiling(nPlots / ncol) * ncol), ncol = ncol, byrow = TRUE))
        par(mar = c(4.6, 4.1, 0.6, 0.1), cex = 1)
        for(distribution in distributions) plotH(distribution)
        barplot(N ~ edgeType, jxns[, .N, keyby = .(edgeType)])
        plot(jxns[, insertSize], jxns[, nInstances], pch = 19, cex = 0.2, xlab = "insertSize", ylab = "nInstances", xlim = c(-50, 50))
        plot(jxns[, insertSize], jxns[, log10(size)], pch = 19, cex = 0.2, xlab = "insertSize", ylab = "log10(size)", xlim = c(-50, 50), ylim = c(2, 7)) 
    }) 

    # return the track's magick image and associated metadata
    list(ylim  = NA, mai   = mai, image = image)
}
