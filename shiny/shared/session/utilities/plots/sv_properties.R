#----------------------------------------------------------------------
# SV junction properties plot
#----------------------------------------------------------------------
svPropertiesPlotUI <- function(ns, width){
    staticPlotBoxUI(
        ns('svProperties'),
        width = width,
        title = "SV Properties"
    )    
}
svPropertiesPlotServer <- function(settings, svPointColors, filteredSvs){
    plot <- staticPlotBoxServer(
        'svProperties', 
        legend = TRUE,
        immediate = TRUE,
        create = function(...){
            svFilters <- settings$SV_Filters()
            stepSettings <- settings$Plot_Settings()   
            svPointColors <- svPointColors()        
            svs <- filteredSvs()[, .(
                plotted = JXN_BASES != "*", # MICROHOM_LEN meaningless if not a sequenced junction
                color = svPointColors$colors,
                MICROHOM_LEN, 
                SV_SIZE,
                JXN_TYPE
            )]
            svs[
                JXN_TYPE == "T", 
                SV_SIZE := rnorm(.N, svFilters$Max_SV_Size$value, svFilters$Max_SV_Size$value / 10)
            ]
            par(mar = c(4.1, 4.1, 0.1, 1.1))
            maxX <- stepSettings$Max_Microhomology$value
            xlim <- c(-maxX, maxX)
            ylim <- log10(c(max(svFilters$Min_SV_Size$value, 1), svFilters$Max_SV_Size$value * 1.5))
            plot(
                NA, 
                NA,
                typ = "n",
                xlim = xlim,
                ylim = ylim,
                xlab = "Microhomology Length (bp)",
                ylab = "log10 SV Size (bp)"
            )
            abline(h = seq(0, 10, 1), col = "grey60")
            abline(v = seq(-100, 100, 10), col = "grey60")
            abline(v = 0)
            points(
                jitter(svs[plotted == TRUE, MICROHOM_LEN], amount = 0.5), # spread points from (i-0.5):(i+0.5)
                log10( svs[plotted == TRUE, SV_SIZE]), 
                pch = 20, 
                cex = stepSettings$Point_Size$value,
                col = svs[plotted == TRUE, color]
            )

            # add text to denote microhomology vs. insertions
            y <- ylim[1] + (ylim[2] - ylim[1]) * 0.035
            # text(xlim[1], y, "insertions", pos = 4, offset = 0.05, cex = 0.95)
            # text(xlim[2], y, expression(paste(mu, "homology")), pos = 2, offset = 0.05, cex = 0.95)
            mtext("insertion", side = 1, line = 2, at = xlim[1], adj = 0, cex = 0.95)
            mtext(expression(paste(mu, "homology")), side = 1, line = 2, at = xlim[2], adj = 1, cex = 0.95)

            # add a legend
            pointColorLegend(stepSettings, plot$settings, svPointColors)
        }
    )
}
