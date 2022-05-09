# ----------------------------------------------------------------------
# plot of outer endpoints of all molecules matching the selected junction
# ----------------------------------------------------------------------
junctionNodesPlotUI <- function(ns, width){
    staticPlotBoxUI(
        ns('junctionNodes'),
        width = width,
        title = "Supporting Molecules"
    )
}
junctionNodesPlotServer <- function(svMols){
    plot <- staticPlotBoxServer(
        'junctionNodes', 
        legend = TRUE,
        points = TRUE,
        margins = TRUE,
        immediate = TRUE,
        create = function(...){
            x <- svMols()
            req(x)
            xAllowedLim <- x$sv$POS_1 + c(-1, 1) * x$maxTLen
            yAllowedLim <- x$sv$POS_2 + c(-1, 1) * x$maxTLen
            x$mols <- x$mols[ # in case molecule has >1 jxn; only plot this junction
                OUT_POS1 %between% xAllowedLim & 
                OUT_POS2 %between% yAllowedLim
            ]
            axisSpan <- max( # ensure a sufficiently wide, and narrow, axis span to show all points well
                max(abs(x$sv$POS_1 - x$mols$OUT_POS1)), 
                max(abs(x$sv$POS_2 - x$mols$OUT_POS2)),
                x$maxTLen / 2
            ) # adjust to show the SV in its appropriate quadrant, plus outliers if present
            xmin <- if(any(x$mols$OUT_POS1 < x$sv$POS_1)) x$sv$POS_1 - axisSpan else x$sv$POS_1
            xmax <- if(any(x$mols$OUT_POS1 > x$sv$POS_1)) x$sv$POS_1 + axisSpan else x$sv$POS_1
            ymin <- if(any(x$mols$OUT_POS2 < x$sv$POS_2)) x$sv$POS_2 - axisSpan else x$sv$POS_2
            ymax <- if(any(x$mols$OUT_POS2 > x$sv$POS_2)) x$sv$POS_2 + axisSpan else x$sv$POS_2
            plot$initializeFrame(
                xlim = c(xmin, xmax) / 1e6,
                ylim = c(ymin, ymax) / 1e6,
                xlab = "Junction Coordinate 1 (Mbp)",
                ylab = "Junction Coordinate 2 (Mbp)"
            )
            abline(h = x$sv$POS_2 / 1e6) # crosshairs at the junction point
            abline(v = x$sv$POS_1 / 1e6)   
            plot$addPoints(
                x = c(x$sv$POS_1, x$mols$OUT_POS1) / 1e6, 
                y = c(x$sv$POS_2, x$mols$OUT_POS2) / 1e6,
                col = c(
                    CONSTANTS$plotlyColors$red, # red dot at the junction point
                    SVX$getMolColors(x$mols$NODE_CLASS)
                ),
                pch = c(
                    19,
                    ifelse(x$mols$IS_DUPLEX == 1, 5, ifelse(x$mols$MOL_STRAND == 0, 2, 6))
                )
            )
            plot$addLegend(
                legend = c(names(SVX$nodeClassColors), "Duplex", "Strand 1", "Strand 2"),
                col    = c(unlist(SVX$nodeClassColors), rep(CONSTANTS$plotlyColors$grey, 3)),
                pch    = c(19, 19, 19, 5, 2, 6)
            )
        }
    )
}
