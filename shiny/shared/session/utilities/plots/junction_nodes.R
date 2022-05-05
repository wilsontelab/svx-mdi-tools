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
            plot$initializeFrame(
                xlim = range(c(x$sv$POS_1, x$mols$OUT_POS1)) / 1e6,
                ylim = range(c(x$sv$POS_2, x$mols$OUT_POS2)) / 1e6,
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
                )
            )
            plot$addLegend(
                legend = names(SVX$nodeClassColors),
                col   = unlist(SVX$nodeClassColors)
            )
        }
    )
}
