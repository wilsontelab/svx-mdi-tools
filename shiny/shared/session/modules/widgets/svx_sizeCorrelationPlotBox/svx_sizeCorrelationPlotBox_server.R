#----------------------------------------------------------------------
# server components for plotting SV insert size vs. SV size
#----------------------------------------------------------------------

#----------------------------------------------------------------------
svx_sizeCorrelationPlotBoxServer <- function(
    id,
    jxns # junctions reactive
) { 
#----------------------------------------------------------------------
plot <- staticPlotBoxServer(
    id,
    margins = TRUE,
    points = TRUE,
    title = TRUE,
    settings = c(list(
        Size_Plot = list(
            Min_X_Value = list(
                type = "numericInput",
                value = -10,
                min = -50, 
                max = 0,
                step = 5
            ),
            Max_X_Value = list(
                type = "numericInput",
                value = 15,
                min = 0, 
                max = 50,
                step = 5
            ),
            Min_Y_Value = list(
                type = "numericInput",
                value = 2,
                min = 0, 
                max = 8,
                step = 0.1
            ),
            Max_Y_Value = list(
                type = "numericInput",
                value = 7,
                min = 0, 
                max = 8,
                step = 0.1
            ),
            Size_Plot_Type = list(
                type = "radioButtons",
                choices = c("levelplot", "scatterplot"),
                value = "levelplot"
            )
        )
    ), mdiLevelPlotSettings), 
    size = "m",
    create = function() {
        jxns <- jxns()
        xlim <- c(
            plot$settings$get("Size_Plot","Min_X_Value"),
            plot$settings$get("Size_Plot","Max_X_Value")
        )
        ylim <- c(
            plot$settings$get("Size_Plot","Min_Y_Value"),
            plot$settings$get("Size_Plot","Max_Y_Value")
        )
        title <- plot$settings$get("Plot_Frame", "Title", NULL)
        totalJxns <- paste(trimws(commify(nrow(jxns))), "Junctions")
        title <- if(is.null(title) || title == "") totalJxns else paste0(title, " (", totalJxns, ")")
        plot$initializeFrame(
            xlim = xlim,
            ylim = ylim,
            xlab = "Insert Size (bp)",
            ylab = "log10 SV Size (bp)",
            title = title,
            cex.main = 0.95
        )
        if(plot$settings$get("Size_Plot","Size_Plot_Type") == "levelplot"){
            mdiLevelPlot(
                dt = jxns[, .(x = insertSize, y = round(log10(size), 1))],
                xlim = xlim,
                xinc = 1, 
                ylim = ylim,
                yinc = 0.1,
                z.fn = length, 
                z.column = "x",
                settings = plot$settings,
                legendTitle = "# Junctions",
                v = c(seq(-50, 50, 5), -1, -2),
                h = 0:10,
                border = NA
            )            
        } else {
            d <- jxns[, .(
                x = insertSize,
                y = log10(size),
                color = svx_jxnType_codeToX(edgeType, "color")
            )][sample.int(.N)]
            plot$addPoints(
                x = jitter(d$x, a = 0.5),
                y = jitter(d$y),
                col = d$color
            )
            edgeTypes <- svx_jxnTypes$code[svx_jxnTypes$code %in% unique(jxns$edgeType)]
            plot$addMarginLegend(
                xlim[2] * 1.1, 
                ylim[2], 
                lty = NULL, 
                lwd = NULL, 
                legend = svx_jxnType_codeToX(edgeTypes, "longName"),
                col = svx_jxnType_codeToX(edgeTypes, "color"),
                bty = "n",
                pt.cex = 1,
                cex = 0.85
            )
        }
        stopSpinner(session)
    }
)

# return the plot
plot 
#----------------------------------------------------------------------
}
#----------------------------------------------------------------------
