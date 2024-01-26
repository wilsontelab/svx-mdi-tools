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
            ),
            Group_Points_By = list(
                type = "radioButtons",
                choices = c("SV Type", "Sample"),
                value = "SV Type"
            )
        )
    ), mdiLevelPlotSettings, mdiXYPlotSettings), 
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
        v <- c(seq(-50, 50, 5), -1, -2)
        h <- 0:10
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
                v = v,
                h = h,
                border = NA
            )            
        } else {
            groupBy <- plot$settings$get("Size_Plot","Group_Points_By")
            if(groupBy == "SV Type"){
                jxns[, group := svx_jxnType_codeToX(edgeType, "longName")]
                edgeTypes <- svx_jxnTypes$code[svx_jxnTypes$code %in% unique(jxns$edgeType)] 
                groupColors <- as.list(svx_jxnType_codeToX(edgeTypes, "color")) 
                names(groupColors) <- svx_jxnType_codeToX(edgeTypes, "longName")
            } else {
                jxns[, group := gsub(",", " ", samples)]
                groupColors <- NULL
            }
            mdiXYPlot(
                plot,
                dt = jxns[, .(
                    x = insertSize, 
                    y = log10(size), 
                    group = group
                )],
                groupingCols = "group",
                groupColors = groupColors,
                xlim = xlim,
                ylim = ylim,
                plotAs = "points",
                legendTitle = groupBy,
                v = v,
                h = h,
                x0Line = TRUE
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
