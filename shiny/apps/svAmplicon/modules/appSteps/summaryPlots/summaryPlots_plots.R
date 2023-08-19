#----------------------------------------------------------------------
# construct a composite paired quality plot, i.e. MAPQ vs. QUAL
#----------------------------------------------------------------------
getJunctionDensity <- function(allValues, values){
    density <- merge(
        data.table(x = allValues),
        data.table(x = values)[, .(y = .N), by = "x"],
        by = "x",
        all.x = TRUE
    )
    density[is.na(y), y := 0]
    density[, .(
        x = x,
        y = y / sum (y)
    )]
}
junctionPlotDataReactive <- function(selectedAmplicons, junctionTypesJunctions) reactive({
    amplicons <- selectedAmplicons()
    junctions <- junctionTypesJunctions()
    req(amplicons, junctions)
    startSpinner(session, message = "loading plot data")
    ampKeys <- getAmpliconKeys(amplicons)
    posLim <- c(min(amplicons$pos1), max(amplicons$pos2))
    sizeMax <- max(amplicons$pos2 - amplicons$pos1 + 1)
    x <- list(
        amplicons = amplicons,
        ampKeys = ampKeys,
        posLim = posLim, # in bp,
        sizeMax = sizeMax,    
        junctions = junctions,
        triangle = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            junctions[getAmpliconKeys(junctions) == ampKeys[i], {
                size <- pos2 - pos1 + 1
                .(
                    x = pos1 + size / 2, # center of SV span
                    y = size             # SV size
                )
            }]
        }),
        circles = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            x <- seq(0, pi, length.out = 25)
            junctions[getAmpliconKeys(junctions) == ampKeys[i] & 
                      edgeType == "D", {
                halfsize <- (pos2 - pos1 + 1) / 2
                center <- pos1 + halfsize
                .(
                    x = cos(x) * halfsize + center,
                    y = sin(x) * halfsize 
                )
            }, by = .(jxnKey)]
        }),
        positionDensity = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            junctions[
                getAmpliconKeys(junctions) == ampKeys[i], 
                getJunctionDensity(posLim[1]:posLim[2], c(pos1, pos2))
            ]
        }),
        sizeDensity = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            junctions[
                getAmpliconKeys(junctions) == ampKeys[i], 
                getJunctionDensity(0:sizeMax, eventSize)
            ]
        })
    )
    stopSpinner(session)
    x
})

#----------------------------------------------------------------------
# render a triangle and a circles plot showing paired SV ends within one or more amplicons
#----------------------------------------------------------------------
svTrianglePlotServer <- function(junctionPlotData) {
    plot <- staticPlotBoxServer(
        "svTrianglePlot",
        maxHeight = "400px",
        points   = TRUE,
        legend  = TRUE,
        margins = FALSE,
        create = function() {
            d <- junctionPlotData()
            par(mar = c(4.1, 4.1, 0.1, 0.1))
            plot$initializeFrame(
                xlim = d$posLim / 1e6,
                ylim = c(0, d$sizeMax),
                xlab = "SV Center (Mb)",
                ylab = "SV Size (bp)"
            )
            for(i in seq_along(d$triangle)){
                plot$addPoints(
                    x = d$triangle[[i]]$x / 1e6,
                    y = d$triangle[[i]]$y,
                    col = CONSTANTS$plotlyColors[[i]]
                )            
            }
            plot$addLegend(
                legend = d$ampKeys,
                col = unlist(CONSTANTS$plotlyColors[1:length(d$ampKeys)])
            )
        }
    )
    plot
}
svCirclesPlotServer <- function(junctionPlotData) {
    plot <- staticPlotBoxServer(
        "svCirclesPlot",
        maxHeight = "400px",
        lines   = TRUE,
        legend  = TRUE,
        margins = FALSE,
        create = function() {
            d <- junctionPlotData()
            par(mar = c(4.1, 4.1, 0.1, 0.1))
            startSpinner(session, message = "rendering circles")
            plot$initializeFrame(
                xlim = d$posLim / 1e6,
                ylim = c(0, d$sizeMax * 0.4),
                xlab = "Endpoint Coordinate (Mb)",
                ylab = "",
                yaxt = "n",
                yaxs = "i"
            )
            for(i in seq_along(d$circles)){
                dd <- d$circles[[i]]  
                col <- col2rgb(CONSTANTS$plotlyColors[[i]]) / 255                      
                plot$addLines(
                    x = dd$x / 1e6,
                    y = dd$y,
                    col = rgb(t(col), alpha = 0.8)
                )   
            }
            plot$addLegend(
                legend = d$ampKeys,
                col = unlist(CONSTANTS$plotlyColors[1:length(d$ampKeys)])
            )
            stopSpinner(session)
        }
    )
    plot
}

#----------------------------------------------------------------------
# plot junction density by position and SV size
#----------------------------------------------------------------------
densityPlotSettings <- list(Distribution = list(
    Bin_Size = list(
        type = "numericInput",
        value = 5,
        min = 1,
        max = 50,
        step = 1
    ),
    Max_Bin_Size = list(
        type = "numericInput",
        value = 0,
        min = 0,
        max = 1000,
        step = 50
    )
))
plotJunctionDensity <- function(plot, d, type, xmax, xlab, scalar = 1){
    binSize    <- plot$settings$get("Distribution","Bin_Size")
    maxBinSize <- plot$settings$get("Distribution","Max_Bin_Size")
    ymax <- 0
    for(i in seq_along(d[[type]])){
        dd <- d[[type]][[i]]
        if(is.null(binSize)) binSize <- 5
        dd[, bin := floor(x / binSize)]
        dd <- dd[, .(
            x = mean(x),
            y = sum(y)
        ), keyby = .(bin)]
        d[[type]][[i]] <- dd
        ymax <- max(ymax, dd$y)          
    }
    par(mar = c(4.1, 4.1, 0.1, 0.1))
    plot$initializeFrame(
        xlim = switch(
            type,
            sizeDensity = c(0, if(!is.null(maxBinSize) && maxBinSize > 0) maxBinSize else xmax),
            positionDensity = d$posLim
        ) / scalar,
        ylim = c(0, ymax),
        xlab = xlab,
        ylab = "Density"
    )
    for(i in seq_along(d[[type]])){
        plot$addLines(
            x = d[[type]][[i]]$x / scalar,
            y = d[[type]][[i]]$y,
            col = CONSTANTS$plotlyColors[[i]]
        )            
    }
    plot$addLegend(
        legend = d$ampKeys,
        col = unlist(CONSTANTS$plotlyColors[1:length(d$ampKeys)])
    )
}
positionDensityPlotServer <- function(junctionPlotData) {
    plot <- staticPlotBoxServer(
        "positionDensityPlot",
        maxHeight = "400px",
        lines   = TRUE,
        legend  = TRUE,
        margins = FALSE,
        settings = densityPlotSettings,
        create = function() {
            d <- junctionPlotData()
            plotJunctionDensity(plot, d, "positionDensity", d$posLim, "Endpoint Coordinate (Mb)", 1e6)
        }
    )    
    plot
}
sizeDensityPlotServer <- function(junctionPlotData) {
    plot <- staticPlotBoxServer(
        "sizeDensityPlot",
        maxHeight = "400px",
        lines   = TRUE,
        legend  = TRUE,
        margins = FALSE,
        settings = densityPlotSettings,
        create = function() {
            d <- junctionPlotData()
            plotJunctionDensity(plot, d, "sizeDensity", d$sizeMax, "SV Size (bp)")
        }
    )
    plot
}
