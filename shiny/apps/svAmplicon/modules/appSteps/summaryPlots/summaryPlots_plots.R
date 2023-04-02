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
    ampKeys <- getAmpliconKeys(amplicons)
    posLim <- c(min(amplicons$pos1), max(amplicons$pos2))
    sizeMax <- max(amplicons$pos2 - amplicons$pos1 + 1)
    x <- as.data.table(t(as.data.table(strsplit(gsub('[/:]', ' ', junctions$nodePair), ' '))))
    setnames(x, c("chrom1","strand1","pos1","chrom2","strand2","pos2"))    
    x[, ":="(pos1 = as.integer(pos1), pos2 = as.integer(pos2))]
    junctions <- cbind(junctions, x)
    list(
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
        positionDensity = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            junctions[
                getAmpliconKeys(junctions) == ampKeys[i], 
                getJunctionDensity(posLim[1]:posLim[2], c(pos1, pos2))
            ]
        }),
        sizeDensity = lapply(seq_along(ampKeys), function(i){ # one series for every sample+amplicon combination
            junctions[
                getAmpliconKeys(junctions) == ampKeys[i], 
                getJunctionDensity(0:sizeMax, pos2 - pos1 + 1)
            ]
        })
    )
})

#----------------------------------------------------------------------
# render a triangle plot showing paired SV ends within one or more amplicons
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

#----------------------------------------------------------------------
# plot junction density by position and SV size
#----------------------------------------------------------------------
positionDensityPlotServer <- function(junctionPlotData) {
    plot <- staticPlotBoxServer(
        "positionDensityPlot",
        maxHeight = "400px",
        lines   = TRUE,
        legend  = TRUE,
        margins = FALSE,
        create = function() {
            d <- junctionPlotData()
            par(mar = c(4.1, 4.1, 0.1, 0.1))
            plot$initializeFrame(
                xlim = d$posLim / 1e6,
                ylim = c(0, 0.05),
                xlab = "Endpoint Coordinate (Mb)",
                ylab = "Frequency"
            )
            for(i in seq_along(d$positionDensity)){
                plot$addLines(
                    x = d$positionDensity[[i]]$x / 1e6,
                    y = d$positionDensity[[i]]$y,
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
sizeDensityPlotServer <- function(junctionPlotData) {
    plot <- staticPlotBoxServer(
        "sizeDensityPlot",
        maxHeight = "400px",
        lines   = TRUE,
        legend  = TRUE,
        margins = FALSE,
        create = function() {
            d <- junctionPlotData()
            par(mar = c(4.1, 4.1, 0.1, 0.1))
            plot$initializeFrame(
                xlim = c(0, d$sizeMax),
                ylim = c(0, 0.05),
                xlab = "SV Size (bp)",
                ylab = "Density"
            )
            for(i in seq_along(d$sizeDensity)){
                plot$addLines(
                    x = d$sizeDensity[[i]]$x,
                    y = d$sizeDensity[[i]]$y,
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
