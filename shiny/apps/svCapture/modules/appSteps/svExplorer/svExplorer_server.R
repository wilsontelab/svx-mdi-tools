#----------------------------------------------------------------------
# reactive components to filter and examine SV locations and junction sequences
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
svExplorerServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'svExplorer' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
    templates = list(file.path(app$sources$suiteGlobalDir, "sv_filters.yml"), id),
    fade = FALSE
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues() # logical failure vectors keyed as [[sampleSet]]

#----------------------------------------------------------------------
# generate the list of all filtered SVs from all selected samples
#----------------------------------------------------------------------
targetClasses <- reactive({ SVX$targetClasses[[settings$SV_Filters()$Target_Class$value]] })
workingSvs <- reactive({
    getWorkingSvs(settings, sampleSelector)
})
selectedSv <- reactive({
    rowI <- svsTable$rows_selected()
    req(rowI)
    svs <- workingSvs()
    svs[rowI]
})
workingMols <- reactive({
    getWorkingMols(selectedSv(), sampleSelector)
})

#----------------------------------------------------------------------
# set SV plot point colors
#----------------------------------------------------------------------
svPointColors <- reactive({
    svs <- workingSvs()
    req(svs)
    filters <- settings$SV_Filters()
    ps <- settings$Plot_Settings()
    setkey(SVX$jxnTypes, code)
    switch(ps$Point_Color$value,
        type = list(
            colors = svs[, SVX$jxnTypes[JXN_TYPE, color]],
            color = SVX$jxnTypes[!is.na(color) & name %in% filters$SV_Type$value, color],
            label = SVX$jxnTypes[!is.na(color) & name %in% filters$SV_Type$value, name]
        ),
        sample = {
            samples <- svs[N_SAMPLES == 1, sort(unique(SAMPLES))]
            sampleColors <- seq_along(samples)
            names(sampleColors) <- samples
            list(
                colors = svs[, ifelse(N_SAMPLES == 1, sampleColors[SAMPLES], "grey")],
                color = sampleColors,
                label = samples
            )
        },
        target = {
            targets <- svs[!grepl(',', TARGET_REGION), sort(unique(TARGET_REGION))]
            targetColors <- seq_along(targets)
            names(targetColors) <- targets
            list(
                colors = svs[, ifelse(grepl(',', TARGET_REGION), "grey", targetColors[TARGET_REGION])],
                color = targetColors,
                label = targets
            )
        }, 
        duplex = list(
            colors = svs[, ifelse(N_DUPLEX_GS > 0, "green3",  "red3")],
            color = c("green3",  "red3"),
            label = c("Duplex", "Single Strand")
        ),
        blue = list(
            colors = "blue", 
            color = NULL,
            label = NULL
        )
    )
})
pointColorLegend <- function(stepSettings, plotSettings, svPointColors, exclude = character()){
    if(is.null(svPointColors$label)) return()
    is <- which(!(svPointColors$label %in% exclude))
    legend(
        plotSettings$get('Plot_Frame', 'Legend_Placement'),
        svPointColors$label[is],
        pch = 20, 
        pt.cex = stepSettings$Point_Size$value * 1.25,
        col = svPointColors$color[is]
    )
}

#----------------------------------------------------------------------
# SV locations plot
#----------------------------------------------------------------------
locationsPlot <- staticPlotBoxServer(
    'svLocations', 
    legend = TRUE,
    immediate = TRUE,
    create = function(...){

        # load targets bed, assumed to be the same for all sample sources
        assignments <- sampleSelector$selectedAssignments() # Source_ID	Project	Sample_ID	Category1	Category2	uniqueId
        req(assignments)
        req(nrow(assignments) > 0)
        sourceId <- assignments[, Source_ID[1]]
        targetsBed <- loadPersistentFile(sourceId = sourceId, contentFileType = "targetsBed") 
        targets <- persistentCache[[targetsBed]]$data

        # initialize the data
        filters <- settings$SV_Filters()
        stepSettings <- settings$Plot_Settings()
        svs <- workingSvs()[, .(TARGET_POS_1, TARGET_POS_2)]
        svs[, ':='(
            size = abs(TARGET_POS_2 - TARGET_POS_1 + 1),
            center = pmin(TARGET_POS_1, TARGET_POS_2) + abs(TARGET_POS_2 - TARGET_POS_1 + 1) / 2
        )]
        svPointColors <- svPointColors()

        # initialize the plot
        par(mar = c(4.1, 4.1, 0.1, 0.1))
        xlim <- c(min(targets$paddedStartI), max(targets$paddedEndI))
        ylim <- c(0, if(any(c("tt", "ta", "aa") %in% targetClasses())){
            targets[, max(paddedEndI)]
        } else {
            targets[, .(x = endI - paddedStartI + 1), by = regionName][, max(x)]
        })
        
        plot(
            NA, NA, typ = "n",
            xlim = xlim,
            ylim = ylim,
            xlab = "SV Center (Mbp)",
            ylab = "SV Size (bp)",
            xaxs = "i", 
            yaxs = "i",
            xaxt = "n"
        )

        # shade and demarcate the capture target regions
        targets[, {
            xinc <- ylim[2] / 2
            polygon(
                x = c(startI, startI + xinc, endI + xinc, endI, startI), 
                y = c(ylim[1], ylim[2], ylim[2], ylim[1], ylim[1]), 
                border = NA, col = "grey90"
            )
            polygon(
                x = c(startI, startI - xinc, endI - xinc, endI, startI), 
                y = c(ylim[1], ylim[2], ylim[2], ylim[1], ylim[1]), 
                border = NA, col = "grey90"
            )
            mtext(paste(regionName, chrom, sep = ","), side = 1, line = 2.25, at = centerI, cex = 1)
            paddedSize <- paddedEndI - paddedStartI + 1
            unit <- paddedSize / 10
            at <- seq(unit, paddedSize - unit, unit)
            axis(1, at = paddedStartI + at, labels = round((start - (startI - paddedStartI) + at) / 1e6, 2))
        }, by = regionName]    
        for(i in unlist(targets[, .(paddedStartI, startI, endI, paddedEndI)])){
            abline(-i * 2,  2, col = "grey60")
            abline( i * 2, -2, col = "grey60")
        }

        # plot the SV points on top
        points(
            svs[, center], 
            svs[, size], 
            pch = 20, 
            cex = stepSettings$Point_Size$value,
            col = svPointColors$colors
        )

        # add a legend
        pointColorLegend(stepSettings, locationsPlot$settings, svPointColors)
    }
)

#----------------------------------------------------------------------
# SV junction properties plot
#----------------------------------------------------------------------
propertiesPlot <- staticPlotBoxServer(
    'svProperties', 
    legend = TRUE,
    immediate = TRUE,
    create = function(...){
        filters <- settings$SV_Filters()
        stepSettings <- settings$Plot_Settings()   
        svPointColors <- svPointColors()        
        svs <- workingSvs()[, .(
            plotted = JXN_BASES != "*", # MICROHOM_LEN meaningless if not a sequenced junction
            color = svPointColors$colors,
            MICROHOM_LEN, 
            SV_SIZE,
            JXN_TYPE
        )]
        svs[JXN_TYPE == "T", SV_SIZE := rnorm(.N, filters$Max_SV_Size$value, filters$Max_SV_Size$value / 10)]
        par(mar = c(4.1, 4.1, 0.1, 1.1))
        xlim <- c(-40, 40)
        # xlim <- c(-100, 100)
        ylim <- log10(c(max(filters$Min_SV_Size$value, 1), filters$Max_SV_Size$value * 1.5))
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
            jitter(svs[plotted == TRUE, MICROHOM_LEN]), 
            log10( svs[plotted == TRUE, SV_SIZE]), 
            pch = 20, 
            cex = stepSettings$Point_Size$value,
            col = svs[plotted == TRUE, color]
        )

        # add text to denote microhomology vs. insertions
        y <- ylim[1] + (ylim[2] - ylim[1]) * 0.035
        text(xlim[1], y, "insertions", pos = 4, offset = 0.05, cex = 0.9)
        text(xlim[2], y, "microhomologies", pos = 2, offset = 0.05, cex = 0.9)
        
        # add a legend
        pointColorLegend(stepSettings, propertiesPlot$settings, svPointColors)
    }
)

# ----------------------------------------------------------------------
# summary table of all filtered SVs from all selected samples
# ----------------------------------------------------------------------
invalidateTable <- reactiveVal(0)
svsTable <- bufferedTableServer(
    id = 'svsTable',
    parentId = id,
    parentInput = input,
    selection = 'single',
    tableData = reactive({
        invalidateTable()
        setkey(SVX$jxnTypes, code)
        workingSvs()[, .(
            svId = SV_ID,
            #---------------
            type = SVX$jxnTypes[JXN_TYPE, name],
            class = TARGET_CLASS,
            target = TARGET_REGION,
            #---------------
            nSmp = N_SAMPLES,
            nTot = N_TOTAL,
            nSplit = N_SPLITS,            
            nGap = N_GAPS,
            nClip = N_OUTER_CLIPS,
            nDpx = N_DUPLEX_GS,
            #---------------
            count = STRAND_COUNT_GS,
            count1 = STRAND_COUNT1,
            count2 = STRAND_COUNT2,
            #---------------
            chr1 = CHROM_1,
            pos1 = POS_1,
            chr2 = CHROM_2,
            pos2 = POS_2,
            #---------------
            uHom = MICROHOM_LEN,
            jxnSeq = JXN_BASES,
            svSize = SV_SIZE
        )]
    })
    # 'SHARED_PROPER' =  'double',
    # 'SHARED_PROPER_GS'  =  'double',
    # #---------------
    # 'N_AMBIGUOUS'   =  'integer',
    # 'N_DOWNSAMPLED' =  'integer',
    # 'N_COLLAPSED'   =  'integer',
    #     lm[, ':='(
    #         FAILED = tableCheckboxes(ns('libraryFailed'), failed ),
    #         Failed = failed,
    #         nReadPairs = commify(nReadPairs),
    #         nSourceMolecules = commify(nSourceMolecules),
    #         onTargetCoverage = commify(round(onTargetCoverage, 0)),
    #         offTargetCoverage = round(offTargetCoverage, 3),
    #         enrichment = commify(round(enrichment, 0)),
    #         efficiency = round(efficiency, 3)
    #     )]
    # editBoxes = list( # handle the libraryFailed checkbox
    #     libraryFailed = list(
    #         type = 'checkbox',
    #         boxColumn = 1,
    #         rawColumn = 2,
    #         handler = function(checked){ # enter the new failure value into our outcomes
    #             ss <- sampleSet$input$sampleSet
    #             failed <- outcomes[[ss]]
    #             failed[checked$selectedRow] <- checked$newValue
    #             outcomes[[ss]] <- failed # must replace the entire value for outcomes invalidation
    #             checked
    #         }
    #     )
    # )
)

# ----------------------------------------------------------------------
# plot of node points of all molecules matching a specific selected junction
# ----------------------------------------------------------------------


# get rightmost mapped read position in reference genome from POS and CIGAR
getEnd <- Vectorize(function(start, cigar) {
    lengths    <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    operations <-            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    allowed <- operations %notin% c('S', 'I')
    start - 1 + sum(lengths[allowed])
})


junctionNodesPlot <- staticPlotBoxServer(
    'junctionNodes', 
    # legend = TRUE,
    immediate = TRUE,
    create = function(...){
        mols <- workingMols()
        req(mols)
        sv <- selectedSv()

        mols <- mols[NODE_CLASS != SVX$nodeClasses$OUTER_CLIP]

        junctionNodesPlot$initializeFrame(
            xlim = sv$POS_1 + c(-1, 1) * 500,
            ylim = sv$POS_2 + c(-1, 1) * 500
        )
        abline(h = sv$POS_2)
        abline(v = sv$POS_1)        
        junctionNodesPlot$addPoints(
            x = c(sv$POS_1, mols$POS_1), 
            y = c(sv$POS_2, mols$POS_2),
            # x = c(sv$POS_1, if(sv$SIDE_1 == "L") mols$POS_1 else getEnd(mols$POS_1, mols$CIGAR_1)), 
            # y = c(sv$POS_2, if(sv$SIDE_2 == "L") mols$POS_2 else getEnd(mols$POS_2, mols$CIGAR_2)), 
            col = c(
                CONSTANTS$plotlyColors$red, 
                ifelse(mols$NODE_CLASS == SVX$nodeClasses$SPLIT, 
                       CONSTANTS$plotlyColors$blue, 
                       CONSTANTS$plotlyColors$orange)
            )
        )
    }
)

# ----------------------------------------------------------------------
# text-based zoom of a single selected junction
# ----------------------------------------------------------------------
zoomInfo <- reactive({
    reportProgress('getZoomInfo')
    rowI <- svsTable$rows_selected()
    req(rowI)
    svs <- workingSvs()
    svs[rowI]
    # smp  <- svTable[rowI,smp]
    # svId <- svTable[rowI,svId]
    # nodes <- getSvMoleculesFromList(smp, svId)
    # nodes[,c('chrom','side','pos') := unpackNodeNames(NODE)]
    # list(
    #     smp  = smp,
    #     svId = svId,
    #     sv = sampleData[[smp]]$svTable[SV_ID == svId],
    #     nodes = nodes[order(side,chrom,pos)] # thus, nodes are in canonical order
    # )
})
isSequencedJunction <- function(sv) sv[, !is.na(JXN_SEQ)] # either sequenced or reconstructed
isReconstructedJunction <- function(sv) sv[, !is.na(JXN_SEQ) & !is.na(MERGE_LEN)]
output$junctionZoom <- renderText({
    reportProgress('makeJunctionZoom')
    zoom <- zoomInfo()
    req(zoom)
    req(zoom$JXN_BASES != "*")

    return(zoom$JXN_SEQ)

    # # determine the requested alignment type
    # if(input$jxnDisplayType == "Genome"){
    #     offset <- zoom$idx$ALIGNS_OFFSET
    #     length <- zoom$idx$ALIGNS_LENGTH
    # } else if(input$jxnDisplayType == "Genome+") {
    #     offset <- zoom$idx$ALIGNS_PLUS_OFFSET
    #     length <- zoom$idx$ALIGNS_PLUS_LENGTH
    # } else if(input$jxnDisplayType == "Reads") {
    #     offset <- zoom$idx$READS_OFFSET
    #     length <- zoom$idx$READS_LENGTH        
    # }
    
    # get the alignment data
    data <- getSvIndexed(zoom$smp, zoom$svId, 'alignments', offset, length)
    
    # parse into readable rows
    nLines <- length(data)
    nChar <- nchar(data[1])
    charPerLine <- 100 # as.numeric(input$charPerLine)
    nChunks <- ceiling(nChar / charPerLine)
    nCharLastLine <- nChar %% charPerLine
    if(nCharLastLine == 0) nCharLastLine <- charPerLine
    output <- character()
    blueLines <- list( # junction lines to denote with highlight coloring
        Genome = 3:4,
        'Genome+' = 3:4,
        Reads = 1  
    )
    for(i in 1:nChunks){
        if(i > 1) output <- c(output, "<br>")        
        start <- 1 + (i - 1) * charPerLine
        for(j in 1:nLines){
            str  <- substr(data[j], start, start + charPerLine - 1)
            str  <- gsub(' ', '&nbsp;', str) 
            str_ <- gsub('~', '', str) 
            if(nchar(str_) > 0){
                col <- if(j %in% blueLines[[input$jxnDisplayType]]) 'blue' else 'black'
                output <- c(
                    output,
                    paste(paste('<span style="color:', col, ';">', sep = ""),
                          str,
                          '</span>', sep = "")
                )                
            }            
        }
    }
    output <- gsub('~', '&nbsp;', output)
    
    # # app/prepend the outermost chromosomal positions
    # firstInPair <- 64
    # svx <- if(bitwAnd(zoom$sv$FLAG1, firstInPair)) {
    #     zoom$sv[,.(RNAME1,PROX_OUT_POS1,RNAME2,PROX_OUT_POS2)]
    # } else {
    #     zoom$sv[,.(RNAME2,PROX_OUT_POS2,RNAME1,PROX_OUT_POS1)]
    # }
    # output <- c(
    #     "<br>",
    #     paste(svx[[1]], format(svx[[2]], big.mark=","), sep=":"),
    #     '&darr;',
    #     output,
    #     paste(paste(rep("&nbsp;", nCharLastLine-1),collapse=""),'&uarr;',sep=""),
    #     paste(svx[[3]], format(svx[[4]], big.mark=","), sep=":")
    # )
    
    # return the hopefully readable sequence block
    paste(output, collapse = "<br>")
})

# ----------------------------------------------------------------------
# define bookmarking actions
# ----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    sampleSet <- bm$input[['sampleSelector-sampleSet']]
    sampleSelector$setSampleSet(sampleSet) 
    if(!is.null(bm$outcomes)) {
        outcomes <<- listToReactiveValues(bm$outcomes)
        sampleSelector$setSelectedSamples(sampleSet, bm$outcomes$samples)
        locationsPlot$settings$replace(bm$outcomes$locationsPlotSettings)
        propertiesPlot$settings$replace(bm$outcomes$propertiesPlotSettings)
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    samples  = sampleSelector$selectedSamples,
    outcomes = reactive({ list(
        samples = sampleSelector$selectedSamples(),
        locationsPlotSettings  = locationsPlot$settings$all_(),
        propertiesPlotSettings = propertiesPlot$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
