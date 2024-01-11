#----------------------------------------------------------------------
# reactive components to plot summary results over all selected samples
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
analyzeInsertionsServer <- function(id, options, bookmark, locks) {
    moduleServer(id, function(input, output, session) {
        ns <- NS(id) # in case we create inputs, e.g. via renderUI
        module <- 'analyzeInsertions' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
settings <- settingsServer( # display settings not stored in the UI, exposed by gear icon click
    id = 'settings',
    parentId = id,
    templates = list(
        file.path(app$sources$suiteGlobalDir, "settings", "svx_filters.yml"), 
        file.path(app$sources$suiteGlobalDir, "settings", "svCapture_filters.yml"),
        id
    ),
    fade = FALSE,
    presets = settingsPresets
)
sampleSelector <- sampleSelectorServer( # selectors to pick one or more samples from a sample set
    id = 'sampleSelector',
    parentId = id
)
outcomes <- reactiveValues()

#----------------------------------------------------------------------
# generate the list of all filtered SVs from all selected samples
#----------------------------------------------------------------------
filteredSvs <- reactive({
    getFilteredSvs(settings, sampleSelector, isCapture = TRUE)
})

#----------------------------------------------------------------------
# analyze the insertions in the filtered SVs
#----------------------------------------------------------------------
insertions <- reactive({
    svs <- filteredSvs()
    svx_analyzeInsertions(svs, settings$Junction_Properties())
})

#----------------------------------------------------------------------
# templated yield plots
#----------------------------------------------------------------------
svCounts <- staticPlotBoxServer(
    'svCounts', 
    points = TRUE,
    margins = TRUE,
    title = TRUE,
    immediate = TRUE,
    create = function(...){
        insertions <- insertions()
        req(insertions)
        count <- insertions$svs[, .N, keyby = MICROHOM_LEN] 
        properties <- settings$Junction_Properties()
        svCounts$initializeFrame(
            xlim = range(-count$MICROHOM_LEN),            
            ylim = c(0, max(count$N) * 1.05),
            xlab = "Insertion Size (bp)",            
            ylab = "# of SVs"
        )       
        abline(h = 0)
        abline(h = c(1, seq(10, 5000, 10)), col = "grey")
        svCounts$addPoints(
            x = -count$MICROHOM_LEN,            
            y = count$N
        )
    }
)
templateYield <- staticPlotBoxServer(
    'templateYield', 
    points = TRUE,
    lines = TRUE,
    margins = TRUE,
    title = TRUE,
    legend = TRUE,
    immediate = TRUE,
    create = function(...){
        insertions <- insertions()
        req(insertions)
        properties <- settings$Junction_Properties()
        searchSpace <- 4 * insertions$refWidth # two genome references searched on both strands
        yield <- insertions$svs[, .(
            nCombinations = nCombinations[1],
            nSvs = .N,           # number of trials ...
            nFound = sum(found), # number of successes ...
            mu = searchSpace * (1 / nCombinations[1]) # Poisson mean, i.e., expected target density ...
        ), keyby = MICROHOM_LEN] # ... by insert size
        yield <- yield[, ":="(
            successRate = nFound / nSvs, # fraction of SVs whose search sequence was found
            trialSuccessProb = 1 - dpois(0, mu) # probably of finding at least one match
        )]
        yield <- yield[, pValue := 1 - pbinom(nFound - 1, nSvs, trialSuccessProb)]
        templateYield$initializeFrame(
            xlim = range(-yield$MICROHOM_LEN),            
            ylim = c(0, min(100, max(yield$successRate, yield$trialSuccessProb) * 110)),
            xlab = "Insertion Size (bp)",            
            ylab = "% Found Templates"
        )  
        templateYield$addLines(
            x = -yield$MICROHOM_LEN,            
            y = yield$trialSuccessProb * 100,
            col = CONSTANTS$plotlyColors$blue
        )     
        pThreshold <- 0.05
        templateYield$addPoints(
            x = -yield$MICROHOM_LEN,            
            y = yield$successRate * 100,
            col = ifelse(yield$pValue < 0.05, 
                         CONSTANTS$plotlyColors$red, 
                         CONSTANTS$plotlyColors$black)
        )
        templateYield$addLegend(
            legend = c(
                paste("p", "<",  pThreshold),
                paste("p", ">=", pThreshold)
            ),
            col = c(
                CONSTANTS$plotlyColors$red, 
                CONSTANTS$plotlyColors$black
            )
        )
    }
)

#----------------------------------------------------------------------
# template locations plot
#----------------------------------------------------------------------
templateLocations <- staticPlotBoxServer(
    'templateLocations', 
    points = TRUE,
    margins = TRUE,
    title = TRUE,
    immediate = TRUE,
    template = stepModuleInfo$analyzeInsertions$locationSettings,
    create = function(...){
        insertions <- insertions()
        req(insertions)

        # initialize shared values
        jxn  <- settings$Junction_Properties()
        dist <- templateLocations$settings$Distance()
        microhomologyLength <- jxn$Flanking_Microhomology$value
        foldback <- CONSTANTS$plotlyColors$blue
        crossJxn <- CONSTANTS$plotlyColors$green
        svs <- insertions$svs

        # calculate the plot points and colors
        getPoints <- function(svId, insSize, match, y, col, multiplier){
            matchPositions <- as.integer(strsplit(match, ",")[[1]])
            if(length(matchPositions) == 1 && matchPositions == -1) return(NULL)
            pos <- matchPositions - insertions$padding - 1 # leftmost position of (uHom)(ins)(uHom) searchSeq
            pos <- ifelse( # convert to the closest base of the insertion itself
                pos < 0, 
                pos + microhomologyLength + insSize - 1, 
                pos + microhomologyLength
            )
            data.table(
                svId = svId,
                distance = abs(pos),
                x = pos,
                y = y,
                col = ifelse(multiplier * pos > 0, col, CONSTANTS$plotlyColors$black)
            )
        }
        d <- do.call(rbind, list(
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match1,   4, crossJxn, -1)
            ),
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match1rc, 3, foldback, -1)
            ),
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match2,   2, crossJxn,  1)
            ),            
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match2rc, 1, foldback,  1)
            )           
        ))

        # for SVs with multiple possible template, prefer the one closest to the junction
        d <- unique(d[order(distance)], by = 'svId')

        # initialize plot frame with genomic molecule tracing and grid
        limit <- dist$Max_Plotted_Distance$value
        templateLocations$initializeFrame(
            xlim = c(-limit, limit),          
            ylim = c(0.5, 6),
            xlab = "Insertion Template Distance from Junction (bp)",            
            ylab = "",
            yaxt = "n",
            bty = "n"
        )
        grid <- dist$Distance_Grid_Spacing$value
        # abline(v = seq(-grid * 1000, grid * 1000, grid), col = CONSTANTS$plotlyColors$grey) 
        for(x in seq(-grid * 1000, grid * 1000, grid)){
            lines(c(x, x), c(0, 4.75), col = CONSTANTS$plotlyColors$grey)
        }
        for(y in 3:4){ # left side of junction, drawn on top two lines
            lines(c(-insertions$padding, 0), c(y, y), col = CONSTANTS$plotlyColors$black, lwd = 2)  
            lines(c(0, insertions$padding),  c(y, y), col = CONSTANTS$plotlyColors$grey)   
        }
        for(y in 1:2){ # right side of junction, drawn on bottom two lines
            lines(c(-insertions$padding, 0), c(y, y), col = CONSTANTS$plotlyColors$grey)  
            lines(c(0, insertions$padding),  c(y, y), col = CONSTANTS$plotlyColors$black, lwd = 2)   
        }        
        lines(c(0, 0), c(1, 4), col = CONSTANTS$plotlyColors$black, lwd = 2)
        text(-limit, 5.5, "Left Side of SV",  pos = 4)
        text(limit,  5.5, "Right Side of SV", pos = 2)

        # add points for the found putative insertion templates
        # positions are the 1st base of the insertion closest to the junction at position 0
        templateLocations$addPoints(
            x = jitter(d$x, amount = 0.25), # jitter to make all point visible      
            y = jitter(d$y, amount = 0.25),
            col = d$col
        )

        # add a point color legend
        cex <- templateLocations$settings$Points_and_Lines()$Point_Size$value
        templateLocations$addLegend(
            legend = c("cross-junction", "fold-back", "other"),
            col =    c(crossJxn,         foldback,    CONSTANTS$plotlyColors$black),
            x = "top",
            horiz = TRUE,
            x.intersp = 0,
            pt.cex = cex * 1.5
        )
    }
)
templateLocations2 <- staticPlotBoxServer(
    'templateLocations2', 
    lines = TRUE,
    margins = TRUE,
    title = TRUE,
    immediate = TRUE,
    template = stepModuleInfo$analyzeInsertions$locationSettings,
    create = function(...){
        insertions <- insertions()
        req(insertions)

        # initialize shared values
        jxn  <- settings$Junction_Properties()
        dist <- templateLocations2$settings$Distance()
        microhomologyLength <- jxn$Flanking_Microhomology$value
        foldback <- CONSTANTS$plotlyColors$blue
        crossJxn <- CONSTANTS$plotlyColors$green
        other    <- CONSTANTS$plotlyColors$black
        cols <- list(foldback = foldback, crossJxn = crossJxn, other = other)
        svs <- insertions$svs

        # calculate the plot points and colors
        getPoints <- function(svId, insSize, match, y, type, multiplier){
            matchPositions <- as.integer(strsplit(match, ",")[[1]])
            if(length(matchPositions) == 1 && matchPositions == -1) return(NULL)
            pos <- matchPositions - insertions$padding - 1 # leftmost position of (uHom)(ins)(uHom) searchSeq
            pos <- ifelse( # convert to the closest base of the insertion itself
                pos < 0, 
                pos + microhomologyLength + insSize - 1, 
                pos + microhomologyLength
            )
            data.table(
                svId = svId,
                distance = abs(pos),
                type = ifelse(multiplier * pos > 0, type, "other")
            )
        }
        d <- do.call(rbind, list(
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match1,   4, "crossJxn", -1)
            ),
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match1rc, 3, "foldback", -1)
            ),
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match2,   2, "crossJxn",  1)
            ),            
            do.call(rbind, 
                mapply(getPoints, svs$SV_ID, -svs$MICROHOM_LEN, svs$match2rc, 1, "foldback",  1)
            )           
        ))

        # for SVs with multiple possible template, prefer the one closest to the junction
        d <- unique(d[order(distance)], by = 'svId')

        # calculate the frequency histogram by type and distance
        limit <- dist$Max_Plotted_Distance$value
        d <- dcast(d, distance ~ type, fun.aggregate = length, fill = 0)
        d <- merge(data.table(distance = 1:limit), d, by = "distance", all.x = TRUE)

        # initialize plot frame with genomic molecule tracing and grid
        plotCols <- c("crossJxn", "foldback", "other")
        templateLocations2$initializeFrame(
            xlim = c(0, limit),          
            ylim = c(0, d[, max(.SD, na.rm = TRUE) * 1.2, .SDcols = plotCols]),
            xlab = "Distance from Junction (bp)",            
            ylab = "# of Insertions",
            # yaxt = "n",
            bty = "n"
        )

        # add location histograms for the found putative insertion templates
        lwd <- templateLocations2$settings$Points_and_Lines()$Line_Width$value
        for(type in plotCols){
            lines(d$distance, d[[type]], col = cols[[type]], lwd = lwd)
        }

        # add a point color legend
        templateLocations$addLegend(
            legend = c("cross-junction", "fold-back", "other"),
            col =    unlist(cols[plotCols]),
            x = "top",
            horiz = TRUE,
            x.intersp = 0,
            pt.cex = lwd * 0.9
        )
    }
)

# ----------------------------------------------------------------------
# summary table of all filtered SVs from all selected samples
# ----------------------------------------------------------------------
foundInsertions <- reactive({ 
    insertions <- insertions()
    req(insertions)
    insertions$svs[found == TRUE] 
})
svsTable <- filteredSvsTableServer(id, input, foundInsertions)

# # ----------------------------------------------------------------------
# # dotplot alignments of junction and genome reference sequences
# # ----------------------------------------------------------------------
# # the one working SV the user has clicked on
# selectedSv <- reactive({ 
#     rowI <- svsTable$rows_selected()
#     if(is.null(rowI) || rowI == 0) return(NULL) # req(rowI)
#     foundInsertions()[rowI]
# })
# # complete supporting molecule evidence on selectedSv()
# svMols <- reactive({ getSVMolecules(sampleSelector, selectedSv) })
# # map of all base positions at and around the selected SV junction
# junctionMap <- reactive({ 
#     map <- getJunctionMap(svMols()) 
#     req(map)
#     startSpinner(session, 'expanding junctionMap combinations') 
#     map$consensus <- getJunctionConsensus(map$text)
#     map$pairs <- expand.grid(
#         x = (map$leftRefI + 1 - 10):(map$rightRefI - 1),
#         y = 1:(length(map$consensus) + map$sv$MICROHOM_LEN)
#     )
#     stopSpinner(session) 
#     map    
# })
# dotPlotBoxServer <- function(id, x, y, xlab, ylab){
#     plot <- staticPlotBoxServer(
#         id, 
#         points = TRUE,
#         margins = TRUE,
#         title = TRUE,
#         immediate = TRUE,
#         create = function(...){
#             map <- junctionMap()
#             req(map)
#             startSpinner(session, paste("dotPlotBoxServer", id)) 
#             xseq <- toupper(sapply(map$pairs$x, function(i) 
#                 paste(map[[x]][i:(i - map$sv$MICROHOM_LEN - 1)], collapse = "") 
#             ))
#             yseq <- toupper(sapply(map$pairs$y, function(i) 
#                 paste(map[[y]][i:(i - map$sv$MICROHOM_LEN - 1)], collapse = "") 
#             ))
#             match   <- xseq ==    yseq
#             matchrc <- xseq == rc(yseq)
#             stopSpinner(session) 
#             plot$initializeFrame(
#                 xlim = range(map$pairs$x),
#                 ylim = range(map$pairs$y),
#                 xlab = xlab,
#                 ylab = ylab
#             )
#             cols <- CONSTANTS$plotlyColors            
#             abline(v = c(map$leftRefI, map$rightRefI), col = cols$grey)
#             abline(h = c(map$leftRefI, map$rightRefI), col = cols$grey)
#             plot$addPoints(
#                 x = map$pairs$x,
#                 y = map$pairs$y,
#                 pch = ifelse(match, ".", NA),
#                 col = ifelse(map$pairs$x == map$pairs$y, cols$grey, cols$blue)
#             )
#             plot$addPoints(
#                 x = map$pairs$x,
#                 y = map$pairs$y,
#                 pch = ifelse(matchrc, ".", NA),
#                 col = ifelse(map$pairs$x == map$pairs$y, cols$grey, cols$red)
#             )
#         }
#     )
# }
# jxn_side1 <- dotPlotBoxServer(
#     'jxn_side1', 
#     "consensus", "GEN_REF_1", 
#     "Junction Sequence", "Genome Reference, Left"
# )
# jxn_side2 <- dotPlotBoxServer(
#     'jxn_side2', 
#     "consensus", "GEN_REF_2", 
#     "Junction Sequence", "Genome Reference, Right"
# )

# ----------------------------------------------------------------------
# define bookmarking actions
# ----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    sampleSet <- bm$input[['sampleSelector-sampleSet']]
    sampleSelector$setSampleSet(sampleSet) 
    if(!is.null(bm$outcomes)) {
        outcomes <<- listToReactiveValues(bm$outcomes)
        sampleSelector$setSelectedSamples(sampleSet, bm$outcomes$samples)
        templateLocations$settings$replace(bm$outcomes$templateLocationsSettings)
        templateLocations2$settings$replace(bm$outcomes$templateLocationsSettings2)
        templateYield$settings$replace(bm$outcomes$templateYieldSettings)
    }
    bookmarkObserver$destroy()
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
        templateLocationsSettings = templateLocations$settings$all_(),
        templateLocationsSettings2 = templateLocations2$settings$all_(),
        templateYieldSettings = templateYield$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
