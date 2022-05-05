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
    templates = list(file.path(app$sources$suiteGlobalDir, "settings", "svx_filters.yml"), id), 
        file.path(app$sources$suiteGlobalDir, "settings", "svCapture_filters.yml"),
    fade = FALSE
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
# analyze the insertsion in the filtered SVs
#----------------------------------------------------------------------
analyzeInsertion <- function(
    microhomologyLength, minRefWidth, padding,
    JXN_BASES, GEN_REF_1, GEN_REF_2, SIDE_1, SIDE_2
){

    # flip one reference strand for inversions to match jxnSeq assembly
    if(SIDE_1 == SIDE_2){ 
        if(SIDE_1 == "L") GEN_REF_2 <- rc(GEN_REF_2)
                     else GEN_REF_1 <- rc(GEN_REF_1)
    }

    # parse the insertion search sequence as (microhomology)(insertion)(microhomology)
    thisRefWidth <- nchar(GEN_REF_1)
    thisPadding <- (thisRefWidth - 1) / 2
    microhomology1 <- substr(GEN_REF_1, 
                             thisPadding - microhomologyLength + 2,     
                             thisPadding + 1)
    microhomology2 <- substr(GEN_REF_2, 
                             thisPadding + 1, 
                             thisPadding + microhomologyLength)
    searchSeq <- paste0(microhomology1, JXN_BASES, microhomology2)
    rcSearchSeq <- rc(searchSeq)

    # trim the genome references to the mininum available sample padding
    if(padding < thisPadding){
        start <- thisPadding - padding + 1
        stop  <- start + minRefWidth - 1
        GEN_REF_1 <- substr(GEN_REF_1, start, stop)
        GEN_REF_2 <- substr(GEN_REF_2, start, stop)
    }

    # return all matches to the search sequence
    c(
        searchSeq = searchSeq,
        match1   = paste(gregexpr(  searchSeq, GEN_REF_1)[[1]], collapse = ","),
        match1rc = paste(gregexpr(rcSearchSeq, GEN_REF_1)[[1]], collapse = ","),        
        match2   = paste(gregexpr(  searchSeq, GEN_REF_2)[[1]], collapse = ","),
        match2rc = paste(gregexpr(rcSearchSeq, GEN_REF_2)[[1]], collapse = ",") 
    )
}
insertions <- reactive({ # get the SVs with accepted insertions
    svs <- filteredSvs()
    req(svs)
    properties <- settings$Junction_Properties()
    svs <- svs[
        -MICROHOM_LEN >= properties$Min_Insertion_Size$value & 
        -MICROHOM_LEN <= properties$Max_Insertion_Size$value
    ]
    req(svs)
    minRefWidth <- svs[, min(nchar(GEN_REF_1))]    
    padding <- (minRefWidth - 1) / 2
    svs <- cbind(svs, t(svs[, mapply(
        analyzeInsertion, 
        properties$Flanking_Microhomology$value, minRefWidth, padding, 
        JXN_BASES, GEN_REF_1, GEN_REF_2, SIDE_1, SIDE_2
    )]))
    list(
        refWidth = minRefWidth,
        padding = padding,
        svs = svs
    )
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
        insertions$svs[, ':='(
            nCombinations = 4 ** (-MICROHOM_LEN + properties$Flanking_Microhomology$value * 2), 
            nMatches = mapply(function(match1, match1rc, match2, match2rc){
                sum(strsplit(match1,   ",")[[1]] != "-1") + # total number of times search sequence was found
                sum(strsplit(match1rc, ",")[[1]] != "-1") + 
                sum(strsplit(match2,   ",")[[1]] != "-1") + 
                sum(strsplit(match2rc, ",")[[1]] != "-1")
            }, match1, match1rc, match2, match2rc)
        )]
        insertions$svs[, found := nMatches > 0] # whether or not each SV's search sequence was found
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
        templateLocations$addLegend(
            legend = c("cross-junction", "fold-back", "other"),
            col =    c(crossJxn,         foldback,    CONSTANTS$plotlyColors$black),
            x = "top",
            horiz = TRUE,
            x.intersp = 0
        )
        
    }
)

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
        templateLocations$settings$replace(bm$outcomes$templateLocationsSettings)
        templateYield$settings$replace(bm$outcomes$templateYieldSettings)
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
        templateLocationsSettings = templateLocations$settings$all_(),
        templateYieldSettings = templateYield$settings$all_()
    ) }),
    isReady  = reactive({ getStepReadiness(options$source, outcomes) })
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
