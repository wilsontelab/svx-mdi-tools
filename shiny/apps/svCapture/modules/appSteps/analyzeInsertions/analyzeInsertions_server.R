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
    templates = list(file.path(app$sources$suiteGlobalDir, "sv_filters.yml"), id),
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
    getFilteredSvs(settings, sampleSelector)
})

#----------------------------------------------------------------------
# analyze the insertsion in the filtered SVs
#----------------------------------------------------------------------
analyzeInsertion <- function(JXN_BASES, GEN_REF_1, GEN_REF_2, SIDE_1, SIDE_2){

    # flip one reference strand for inversion
    if(SIDE_1 == SIDE_2){ 
        if(SIDE_1 == "L") GEN_REF_2 <- rc(GEN_REF_2)
                     else GEN_REF_1 <- rc(GEN_REF_1)
    }

    # parse the search sequence as (2bp microhomology)(insertion)(2bp microhomology)
    refWidth <- nchar(GEN_REF_1)
    faidxPadding <- (refWidth - 1) / 2
    microhomology1 <- substr(GEN_REF_1, faidxPadding,     faidxPadding + 1)
    microhomology2 <- substr(GEN_REF_2, faidxPadding + 1, faidxPadding + 2)
    searchSequence <- paste0(microhomology1, JXN_BASES, microhomology2)

    # return all matches to the search sequence
    c(
        faidxPadding = faidxPadding,
        searchSequence = searchSequence,
        match1   = paste(gregexpr(searchSequence,    GEN_REF_1 )[[1]], collapse = ","),
        match2   = paste(gregexpr(searchSequence,    GEN_REF_2 )[[1]], collapse = ","),
        match1rc = paste(gregexpr(searchSequence, rc(GEN_REF_1))[[1]], collapse = ","),
        match2rc = paste(gregexpr(searchSequence, rc(GEN_REF_2))[[1]], collapse = ",") 
    )
}
insertionSvs <- reactive({ 

    # get the SVs with accepted insertions
    svs <- filteredSvs()
    req(svs)
    properties <- settings$Junction_Properties()
    svs <- svs[
        -MICROHOM_LEN >= properties$Min_Insertion_Size$value & 
        -MICROHOM_LEN <= properties$Max_Insertion_Size$value 
        # & 
        # SV_ID == "6773:1" ########################### CCttaccTT
    ]
    req(svs)    
    svs <- cbind(svs, t(svs[, mapply(analyzeInsertion, JXN_BASES, GEN_REF_1, GEN_REF_2, SIDE_1, SIDE_2)]))
    svs[, faidxPadding := as.integer(faidxPadding)]
    svs
})

#----------------------------------------------------------------------
# templated locations plot
#----------------------------------------------------------------------
templateLocations <- staticPlotBoxServer(
    'templateLocations', 
    legend = TRUE,
    immediate = TRUE,
    points = TRUE,
    margins = TRUE,
    title = TRUE,
    create = function(...){

        insertionsSvs <- insertionSvs()
        req(insertionsSvs)
        # dstr(insertionsSvs)
        # dprint(insertionsSvs[, .(SV_ID, JXN_BASES, MICROHOM_LEN, search)])
        # dprint(dim(insertionsSvs))


        # svRates <- copy(svRates())[order(rev(1:.N))]
        # req(svRates)
        # svRates[, groupFactor := factor(group, unique(group))]
        # xi <- as.integer(svRates$groupFactor)
        templateLocations$initializeFrame(
            xlim = c(-10, 10),            
            ylim = c(-10, 10),
            xlab = "Distance from Junction (bp)",            
            ylab = "Strand"
        )
        # abline(v = seq(0, 1, 0.1), col = "grey80")
        # abline(h = c(0, xi) + 0.5, col = "grey80")
        templateLocations$addPoints(
            x = 1,            
            y = jitter(1, amount = 0.5)
        )
        # axis(
        #     2, 
        #     at = xi, 
        #     labels = svRates$group,
        #     las = 1
        # )
        # # barplot(as.matrix(x), beside = TRUE, horiz = TRUE)
    }
)

#----------------------------------------------------------------------
# templated yield plots
#----------------------------------------------------------------------
svCounts <- staticPlotBoxServer(
    'svCounts', 
    legend = TRUE,
    immediate = TRUE,
    points = TRUE,
    margins = TRUE,
    title = TRUE,
    create = function(...){
        insertionsSvs <- insertionSvs()
        req(insertionsSvs)
        count <- insertionsSvs[, .N, keyby = MICROHOM_LEN] 
        properties <- settings$Junction_Properties()        
        templateYield$initializeFrame(
            xlim = range(-count$MICROHOM_LEN),            
            ylim = c(0, max(count$N) * 1.05),
            xlab = "Insertion Size (bp)",            
            ylab = "# of SVs"
        )       
        abline(h = 0)
        abline(h = c(1:5, seq(5, 100, 5)), col = "grey")
        templateYield$addPoints(
            x = -count$MICROHOM_LEN,            
            y = count$N
        )
    }
)
templateYield <- staticPlotBoxServer(
    'templateYield', 
    legend = TRUE,
    immediate = TRUE,
    points = TRUE,
    margins = TRUE,
    title = TRUE,
    create = function(...){
        insertionsSvs <- insertionSvs()
        req(insertionsSvs)
        insertionsSvs[, ':='(
            searchSpace = 4 * (faidxPadding * 2 + 1), # two genome references searched on both strands
            nCombinations = 4 ** (-MICROHOM_LEN + 4), # inserted sequence plus two 2bp flanking microhomologies
            nTemplates = mapply(function(match1, match1rc, match2, match2rc){
                sum(strsplit(match1,   ",")[[1]] != "-1") + # total number of times search sequence was found
                sum(strsplit(match1rc, ",")[[1]] != "-1") + 
                sum(strsplit(match2,   ",")[[1]] != "-1") + 
                sum(strsplit(match2rc, ",")[[1]] != "-1")
            }, match1, match1rc, match2, match2rc)
        )]
        insertionsSvs[, found := nTemplates > 0] # whether or not each SV's search sequence was found
        yield <- insertionsSvs[, .(
            nCombinations = nCombinations[1],
            nSvs = .N,           # number of trials ...
            nFound = sum(found), # number of successes ...
            mu = searchSpace[1] * (1 / nCombinations[1]), # Poisson mean, i.e., expected target density ...
            successRate = sum(found) / .N # fraction of SVs whose search sequence was found, i.e., positive trials ...
        ), keyby = MICROHOM_LEN]      # ... by insert size

        yield <- yield[, trialSuccessProb := 1 - dpois(0, mu)] # probably of finding at least one match
        yield <- yield[, pValue := pbinom(nFound - 1, nSvs, trialSuccessProb, lower.tail = FALSE)] # nFound or more


        dprint(yield)

        # dbinom(nFound, N, trialSuccessProb)
        # pbinom(nFound, N, trialSuccessProb, lower.tail = TRUE)

        properties <- settings$Junction_Properties()        
        templateYield$initializeFrame(
            xlim = range(-yield$MICROHOM_LEN),            
            ylim = c(0, min(100, max(yield$successRate, yield$trialSuccessProb) * 110)),
            xlab = "Insertion Size (bp)",            
            ylab = "% Found Templates"
        )
        # templateYield$addLines(
        #     x = -yield$MICROHOM_LEN,            
        #     y = yield$expectedHitRate * 100,
        #     col = "blue"
        # )    
        templateYield$addLines(
            x = -yield$MICROHOM_LEN,            
            y = yield$trialSuccessProb * 100,
            col = "blue"
        )     
        templateYield$addPoints(
            x = -yield$MICROHOM_LEN,            
            y = yield$successRate * 100,
            col = ifelse(yield$pValue < 0.05, "red3", "black")
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
