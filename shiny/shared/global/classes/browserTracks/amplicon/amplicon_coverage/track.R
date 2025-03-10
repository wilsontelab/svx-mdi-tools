#----------------------------------------------------------------------
# amplicon_coverage trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
amplicon_coverageBuffer  <- reactiveVal(NULL)
amplicon_highCopyRegions <- reactiveVal(NULL)
amplicon_hcrHighlightsFn <- function(track, reference, coord, sampleName, sample){
    hcr <- amplicon_highCopyRegions()
    if(!isTruthy(hcr) || is.null(hcr[[sampleName]])) return(NULL)
    if(coord$chromosome == "all") hcr$highCopyRegions[[sampleName]][
        genomeStart <= coord$end & 
        coord$start <= genomeEnd,
        .(
            x1 = genomeStart,
            x2 = genomeEnd,
            color = CONSTANTS$plotlyColors$RED
        )
    ] else hcr$highCopyRegions[[sampleName]][
        chrom == coord$chromosome &
        start <= coord$end & 
        coord$start <= end,
        .(
            x1 = start,
            x2 = end,
            color = CONSTANTS$plotlyColors$RED
        )
    ]
}

# constructor for the S3 class
new_amplicon_coverageTrack <- function(...) {
    new_svx_coverageTrack(..., click = TRUE)
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.amplicon_coverageTrack <- function(...) showTrackSamplesDialog(...)

# build method for the S3 class; REQUIRED
build.amplicon_coverageTrack <- function(track, reference, coord, ...){
    amplicon_coverageBuffer(list(
        reference = reference,
        coord = coord
    ))
    build.svx_coverageTrack(
        track, reference, coord,
        ..., 
        loadFn = svWGS_loadSampleCoverage,
        highlightsFn = amplicon_hcrHighlightsFn
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
# regionI indicates the region plot the user interacted with, the value must be passed to app$browser$jumpToCoordinates, etc.
click.amplicon_coverageTrack <- function(track, click, regionI){
    buffer <- amplicon_coverageBuffer()
    minAmpliconCN <- click$coord$y
    req(buffer, minAmpliconCN)
    selectedSources <- getSourcesFromTrackSamples(track$settings$items())
    power <- 0
    medianPloidy <- 2

    # TODO: here or elsewhere, resolve chains and cycles in junction graph
    # TODO: here or elsewhere, find maxima in a coverage distribution over all high copy regions
    #       this cannot be done by taking median CN over regions, as regions may have multiple CN segments
    # runs of these discrete CN become graph segments, for which we expect terminal jxns
    # so:
    #   find CN maxima
    #   find runs of CN that match the same maximum, i.e., CN state
    #   infer that runs with the same CN_state are on the same amplicon
    #   HOWEVER, it is NOT that simple, since step-down regions occur
    #   need to fuse CN_segments that are:
    #       in the same highCopyRegion
    #       share the same CN state
    #       are interrupted by a region of (higher?) CN

    do.call(c, lapply(names(selectedSources), function(sourceId_){
        sapply(selectedSources[[sourceId_]]$Sample_ID, function(sampleName_){
            coverage <- svWGS_loadSampleCoverage(sourceId, sampleName_)[[paste0("x", power)]]
            # x0 <- coverage[, .(
            #     chrom = chrom,
            #     start = start + 1,
            #     genomeStart = getSignedNode(chroms$chromSizes, unlist(chroms$chromIndex[chrom]), start, "+", 1),
            #     coverage = ifelse(excluded > 0, NA, sampleCoverage),
            #     gc = gc,
            #     excluded = excluded / (start[2] - start[1])
            # )]            
            binSize <- coverage[1:2, diff(start)] 
            normalized <- if(is.null(app$normalizeGC)) NULL else app$normalizeGC$getBinNormalizedCN(
                sourceId_, 
                sampleName_,
                buffer$reference, # this and below only used to parse the genome, not to filter by browser window
                buffer$coord, 
                coverage$gc, 
                coverage$coverage
            )
            coverage$CN <- if(is.null(normalized)) coverage$coverage / median(coverage$coverage, na.rm = TRUE) * medianPloidy else normalized
            highCopyRegions <- coverage[,
                .(
                    start = min(start - 1, na.rm = TRUE),
                    end = max(start - 1 + binSize, na.rm = TRUE),
                    genomeStart = min(genomeStart, na.rm = TRUE),
                    genomeEnd = max(genomeStart + binSize, na.rm = TRUE),
                    coverage = median(coverage, na.rm = TRUE),
                    CN = median(CN, na.rm = TRUE),
                    isAmplified = CN[1] > minAmpliconCN
                ),
                by = .(chrom, rleidv(CN > minAmpliconCN))
            ]
            CN_states <- coverage[!is.na(CN) & CN > minAmpliconCN, {
                CN_maxima <- getLocalMaxima(CN) # must return peak value and SD
                CN_states <- solveMixtureModel(CN, CN_maxima)
                .(
                    cnState = 1:length(CN_maxima),
                    CN = CN_states$mean,
                    sd = CN_states$sd
                )
            }]
            coverage[, CN_state := solveSimpleHMM(CN, CN_states)]
            list(
                highCopyRegions = highCopyRegions,
                CN_states = NULL
            )
        }, simplify = FALSE, USE.NAMES = TRUE)
    })) %>% amplicon_highCopyRegions()
}

# expand method for the S3 class
# one expansion image can be shown per region, with same width as the main plots
# regionI must be passed to app$browser$expandingTrack
expand.amplicon_coverageTrack <- function(track, reference, coord, layout, regionI){
    # build a track image the same way as for the main track build
    # typically, need to pass data from "build" to "expand" via a variable scoped to the app
    # expansion tracks may follow the browser genome coordinate, or have an entirely different layout
}

# expand2 method for the S3 class
# a single output at browser bottom, replaced on call to any expand2 function
expand2.amplicon_coverageTrack <- function(track, browser, selectedRowData){
    # return a tagList of arbitrary UI content in response to an expansion table click
}

# method for the S3 class to populate one or more trackNav inputs above the browser output
# only one navigation set is shown per track, your navigation should decide how to handle multiple regions
navigation.amplicon_coverageTrack <- function(track, session, id, browser){

    # initialize the trackNavs, including input observers to handle user actions
    # initTrackNav will fail silenty if setting Track/Show_Navigation is not set or =="hide"
    navName1 <- initTrackNav(track, session, "inputName1", function(inputValue1){
        # do work as needed based on the input value, e.g., make a call to app$browser$jumpToCoordinates()
    })
    navName2 <- initTrackNav(track, session, "tableName2") # table reactive functions are provided below
    # etc.

    # as needed, create a reactive with the data to enumerate in the trackNavs
    trackNavData <- reactive({
        data.table( # use track settings or other information to populate dynamically
            chrom = c("chr1","chr2","chr3"),
            start = 1e8,
            end = 2e8
        )
    })

    # return the associated navigation UI elements
    tagList(
        trackNavInput(
            track, 
            session, 
            navName1, # the name as provided by initTrackNav
            radioButtons, # any valid Shiny UI input function
            label = "Input Name 1", # all further arguments must be named and are passed to the UI function
            choices = c("foo", "bar"),            
            selected = "bar", # the default value
            inline = TRUE,
            width = "150px"
            # add other argument to pass to the Shiny UI input function
        ),
        trackNavTable(
            track, 
            session, 
            browser$id,
            navName2, # the name as provided by initTrackNav
            tableData = trackNavData, # populate a table based on track settings, etc.
            actionFn = function(selectedRow){
                req(selectedRow)
                d <- trackNavData()[selectedRow]
                # do other work as needed based on the input value (e.g., open a modal, navigate)
                # you should honor the value of trackNavCanNavigate(track) and trackNavCanExpand(track)
                # the easiest way is to use handleTrackNavTableClick(regionI, track, chrom, start, end, expandFn)
            }
            # add other argument to pass to bufferedTableServer, but omit "selection" and "width"
        )
        # etc.
    )
}
