
#----------------------------------------------------------------------
# plot_data.R controls the data visualization actions
#----------------------------------------------------------------------

# some common plot values
maxJxnDist <- 1000

# global data objects
SVs        <- NA 
molecules  <- NA
alignments <- NA
 
# mapping of toDistribute to plot actions
plotActions <- list(
    insSize = list(
        level = 'molecule',
        columns = 'insSize',
        xlab = 'Molecule Insert Size (bp)',
        plotType = 'ecdf',
        xlim = c(0,400)
    ),
    #readPairDepth = list(
    #    level = 'molecule',
    #    columns = 'readPairDepth',
    #    xlab = 'Reads Pairs per Strand',
    #    plotType = 'ecdf',
    #    xlim = c(0,50)
    #),
    #---------------------------------
    mapQ = list(
        level = 'alignment',
        #columns = c('MAPQ1','MAPQ2'),
        columns = 'mapQ',
        xlab = 'Average Read Mapping Quality',
        plotType = 'bar',
        categories = 0:60
    ),    
    #sharedProper = list(
    #    level = 'molecule',
    #    columns = 'SHARED_PROPER',
    #    xlab = 'Ends Shared With Proper Molecules',
    #    plotType = 'bar',
    #    categories = 0:2
    #),      
    #umiUsage = list(
    #    level = 'molecule',
    #    columns = c('UMI1','UMI2'),
    #    xlab = 'UMI Number',
    #    plotType = 'bar'
    #),
    #duplex = list(
    #    level = 'molecule',
    #    columns = 'IS_DUPLEX',
    #    xlab = 'Molecule is Duplex',
    #    plotType = 'bar',
    #    categories = 0:1,
    #    labels = c('Strand Orphan','Duplex')
    #),
    #---------------------------------
    svSize = list(
        level = 'sv',
        columns = 'svSize',
        xlab = 'SV Size, log10(bp)',
        plotType = 'ecdf',
        xlim = log10(c(1e2,1e7)),
        log = TRUE
    ),
    fracDuplex = list(
        level = 'sv',
        columns = 'fracDuplex',
        xlab = 'Fraction Duplex Molecules per SV',
        plotType = 'ecdf',
        xlim = c(0,1)
    ),
    fracSharedProper = list(
        level = 'sv',
        columns = 'FRAC_SHARED_PROPER',
        xlab = 'Fraction Ends Shared With Proper Molecules',
        plotType = 'ecdf',
        xlim = c(0,1)
    ),
    readPairDepth = list(
        level = 'sv',
        columns = 'NET_STRAND_COUNT',
        xlab = 'Total read pairs, both strands',
        plotType = 'ecdf',
        xlim = c(0,30)
    ),
    #---------------------------------
    svType = list(
        level = 'sv',
        columns = 'JXN_TYPE',
        xlab = 'Structural Variant Type',
        plotType = 'bar',
        categories = c('D','I','L','T'),
        labels = c('Dup','Inv','Del','Trans')
    ),
    nMol = list(
        level = 'sv',
        columns = 'N_TOTAL',
        xlab = 'Number of Molecules per SV',
        plotType = 'bar',
        categories = 1:250
    ),
    uHomLen = list(
        level = 'sv',
        columns = 'uHomLen',
        xlab = 'Microhomology Length',
        plotType = 'bar',
        categories = -30:30,
        expected = sapply(-30:30, function(x) {
            if(x<0) NA
            else (x + 1) * (1/4)**x * (3/4)**2
        })
    ),
    nSamples = list(
        level = 'sv',
        columns = 'nSamples',
        xlab = 'Number of Samples with a Matching SV',
        plotType = 'bar',
        categories = 1:nSamples
    )
)
 
# retrieve SV data
eventCount_ <- reactiveValues(SVs=0, molecules=0)
setDistSVs <- function(){
    reportProgress('setDistSVs')
    
    # get all matching SVs
    SVs <<- getFilteredSVs(c(
        'SAMPLE','SV_ID','SV_SIZE', 
        'MICROHOM_LEN',
        'JXN_TYPE','JXN_BASES',
        'MICROHOM_LEN','POS_1','POS_2',
        'FRAC_SHARED_PROPER','N_DUPLEX','N_TOTAL',
        'IS_MERGED','SEQ_LEN','UMI','MAPQ'
    ))
    
    # calculate SV derived values
    SVs$uHomLen <<- ifelse(SVs[,'JXN_BASES'] == "0", NA, SVs$MICROHOM_LEN)
    SVs$svSize <<- ifelse(SVs[,'JXN_TYPE'] == "T", NA, SVs$SV_SIZE)
    SVs$fracDuplex <<- SVs$N_DUPLEX / SVs$N_TOTAL
    if(input$toDistribute == 'nSamples'){
        SVs$nSamples <<- if(!is.null(compareInfo)){ # number of samples claiming the same SV (via compare action)
            SV_keys <- paste(SVs$SAMPLE, SVs$SV_ID)
            CI_keys <- paste(compareInfo$SAMPLE, compareInfo$SV_ID)
            sapply(SV_keys, function(key) compareInfo[CI_keys==key,'N_MATCHING_SAMPLES'])
        } else {
            NA
        }
    }

    # update counts
    eventCount_$SVs <- nrow(SVs)
    eventCount_$gapSplit <- sum(SVs$N_GAP_SPLIT)
    eventCount_$molecules <- sum(SVs$N_TOTAL)
    eventCount_$alignments <- length(unlist(strsplit(SVs$MAPQ, '[,:]')))
}

# retrieve molecule data
setDistSVMolecules <- function(){
    reportProgress('setDistSVMolecules')
    molecules <<- data.frame(
        insSize = unlist(mapply(function(isMerged, seqLen){
            #isMerged <- strsplit(isMerged, ',')[[1]]
            #seqLen   <- strsplit(seqLen, ',')[[1]]
            strsplit(seqLen, ',')[[1]]
        }, SVs$IS_MERGED, SVs$SEQ_LEN))
    )
        
    ## define the molecule-level information we need
    #molecules <<- data.frame(
    #    SV_ID = integer(),
    #    UMI1 = integer(),     
    #    UMI2 = integer(),     
    #    IS_DUPLEX = integer(),     
    #    STRAND_COUNT1 = integer(),
    #    STRAND_COUNT2 = integer(),
    #    IS_MERGED = integer(),
    #    SHARED_PROPER = integer(),   
    #    #-------------
    #    JXN_CLASS = character(),          # junction-level information (class = GAP/SPLIT, type=TDIL)
    #    JXN_TYPE = character(),
    #    #-------------
    #    POS1 = integer(), 
    #    OUT_POS1 = integer(),   
    #    MAPQ1 = integer(),   
    #    #-------------
    #    POS2 = integer(), 
    #    OUT_POS2 = integer(),   
    #    MAPQ2 = integer(),   
    #    stringsAsFactors=FALSE,
    #    #-------------
    #    MICROHOM_LEN = integer(),
    #    PROX_JXN_POS1 = integer(),
    #    PROX_JXN_POS2 = integer(),
    #    WAS_SEQUENCED = logical(),
    #    #-------------
    #    MOL_OUT_POS1 = integer(),
    #    MOL_OUT_POS2 = integer()
    #)
    #molCols <- colnames(molecules)
    #
    ## get all molecules that were found as evidence supporting the filtered SVs
    #for (smp in unique(SVs$SAMPLE)){
    #    svIds <- SVs[SVs$SAMPLE==smp,'SV_ID']
    #    data <- getSvMoleculesFromList(smp, svIds)
    #    
    #    # add required SV-level values to each molecule
    #    for(svId in svIds){ 
    #        dataRow <- which(data$SV_ID == svId)
    #        svRow <- which(SVs$SAMPLE == smp & SVs$SV_ID == svId)
    #        data[dataRow,'SV_ID'] <- svId
    #        data[dataRow,'MICROHOM_LEN']  <- SVs[svRow,'MICROHOM_LEN']
    #        data[dataRow,'POS_1'] <- SVs[svRow,'POS_1']
    #        data[dataRow,'POS_2'] <- SVs[svRow,'POS_2']
    #        data[dataRow,'WAS_SEQUENCED'] <- SVs[svRow,'JXN_BASES'] != "0"   
    #    }
    #    molecules <<- rbind(molecules, data[,molCols])
    #}
    #
    ## calculate molecule derived values
    #nReadPairs <- as.integer(molecules$STRAND_COUNT1) + as.integer(molecules$STRAND_COUNT2)
    #nStrands   <- as.integer(molecules$IS_DUPLEX) + 1
    #molecules$readPairDepth <<- nReadPairs / nStrands
    #molecules$insSize <<- mapply(function(mop1, op1, pjp1,
    #                mop2, op2, pjp2,
    #                seq, uHomLen){
    #    if(seq){
    #        if(mop1 == 0) mop1 <- op1
    #        if(mop2 == 0) mop2 <- op2
    #        d1 <- max(mop1,pjp1) - min(mop1,pjp1) + 1 - uHomLen / 2
    #        d2 <- max(mop2,pjp2) - min(mop2,pjp2) + 1 - uHomLen / 2
    #        d <- d1 + d2
    #        if(d <= maxJxnDist) d else NA
    #    } else {
    #        NA
    #    }
    #    
    #}, molecules$MOL_OUT_POS1, molecules$OUT_POS1, molecules$PROX_JXN_POS1,
    #   molecules$MOL_OUT_POS2, molecules$OUT_POS2, molecules$PROX_JXN_POS2,
    #   molecules$WAS_SEQUENCED, molecules$MICROHOM_LEN)
    #
    ## update counts
    #eventCount_$molecules <- nrow(molecules)
}


# retrieve alignment data
setDistSVAlignments <- function(){
    reportProgress('setDistSVAlignments')
    alignments <<- data.frame(
        mapQ = unlist(sapply(SVs$MAPQ, function(mapqs){
            as.integer(strsplit(mapqs, '[,:]')[[1]])
        }))
    )
}

# plot of fragment sizes and junction positions to judge evidence uniqueness
getValuesVector <- function(cols=NULL, isCharacter=FALSE){
    pa <- plotActions[[input$toDistribute]]
    if(is.null(cols)) cols <- pa$columns
    values <- if(pa$level == 'molecule'){
        sapply(cols, function(col) molecules[[col]] )
    } else if(pa$level == 'alignment'){
        sapply(cols, function(col) alignments[[col]] )
    } else {
        sapply(cols, function(col) SVs[[col]] ) 
    }
    values <- if(isCharacter) as.character(values) else as.numeric(values)
    if(!is.null(pa$log)) values <- log10(values)
    values[!is.na(values)]   
}
#getXlim <- function(values){
#    xMin <- if(input$xMin == '') min(values) else as.numeric(input$xMin)
#    xMax <- if(input$xMax == '') max(values) else as.numeric(input$xMax)
#    c(xMin, xMax)
#}
makeDistPlot <- function(){
    reportProgress('makeDistPlot')
    
    # get the filtered SVs (do first so counts show)
    if(input$sample_1 == '-' | input$sample_1 == '') return( plot(1:10, typ="n"))
    setDistSVs()
    
    # get the plot instructions    
    if(input$toDistribute == 'none') return( plot(1:10, typ="n"))    
    pa <- plotActions[[input$toDistribute]]
    
    # get the associated molecules
    if(pa$level == 'molecule') setDistSVMolecules()

    # get the associated alignments
    if(pa$level == 'alignment') setDistSVAlignments()
    
    # manual ecdf, nicer
    if(pa$plotType == 'ecdf'){
        values <- sort( getValuesVector() )
        n <- length(values)
        plot(values, (1:n)/n, type='s',
             xlim=pa$xlim, ylim=c(0, 1),
             xlab = pa$xlab,
             ylab = 'Cumulative Frequency')
        abline(v=median(values))
        abline(h=0.5)
        
    # categorical bar plot
    } else if(pa$plotType == 'bar'){
        values <- getValuesVector(isCharacter=is.character(pa$categories))
        categories <- if(is.null(pa$categories)) unique(values) else pa$categories
        values <- c(pa$categories, values[values %in% categories])
        agg <- aggregate(values, list(values), length)
        agg[[2]] <- agg[[2]] - 1
        y <- agg[[2]]/sum(agg[[2]])
        #ylim <- c(0, if(!is.null(pa$expected)) max(pa$expected, y, na.rm=TRUE) else max(y, na.rm=TRUE))
        bp <- barplot(y, space=0.1,
                names.arg=if(is.null(pa$labels)) agg[[1]] else pa$labels,
                xlab=pa$xlab,
                ylab='Frequency')
        
        str(bp)
        
        if(!is.null(pa$expected)) points(bp, pa$expected, typ="b", pch=16, col="red")
    }
}

