
#----------------------------------------------------------------------
# data_sources_session.R defines data objects and common functions to load sample data
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# session-specific object for storing sample SV data
# known sample data does NOT reload, but new samples will
#----------------------------------------------------------------------
sampleData <- list()

#----------------------------------------------------------------------
# information specific to the selected project
#----------------------------------------------------------------------
# all project must have a project-level server config file 'project.info.txt'
# they may have a 'compare' file with uniquness of SVs across project samples
projectInfo <- list()
compareInfo <- NULL
loadProjectInfo <- function(project){
    reportProgress('setProjects')
    projectInfo <<- list()
    compareInfo <<- NULL
    if(project != '-'){
        
        # load project info
        infoFile <- paste(getProjectDir(project), 'project.info.txt', sep="/")
        pi <- read.table(infoFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
        for(i in 1:nrow(pi)){
            projectInfo[[ pi[i,'name'] ]] <<- pi[i,'value']
        }
        projectInfo$READ_LENGTH <<- as.numeric(projectInfo$READ_LENGTH)
        projectInfo$REGION_PADDING <<- as.numeric(projectInfo$REGION_PADDING)
        projectInfo$PROJECT <<- project

        # load compare data if present
        compareRDataFile <- getProjectCompareFile(projectInfo, 'RData')
        if(file.exists(compareRDataFile)){
            load(compareRDataFile)
            compareInfo <<- compareInfo
        } else {
            compareGzFile <- getProjectCompareFile(projectInfo)
            if(file.exists(compareGzFile)){
                compareInfo <<- read.table(
                    pipe(paste('zcat', compareGzFile)),
                    header = FALSE,
                    sep = "\t",
                    col.names = names(compare$structural_variants),
                    colClasses = compare$structural_variants,
                    stringsAsFactors = FALSE
                )
                save(compareInfo, file=compareRDataFile)
            } 
        }    
    }
}

# capture target regions and establish indexing values
captureTargets <- data.frame()
maxI <- 0
nCaptureTargets <- 0
captureTargetNames <- c('-')
loadCaptureTargets <- function(){
    reportProgress('loadCaptureTargets')
    captureTargets <<- data.frame()
    maxI <<- 0
    nCaptureTargets <<- 0
    captureTargetNames <<- c('-')
    if(is.null(projectInfo$CAPTURE_BED)) return()
    if(!is.na(projectInfo$CAPTURE_BED)){
        captureTargetsBed <- getCaptureTargetsBed(projectInfo) 
        d <- read.table(captureTargetsBed, sep="\t", header=FALSE, stringsAsFactors=FALSE)
        d <- d[,1:4]
        names(d) <- c('chrom','start','end','region')        
        nCaptureTargets <<- nrow(d)
        captureTargetNames <<- c('-', sort(unique(d$region)))           
        d <- rbind(d, data.frame( # add a blank capture target
            chrom  = NA,
            start  = NA,
            end    = NA,
            region = NA
        ))
        d$length  <- d$end - d$start # BED half-open format coordinates
        d$center  <- round(d$start + d$length / 2, 0)
        d$paddedStart  <- d$start - projectInfo$REGION_PADDING
        d$paddedEnd    <- d$end   + projectInfo$REGION_PADDING
        d$paddedLength <- d$length + 2 * projectInfo$REGION_PADDING
        plcs <- cumsum(d$paddedLength)
        d$iOffset <- c(0, plcs[1:(nCaptureTargets-1)], NA)
        maxI <<- plcs[nCaptureTargets]
        #for(i in 1:nCaptureTargets){
        #    d[i,'iOffset'] <- maxI
        #    maxI <<- maxI + d[i,'paddedLength']
        #}
        captureTargets <<- d         
    }
}

# the set of expected CNVs based on clone composition
expectedCnvs <- data.frame()
loadExpectedCnvs <- function(){
    reportProgress('loadExpectedCnvs')
    expectedCnvs <<- data.frame()
    if(is.null(projectInfo$EXPECTED_CNVS_BED)) return()
    if(!is.na(projectInfo$EXPECTED_CNVS_BED)){
        expectedCnvsBed <- getExpectedCnvsBed(projectInfo)
        d <- read.table(file=expectedCnvsBed, sep="\t", header=TRUE, stringsAsFactors=FALSE)
        d$ctI <- sapply(d$Gene, function(region) which(captureTargets$region==region))
        d$startI <- d$Start - captureTargets[d$ctI,'start'] + projectInfo$REGION_PADDING + captureTargets[d$ctI,'iOffset']
        d$size   <- abs(d$End - d$Start)
        d$center <- d$startI + d$size / 2
        expectedCnvs <<- d
    }
}

#----------------------------------------------------------------------
# load and pre-configure all SV data for a given sample
#----------------------------------------------------------------------
# support functions to transform data
getSizePercentile <- function(sample, percentile){
    sizes <- sampleData[[sample]]$insertSizeDistribution
    sizes[sizes$CUM_FREQ >= percentile/100,'SIZE'][1]
}

# function for loading a sample, called by plot_data function
loadSamples <- function(samples=NULL){
    reportProgress('loadSamples')
    if(is.null(samples)) samples <- getUserSamples()
    if(length(samples) == 0) return()
    for (sample in samples){
        if(is.null(sampleData[[sample]]) & sample != "-" & sample != "" & !is.null(sample) & !is.na(sample)){
            assign('isLoadingSample', TRUE, envir=.GlobalEnv)            
            sampleRDataFile <- getSampleRDataFile(projectInfo, sample)
            if(file.exists(sampleRDataFile)){
                load(sampleRDataFile)
                sampleData[[sample]] <<- savedSampleData
            } else {
                sampleData[[sample]] <<- list(
                    insertSizeDistribution = read.table(
                        getExtractFile(projectInfo, sample, 'insertSizeDistribution.txt'),
                        header = FALSE,
                        sep = "\t",
                        col.names = names(extract$insertSizeDistribution),
                        colClasses = extract$insertSizeDistribution
                    ) 
                )                
                sampleData[[sample]]$p5  <<- getSizePercentile(sample, 5)
                sampleData[[sample]]$p50 <<- getSizePercentile(sample, 50)
                sampleData[[sample]]$p95 <<- getSizePercentile(sample, 95)
                sampleData[[sample]]$p99 <<- getSizePercentile(sample, 99)                

                findRDataFile <- getFindRDataFile(projectInfo, sample)
                load(findRDataFile)
                svTable <- svTable[
                    CHROM_1 %in% CHROMS &
                    CHROM_2 %in% CHROMS
                ] # only use canonical chromosomes                
                svTable <- addPlotValues(svTable, 'POS_', sample)
                setkey(svTable, SV_ID)
                setindex(svTable, 'SV_ID')
                sampleData[[sample]]$svTable <<- svTable
                savedSampleData <- sampleData[[sample]]
                save(savedSampleData, file=sampleRDataFile)
            }
            assign('isLoadingSample', FALSE, envir=.GlobalEnv)
        }
    }
    samples
}

#----------------------------------------------------------------------
# functions for filtering and retrieving specific data based on user's selections
#----------------------------------------------------------------------
# the samples selected by the user
getUserSamples <- function(){
    reportProgress('getUserSamples')
    samples <- sapply(1:nSamples, function(i) input[[paste0('sample_', i)]])
    unique(samples[samples != "-" & samples != ""])
}

# get and filter required SV data for a multi-SV plot or table
addCommonFilters <- function(boo, sv, sample){
    reportProgress('addCommonFilters')
    # filters common to Targets and Peaks
    if(input$duplexFilter != 'all'){ # whether an SV has at least one duplex molecule
        boo <- boo & if(input$duplexFilter == 'yes') sv$N_DUPLEX > 0 else sv$N_DUPLEX == 0
    }
    if(input$netReadPairFilter != 'all'){ # number of unique evidence molecules for SV
        boo <- boo & switch(input$netReadPairFilter,
            '1'  = sv$NET_STRAND_COUNT == 1,
            '>3' = sv$NET_STRAND_COUNT > 3,
            '>10' = sv$NET_STRAND_COUNT > 10,
            TRUE           
        )
    }
    #if(input$moleculeCountFilter != 'all'){ # number of unique evidence molecules for SV
    #    boo <- boo & switch(input$moleculeCountFilter,
    #        '1'  = sv$N_GAP_SPLIT == 1,
    #        '>1' = sv$N_GAP_SPLIT > 1,
    #        '>3' = sv$N_GAP_SPLIT > 3,
    #        TRUE           
    #    )
    #}
    if(input$moleculeCountFilter != 'all'){ # number of unique evidence molecules for SV
        boo <- boo & switch(input$moleculeCountFilter,
            '1'  = sv$N_TOTAL == 1,
            '>1' = sv$N_TOTAL > 1,
            '>3' = sv$N_TOTAL > 3,
            TRUE           
        )
    }
    if(input$sampleCountFilter != 'all' &
       !is.null(sample) &
       !is.null(compareInfo)){ # number of samples claiming the same SV (via compare action)
        nms <- compareInfo[compareInfo$SAMPLE==sample,
                           c('SV_ID','N_MATCHING_SAMPLES')]
        nms_id <- switch(input$sampleCountFilter,
            '1'  = nms[nms$N_MATCHING_SAMPLES == 1,'SV_ID'],
            '>1' = nms[nms$N_MATCHING_SAMPLES  > 1,'SV_ID'],
            TRUE           
        )
        boo <- boo & sv$SV_ID %in% nms_id
    }  
    if(input$fracProperEnds != 'all'){ # extent to which SV molecules shared ends with proper molecules
        boo <- boo & !is.nan(sv$FRAC_SHARED_PROPER) & switch(input$fracProperEnds,
            '<1' = sv$FRAC_SHARED_PROPER < 1,
            '>=1' = sv$FRAC_SHARED_PROPER >= 1,
            TRUE
        )
    }
    if(input$svTypeFilter != 'all'){ # SV type, deletion, inversion etc.
        boo <- boo & sv$JXN_TYPE == input$svTypeFilter
    }
    if(!is.null(input$minSvSizeFilter) && input$minSvSizeFilter != ''){ 
        boo <- boo & sv$size >= as.numeric(input$minSvSizeFilter)
    }
    if(!is.null(input$maxSvSizeFilter) && input$maxSvSizeFilter != ''){ 
        boo <- boo & sv$size <= as.numeric(input$maxSvSizeFilter)
    }
    boo
}
getFilteredSVs <- function(outColumns, samples=NULL){
    reportProgress('getFilteredSVs')
    samples <- loadSamples(samples)
    svs <- NULL
    outColumns <- union(
        outColumns,
        c('ct1','ct2','size', # be sure that addPlotValues function sets ct1, ct2 and N_GAP_SPLIT
          'N_DUPLEX','NET_STRAND_COUNT','N_TOTAL','N_GAP_SPLIT',
          'TARGET_CLASS','JXN_TYPE','JXN_SEQ','JXN_BASES',
          'SV_ID')
    )
    for (sample in samples){
        sv <- sampleData[[sample]]$svTable[,..outColumns]
        boo <- rep(TRUE, nrow(sv))
        boo <- addCommonFilters(boo, sv, sample)
        boo <- addSpecificFilters(boo, sv)
        if(is.null(svs)) svs <- sv[boo,] else svs <- rbind(svs, sv <- sv[boo,])
    }
    svs
}
 
# retrieve a portion of an indexed file for a single SV
getSvIndexed <- function(sample, svId, fileType, offset, length){
    reportProgress('getSvIndexed')
    if(length <= 1) return(NULL) # expect no data
    thread <- unlist(strsplit(as.character(svId), '\\.')[[1]])[2]
    file <- getIndexedFile(projectInfo, sample, fileType, thread) 
    conn <- file(file, open = "rb")
    seek(conn, where=offset, origin="start", rw="r")
    data <- readChar(conn, length)
    close(conn)
    strsplit(data, "\n")[[1]]
}

# load the entire set of molecules for an svId
getSvMoleculesFromList <- function(sample, svId){
    reportProgress('getSvMoleculesFromList')
    chunk <- sampleData[[sample]]$svTable[
        SV_ID == svId,
        .(offset=CHUNK_OFFSET,size=CHUNK_SIZE)
    ]
    nodesFile <- getFindAllNodesFile(projectInfo, sample)    
    con = file(nodesFile, open = "rb")
    seek(con, where=chunk$offset, origin="start", rw="r")
    nodes <- readChar(con, chunk$size)
    close(con)    
    nodes <- fread(
        text=nodes,
        sep="\t",        
        header=FALSE,
        stringsAsFactors=FALSE,
        col.names = names(find$all_nodes),
        colClasses = unname(unlist(find$all_nodes))
    )
    nodes[,IS_REPEAT := as.logical(IS_REPEAT)]
    nodes
}

# P(X) = (X + 1) * (1/4)**X * (3/4)**2

