#----------------------------------------------------------------------
# handle junction loading and filtering, for molecule plots, etc.
#----------------------------------------------------------------------

# load all edges for a specific sourceId (not filtered yet)
loadSourceEdges <- function(sourceId){
    req(sourceId)
    startSpinner(session, message = "loading edges")
    svPoreCache$get(
        'svEdges', 
        key = sourceId, 
        permanent = TRUE, 
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            startSpinner(session, message = "loading edges from disk")
            edges <- readRDS(getSourceFilePath(sourceId, "edgesFile")) 
#  [1] "junctionKey"           "sample"                "readI"
#  [4] "blockN"                "edgeN"                 "qName"
#  [7] "node1"                 "qStart"                "node2"
# [10] "qEnd"                  "mapQ"                  "cigar"
# [13] "gapCompressedIdentity" "edgeType"              "eventSize"            
# [16] "insertSize"            "nStrands"              "baseQual"
# [19] "alnBaseQual"           "alnSize"               "sStart"
# [22] "sEnd"                  "channel"               "pod5File"
# [25] "passedBandwidth"       "hasAdapter5"           "hasAdapter3"
# [28] "matchable"             "chrom1"                "chromIndex1"
# [31] "refPos1"               "strand1"               "chrom2"
# [34] "chromIndex2"           "refPos2"               "strand2"
# [37] "isCanonical"           "cChromIndex1"          "cChromIndex2"
# [40] "cStrand1"              "cStrand2"              "cRefPos1"
# [43] "cRefPos2"              "clusterN"              "isIndexJunctionKey"       
# [46] "nCanonical"            "nNonCanonical"         "icRefPos1"
# [49] "icRefPos2"             "iInsertSize"           "clustered"
# [52] "iRefPos1"              "iRefPos2"              "segmentN"
            edges <- edges[, .SD, .SDcols = c(
                "sample","readI","channel",
                "clusterN","nStrands",
                "edgeN","edgeType",
                "qStart","qEnd",
                "node1","node2",
                "chrom1","refPos1","strand1",
                "chrom2","refPos2","strand2",
                "cigar","eventSize","insertSize"
            )]            
            edges[, ":="(
                readKey = paste(sample, readI, sep = ":")
            )]
            setkey(edges, readKey, edgeN) 
            edges
        }
    )$value
}

# pull the edges for a specific SV cluster
loadClusterEdges <- function(sourceId, clusterN_){
    svPoreCache$get(
        'svEdges', 
        keyObject = list(sourceId = sourceId, clusterN = clusterN_), 
        permanent = FALSE, 
        from = "ram", 
        create = "asNeeded", 
        createFn = function(...) {
            edges <- loadSourceEdges(sourceId)
            startSpinner(session, message = "loading cluster")
            readKeys <- edges[clusterN == clusterN_, unique(readKey)]
            edges[readKeys]
        }
    )$value
}

# set the genomic locus of all alignments; different loci plot in different plot "rows"
setAlignmentLoci <- function(edges, locusPadding){
    startSpinner(session, message = "setAlignmentLoci")
    maxQSize <- edges[edgeType == edgeTypes$ALIGNMENT, max(qEnd)]
    getAlnRefPos <- function(readKey_, edgeN_) {
        if(edgeN_ == 1) edges[readKey_][edgeN_ + 1, refPos1] 
                   else edges[readKey_][edgeN_ - 1, refPos2] 
    }
    nLoci <- 0
    fillLocus <- function(edges){
        pending <- edges[, edgeType == edgeTypes$ALIGNMENT & is.na(locus)]
        i <- edges[, min(which(pending), na.rm = TRUE)]
        chrom <- edges[i, chrom1]  
        refPos  <- edges[i, getAlnRefPos(readKey, edgeN)]
        locus_ <- abs(edges[i, node2])
        rRange <- c(
            refPos - locusPadding - maxQSize,
            refPos + locusPadding + maxQSize 
        )
        edges[pending, locus := {
            if(chrom1 == chrom && getAlnRefPos(readKey, edgeN) %between% rRange) locus_
            else NA_integer64_
        }, by = .(readKey, edgeN)]
        nLoci <<- nLoci + 1
        edges
    }
    edges[, locus := NA_integer64_]    
    while(edges[edgeType == edgeTypes$ALIGNMENT, any(is.na(locus))]) edges <- fillLocus(edges)
    edges
}

# retabulate reads for plotting
parseEdgeForMolPlot <- function(edges){
    startSpinner(session, message = "parseEdgeForMolPlot")
    maxI <- nrow(edges)
    dt <- do.call(rbind, lapply(1:maxI, function(i){
        if(edges[i, edgeType == edgeTypes$ALIGNMENT]) data.table(
            readKey = edges[i, readKey],
            edgeN  = edges[i, edgeN],
            edgeType = edgeTypes$ALIGNMENT,
            cigar   = edges[i, cigar],
            locus1  = edges[i, locus],
            locus2  = edges[i, locus],
            chrom   = edges[i, chrom1],
            strand  = edges[i, if(strand1 == "+") 1 else -1],
            qStart  = edges[i, qStart],
            qEnd    = edges[i, qEnd],
            rStart  = edges[i, refPos1],
            rEnd    = edges[i, refPos2],
            qOffset = NA_integer_,
            rOffset = NA_integer_,
            label   = edges[i, commify(eventSize)]
        ) else data.table(
            readKey = edges[i, readKey],
            edgeN  = edges[i, edgeN],
            edgeType = edges[i, edgeType],
            cigar    = NA_character_,
            locus1  = edges[i - 1, locus],
            locus2  = NA_integer64_,
            chrom   = NA_character_,
            strand  = NA_integer_,
            qStart  = edges[i - 1, qEnd],
            qEnd    = edges[i + 1, qStart],
            rStart  = edges[i, refPos1],
            rEnd    = edges[i, refPos2],
            qOffset = edges[i, insertSize],
            rOffset = edges[i, eventSize],
            label   = NA_character_
        )
    }))
    for(i in 1:maxI) if(dt[i, edgeType != edgeTypes$ALIGNMENT]){
        ln2 <- dt[i + 1, locus1]
        dt[i, locus2 := ln2]
    }
    uniqueLoci <- dt[edgeType == edgeTypes$ALIGNMENT, unique(locus1)]
    list(
        readKeys = edges[, unique(readKey)],
        loci = dt[edgeType == edgeTypes$ALIGNMENT, locus1],
        uniqueLoci = uniqueLoci,
        nLoci = length(uniqueLoci),
        dt = dt,
        maxI = maxI
    )
}
