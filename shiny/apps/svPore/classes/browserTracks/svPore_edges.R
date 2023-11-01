#----------------------------------------------------------------------
# handle junction loading and filtering, for molecule plots, etc.
#----------------------------------------------------------------------
svPoreEdgeCreate <- "asNeeded"

# load all edges for a specific sourceId (not filtered yet)
svPore_loadSourceEdges <- function(targetId){
    req(targetId)
    startSpinner(session, message = "loading edges")
    sessionCache$get(
        'svPore_edges', 
        key = targetId, 
        permanent = TRUE, 
        from = "ram", 
        create = svPoreEdgeCreate, 
        createFn = function(...) {
            startSpinner(session, message = "loading edges .")
            edges <- readRDS(getSourceFilePath(targetId, "edgesFile")) 
            startSpinner(session, message = "loading edges ..")
            edges <- edges[, .SD, .SDcols = c(
                "sample","readI","channel",
                "clusterN",
                "duplex","duplex2","duplexCluster","duplexCluster2",
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
            startSpinner(session, message = "loading edges ...")
            setkey(edges, readKey, edgeN) 
            edges
        }
    )$value
}

# pull the edges for a specific junction
svPore_loadEdges <- function(targetId, clusterN_){
    sessionCache$get(
        'svPore_edges', 
        keyObject = list(targetId = targetId, clusterN = clusterN_), 
        permanent = FALSE, 
        from = "ram", 
        create = svPoreEdgeCreate, 
        createFn = function(...) {
            edges <- svPore_loadSourceEdges(targetId)
            startSpinner(session, message = "loading junction")
            readKeys <- edges[clusterN == clusterN_, unique(readKey)]
            edges[readKeys]
        }
    )$value
}


# set the genomic locus of all alignments; different loci plot in different plot "rows"
setAlignmentLoci <- function(edges, locusPadding){
    startSpinner(session, message = "setAlignmentLoci")
    maxQSize <- edges[edgeType == svx_edgeTypes$ALIGNMENT, max(qEnd)]
    getAlnRefPos <- function(readKey_, edgeN_) {
        if(edgeN_ == 1) edges[readKey_][edgeN_ + 1, refPos1] 
                   else edges[readKey_][edgeN_ - 1, refPos2] 
    }
    nLoci <- 0
    fillLocus <- function(edges){
        pending <- edges[, edgeType == svx_edgeTypes$ALIGNMENT & is.na(locus)]
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
    while(edges[edgeType == svx_edgeTypes$ALIGNMENT, any(is.na(locus))]) edges <- fillLocus(edges)
    edges
}

# retabulate reads for plotting
parseEdgeForMolPlot <- function(edges){
    startSpinner(session, message = "parseEdgeForMolPlot")
    maxI <- nrow(edges)
    dt <- do.call(rbind, lapply(1:maxI, function(i){
        if(edges[i, edgeType == svx_edgeTypes$ALIGNMENT]) data.table(
            readKey = edges[i, readKey],
            edgeN  = edges[i, edgeN],
            edgeType = svx_edgeTypes$ALIGNMENT,
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
    for(i in 1:maxI) if(dt[i, edgeType != svx_edgeTypes$ALIGNMENT]){
        ln2 <- dt[i + 1, locus1]
        dt[i, locus2 := ln2]
    }
    uniqueLoci <- dt[edgeType == svx_edgeTypes$ALIGNMENT, unique(locus1)]
    list(
        readKeys = edges[, unique(readKey)],
        loci = dt[edgeType == svx_edgeTypes$ALIGNMENT, locus1],
        uniqueLoci = uniqueLoci,
        nLoci = length(uniqueLoci),
        dt = dt,
        maxI = maxI
    )
}

# Classes ‘data.table’ and 'data.frame':  8803 obs. of  58 variables:
#  $ junctionKey          : chr  NA "25:25:1:1:16565:13:11" NA NA ...
#  $ sample               : chr  "HCT_P24_rapid_102223" "HCT_P24_rapid_102223" "HC
# T_P24_rapid_102223" "HCT_P24_rapid_102223" ...
#  $ readI                : int  3 3 3 8 8 8 9 9 9 11 ...
#  $ blockN               : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ edgeN                : int  1 2 3 1 2 3 1 2 3 1 ...
#  $ edgeType             : chr  "A" "U" "A" "A" ...
#  $ qName                : chr  "0bf747b3-f712-472b-8674-87208eed89a7" "0bf747b3-
# f712-472b-8674-87208eed89a7" "0bf747b3-f712-472b-8674-87208eed89a7" "ef7045f9-5f
# 47-43f2-8466-89b686db3b57" ...
#  $ node1                :integer64 3117289188 3117292066 3117275514 -3117279237 
# -3117275502 -3117292070 3117283096 3117292067 ...
#  $ qStart               : int  153 3014 3025 123 3851 3851 113 9055 9060 117 ...
#  $ node2                :integer64 3117292066 3117275514 3117289128 -3117275502 
# -3117292070 -3117286370 3117292067 3117275502 ...
#  $ qEnd                 : int  3014 3025 16608 3851 3851 9550 9055 9060 16614 16
# 84 ...
#  $ mapQ                 : int  60 60 60 60 60 60 60 60 60 60 ...
#  $ cigar                : chr  "4M1D43M2I502M1D50M1D35M1I56M1D83M1D3M1D8M1D17M1D
# 369M1D123M6D188M1I4M2D273M1D63M1D85M1D361M1I213M1I116M2D202M2D57M" NA "475M1D23M
# 1I20M1I12M1D7M1D8M2D84M1I69M3I4M1I3M1I76M1D102M1I70M1D60M1D1M1D1M2D9M1I273M1D219
# M1I33M1D140M2D204M1D33"| __truncated__ "232M1D73M1D726M1I3M1D293M1D870M1I95M1D42
# M1D3M1D266M1D107M1D499M1I94M1I1M2D173M1D3M1D233M1I10M" ...
#  $ gapCompressedIdentity: num  0.985 0.985 0.986 0.991 0.991 ...
#  $ eventSize            : int  2879 16552 13615 3736 16568 5701 8972 16565 7578 
# 1584 ...
#  $ insertSize           : int  NA 11 NA NA 0 NA NA 5 NA NA ...
#  $ foldback             : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#  $ jxnSeq               : chr  NA "ACAGTTTATAG" NA NA ...
#  $ baseQual             : num  21.2 7.2 28.9 24.8 NA 21.5 25.8 18.4 25.3 17.9 ..
# .
#  $ alnBaseQual          : num  NA 21.2 NA NA 21.5 ...
#  $ alnSize              : int  NA 2879 NA NA 3736 NA NA 7578 NA NA ...
#  $ channel              : int  426 426 426 567 567 567 669 669 669 792 ...      
#  $ pod5File             : chr  "PAK58770_pass_barcode01_1c6f0d53_8922c895_138.po
# d5" NA NA "PAK58770_pass_barcode01_1c6f0d53_8922c895_123.pod5" ...
#  $ duplex               : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ split                : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
#  $ passedBandwidth      : logi  NA TRUE NA NA TRUE NA ...
#  $ hasAdapter5          : logi  NA FALSE NA NA FALSE NA ...
#  $ hasAdapter3          : logi  NA FALSE NA NA FALSE NA ...
#  $ duplexCluster        : chr  "3" "3" "3" "8" ...
#  $ duplex2              : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ duplexCluster2       : chr  "3" "3" "3" "8" ...
#  $ matchable            : logi  NA TRUE NA NA TRUE NA ...
#  $ chrom1               : chr  "chrM" "chrM" "chrM" "chrM" ...
#  $ chromIndex1          : int  25 25 25 25 25 25 25 25 25 25 ...
#  $ refPos1              : int  13687 16565 13 3736 1 16569 7595 16566 1 1592 ...
#  $ strand1              : chr  "+" "+" "+" "-" ...
#  $ chrom2               : chr  "chrM" "chrM" "chrM" "chrM" ...
#  $ chromIndex2          : int  25 25 25 25 25 25 25 25 25 25 ...
#  $ refPos2              : int  16565 13 13627 1 16569 10869 16566 1 7578 9 ...  
#  $ strand2              : chr  "+" "+" "+" "-" ...
#  $ isCanonical          : logi  NA TRUE NA NA FALSE NA ...
#  $ cChromIndex1         : int  25 25 25 25 25 25 25 25 25 25 ...
#  $ cChromIndex2         : int  25 25 25 25 25 25 25 25 25 25 ...
#  $ cStrand1             : int  NA 1 NA NA 1 NA NA 1 NA NA ...
#  $ cStrand2             : int  NA 1 NA NA 1 NA NA 1 NA NA ...
#  $ cRefPos1             : int  NA 16565 NA NA 16569 NA NA 16566 NA NA ...       
#  $ cRefPos2             : int  NA 13 NA NA 1 NA NA 1 NA NA ...
#  $ clusterN             : num  NA 8 NA NA 8 NA NA 8 NA NA ...
#  $ isIndexJunctionKey   : logi  NA FALSE NA NA TRUE NA ...
#  $ nCanonical           : int  NA 1154 NA NA 1154 NA NA 1154 NA NA ...
#  $ nNonCanonical        : int  NA 1038 NA NA 1038 NA NA 1038 NA NA ...
#  $ icRefPos1            : int  NA 16569 NA NA 16569 NA NA 16569 NA NA ...       
#  $ icRefPos2            : int  NA 1 NA NA 1 NA NA 1 NA NA ...
#  $ iInsertSize          : int  NA 0 NA NA 0 NA NA 0 NA NA ...
#  $ clustered            : logi  NA TRUE NA NA TRUE NA ...
#  $ iRefPos1             : int  NA 16569 NA NA 1 NA NA 16569 NA NA ...
#  $ iRefPos2             : int  NA 1 NA NA 16569 NA NA 1 NA NA ...
#  $ segmentN             : num  1 1 1 1 1 1 1 1 1 1 ...
