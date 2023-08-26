#----------------------------------------------------------------------
# handle unique junction loading and filtering
#----------------------------------------------------------------------

# populate the items list with all available amplicons from all loaded samples
svAmplicon_showAmpliconsDialog <- function(track, session, ...){
    showTrackItemsDialog(
        track$settings,
        session,
        title = "Select Amplicons",
        itemTypePlural = "Amplicons",
        tableData = reactive({
            uploadName <- appStepNamesByType$upload
            samples <- as.data.table(app[[uploadName]]$outcomes$samples())
            req(samples)
            amplicons <- assembleDataTable(function() samples, "amplicons", fread, CONSTANTS$ampliconsColNames) %>%
                         fillAmpliconDerivedColumns()
            amplicons[, .(
                ampliconKey = paste(sampleUniqueId, amplicon, sep = ":"),
                sampleName,
                chrom1, pos1, strand1,
                chrom2, pos2, strand2, 
                size
            )]
        }),
        keyColumn = "ampliconKey",
        extraColumns = c("sampleName","chrom1","pos1","strand1","chrom2","pos2","strand2","size"),
        size = "xl" # xl
    )
}

# app-specific junction loading\
# svAmplicon support one library per sample, but it might have pooled multiple discovered amplicons
getAmpliconSample <- function(ampliconKey){
    x <- as.list(strsplit(ampliconKey, ":")[[1]])
    names(x) <- c("Project","Sample_ID","amplicon")
    samples <- data.table(app$upload$outcomes$samples())
    c(x, list(sample = samples[Project == x$Project & Sample_ID == x$Sample_ID]))
}
svAmplicon_loadJunctions <- function(ampliconKey){ # Project_08673:Sample_Chr1_SH_2_IGO_08673_2:1
    x <- getAmpliconSample(ampliconKey)
    sourceId <- x$sample$Source_ID
    sampleUniqueId <- paste(x$Project, x$Sample_ID, sep = ":")
    jxns <- readRDS(getSourceFilePath(sourceId, "junctions")) 
    mts  <- readRDS(getSourceFilePath(sourceId, "moleculeTypes")) 
    jxns[, ":="(
        sourceId = sourceId,
        jxnKey = paste(ampliconKey, node1, node2, insertSize, insertBases, sep = ":"),
        node1 = pos1, # this may cause problem, not a true genomic node position
        node2 = pos2,
        cChrom1 = chrom1, # amplicons have de facto established the canonical strand
        cChrom2 = chrom2,        
        cRefPos1 = pos1,
        cRefPos2 = pos2,
        cStrand1 = ifelse(edgeType == "V" & strand1 == "+", 1, -1),
        cStrand2 = ifelse(edgeType == "V" & strand2 == "+", 1, -1),
        sampleUniqueId = sampleUniqueId,  
        samples = x$Sample_ID,
        nSamples = 1, # svAmplicon only has find, not multiSample compare
        nInstances = nMolTypes,
        nSequenced = nMolTypes,
        nLinkedJunctions = 0,
        flankLength = 0
    )]
    jxns <- merge(
        jxns[, .(sampleMolTypeKey = paste(sampleUniqueId, unlist(molTypeKeys), sep = ":")), by = .(jxnKey)],        
        jxns,
        by = "jxnKey",
        all.x = TRUE
    )
    pc <- CONSTANTS$groupedPathClasses
    mts[, sampleMolTypeKey := paste(sampleUniqueId, molKey, sep = ":")]
    mts[, pathClass := if(pathEdgeTypes %in% names(pc)) pc[[pathEdgeTypes]] else "other", by = "pathEdgeTypes"]  
    setkey(mts, sampleMolTypeKey)  
    jxns[, ":="(
        jxnUniqueKey = paste(sampleMolTypeKey, jxnKey, sep = ":"), # for plotting individual junction encounters
        pathClass = mts[sampleMolTypeKey, pathClass],
        molTypeKeys = NULL
    )]
    app$keepReject$getKeptJunctions(ampliconKey, jxns) # filter by keepReject app step
}
svAmplicon_filterJunctions <- function(jxns, track){ # filter by app-specific browser track settings
    Mol_Path_Classes <- getBrowserTrackSetting(track, "Filters", "Mol_Path_Classes", "deletion (D)") 
    Min_Reads <- getBrowserTrackSetting(track, "Filters", "Min_Reads", 1) 
    jxns[
        pathClass %in% Mol_Path_Classes & 
        nReadPairs >= Min_Reads
    ]
}

# convert sampleId to sampleName for legend
svAmplicon_legendSampleNames <- function(sampleIds, selectedTargets, ...){
    x <- do.call(rbind, lapply(names(selectedTargets), function(ampliconKey){
        x <- getAmpliconSample(ampliconKey)
        data.table(
            sampleId = x$sample$Sample_ID,
            sampleName = getSampleName(x$sample)
        )
    }))
    setkey(x, sampleId)
    x[gsub(",", "", sampleIds), sampleName]
}

# # get a single junction to expand
# svAmplicon_getJunction <- function(x){
#     svx_loadJunctions(x$sourceId, svAmplicon_loadJunctions)[clusterN == x$clusterN]
# }

# Classes 'data.table' and 'data.frame':  279 obs. of  22 variables:
#  $ ampliconId : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ node1      : chr  "chr1/+/33765060" "chr1/+/33764999" "chr1/+/33765021" "chr1
# /+/33764940" ...
#  $ node2      : chr  "chr1/+/33765090" "chr1/+/33765176" "chr1/+/33765103" "chr1
# /+/33765277" ...
#  $ insertSize : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ insertBases: chr  "*" "*" "*" "*" ...
#  $ nReadPairs : int  1 3 3 1 3 1 2 2 2 3 ...
#  $ nMolecules : int  1 1 3 1 1 1 2 1 2 2 ...
#  $ nMolTypes  : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ chrom1     : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ strand1    : chr  "+" "+" "+" "+" ...
#  $ pos1       : int  33765060 33764999 33765021 33764940 33765055 33764996 33765
# 001 33764818 33765022 33765001 ...
#  $ chrom2     : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ strand2    : chr  "+" "+" "+" "+" ...
#  $ pos2       : int  33765090 33765176 33765103 33765277 33765330 33765013 90139
# 353 33764834 33765139 33765174 ...
#  $ mapQ       : int  60 60 60 60 60 60 60 60 60 60 ...
#  $ alnQ       : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ baseQual   : int  37 38 38 38 38 35 36 37 37 38 ...
#  $ edgeType   : chr  "D" "D" "D" "D" ...
#  $ eventSize  : int  29 176 81 336 274 16 0 15 116 172 ...
#  $ jxnBases   : chr  NA "*" NA "*" ...
#  $ mergeLevel : int  2 2 2 2 2 2 0 2 2 2 ...
#  $ molTypeKeys:List of 279
#   ..$ : chr "1:1031503"
#   ..$ : chr "1:1038417"
