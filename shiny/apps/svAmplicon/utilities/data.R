# support for svAmplicon data loading

# samples from a selected sampleSet
samplesReactive <- function(sampleSet) reactive({
    x <- sampleSet$assignments()
    req(x, is.data.table(x), nrow(x) > 0)   
    x 
})

# general function for loading a data.table from the pipeline output
# prepend columns sourceId and sample to the pipeline data.table
assembleDataTable <- function(samples, fileType, loadFn, colNames = NULL){
    x <- samples()
    do.call(rbind, lapply(seq_len(nrow(x)), function(i){
        print(getSourceFilePath(x$Source_ID[i], fileType))
        dt <- loadFn(getSourceFilePath(x$Source_ID[i], fileType))
        if(!is.null(colNames)) setnames(dt, colNames)
        cbind(
            sourceId = x$Source_ID[i],
            sample = getSampleNames(sampleUniqueIds = paste(x$Project[i], x$Sample_ID[i], sep = ":")),
            dt
        )
    }))
}

# amplicons; could be more than one for a given sample if a mixed library
getAmpliconKeys <- function(dt){
    paste(dt$sample, dt$amplicon, sep = ":")
}
ampliconsReactive <- function(samples) reactive({
    x <- assembleDataTable(samples, "amplicons", fread, colNames = c(
        "amplicon","type","nReadPairs",
        "chrom1","side1","pos1","ref1","primer1",
        "chrom2","side2","pos2","ref2","primer2"
    ))
    x[, size := pos2 - pos1 + 1]
})
ampliconsTableServer <- function(parentId, input, amplicons, selection = "single") bufferedTableServer(
    "ampliconsTable",
    parentId,
    input,
    tableData = reactive({
        amplicons()[, .SD, .SDcols = c(
            "sample",
            "amplicon","type",#"nReadPairs",
            "chrom1","side1","pos1",
            "chrom2","side2","pos2",
            "size"
        )]
    }),
    selection = selection
)
selectedAmpliconsReactive <- function(amplicons, ampliconsTable) reactive({
    I <- ampliconsTable$rows_selected()
    req(I)
    amplicons()[I]    
})

# moleculeTypes, i.e., one row for each distinct path through an amplicon molecule
moleculeTypesReactive <- function(samples) reactive({ assembleDataTable(samples, "moleculeTypes", readRDS) })
ampliconMoleculeTypesReactive <- function(moleculeTypes, selectedAmplicons) reactive({
    moleculeTypes <- moleculeTypes()
    moleculeTypes[
        getAmpliconKeys(moleculeTypes) %in% getAmpliconKeys(selectedAmplicons())
    ]
})

# moleculeType, i.e., metadata on a single selected moleculeType
moleculeMetadataReactive <- function(amplicons, moleculeTypes) reactive({
    amplicons <- amplicons()
    moleculeTypes <- moleculeTypes()
    req(amplicons, moleculeTypes)

    # TODO: create molecule list scroll

    dprint(nrow(moleculeTypes))

    moleculeType <- moleculeTypes[sample(nrow(moleculeTypes), 1)]
    # moleculeType <- moleculeTypes[molId == 2473611] # translocation MAPQ = 0
    # moleculeType <- moleculeTypes[molId == 3602861] # translocation MAPQ = 60, 122
    # moleculeType <- moleculeTypes[molId == 813957]  # large deletion, split
    # moleculeType <- moleculeTypes[molId == 2972341] # deletion -3 microhomology
    # moleculeType <- moleculeTypes[molId == 1776032] # internal rearrangement
#  moleculeType <- moleculeTypes[molId == 3470234]
    # moleculeType <- moleculeTypes[molId == 596477]


    

    # molecule-level information
    amplicon <- amplicons[ # currently only works for expectOverlap
        getAmpliconKeys(amplicons) == getAmpliconKeys(moleculeType)
    ]
    alnReadNs <- as.integer(strsplit(moleculeType$readNs, ":")[[1]])
    jxnReadNs <- integer()
    for(i in 1:(length(alnReadNs) - 1)) if(alnReadNs[i] == alnReadNs[i + 1]) jxnReadNs <- c(jxnReadNs, alnReadNs[i])
    moleculeType$tLen1 <- nchar(moleculeType$seq1)
    moleculeType$tLen2 <- if(moleculeType$mergeLevel > 0) NA else nchar(moleculeType$seq2)
    moleculeType$alnReadNs <- paste(alnReadNs, collapse = ":")
    moleculeType$jxnReadNs <- paste(jxnReadNs, collapse = ":")

    # out-of-amplicon segments
    chroms <- strsplit(moleculeType$chroms, ":")[[1]]
    poss   <- as.integer(strsplit(moleculeType$poss, ":")[[1]])
    ampliconic <- chroms == amplicon$chrom1 & poss >= amplicon$pos1 & poss <= amplicon$pos2

    list(
        moleculeType = moleculeType,
        amplicon = amplicon,
        ampliconSize = amplicon$pos2 - amplicon$pos1 + 1,
        tLen = moleculeType[,
            if(mergeLevel == 0) nchar(seq1) + nchar(seq2) - (if(is.na(overlap)) 0 else overlap) # sometimes overlap is NA
            else nchar(seq1)
        ],
        isMerged = moleculeType$mergeLevel > 0,
        ampliconic = ampliconic,

        # alignment-level information
        chroms = chroms,
        strands = strsplit(moleculeType$strands, ":")[[1]],
        poss   = poss,
        cigars = strsplit(moleculeType$cigars, ":")[[1]],
        mapQs = strsplit(moleculeType$mapQs, ":")[[1]],
        baseQuals = strsplit(moleculeType$baseQuals, ":")[[1]],
        alnReadNs = alnReadNs,
        jxnReadNs = jxnReadNs,
        sizes  = as.integer(strsplit(moleculeType$sizes, ":")[[1]]),
        insSizes = as.integer(strsplit(moleculeType$insSizes, ":")[[1]]),
        types = strsplit(moleculeType$pathClass, "")[[1]]
    )
})
moleculeMetadataUI <- function(moleculeMetadata) renderUI({
    mmd <- moleculeMetadata()
    req(mmd)
    mt <- mmd$moleculeType
    keyValuePair <- function(key) {
        tags$div(
            tags$div(paste(key, ":"), style = "display: inline-block; width: 100px; text-align: right; margin-right: 5px;"),
            tags$div(gsub(":", " ", mt[[key]]), style = "display: inline-block; margin-left:  5px;")
        )
    }
    tags$div(
        lapply(c(
            "ampliconId","molId",
            "nReadPairs","nMols",
            "mergeLevel","overlap","tLen1","tLen2",
            "pathClass",
            "alnReadNs",
            "strands","chroms","poss","cigars","mapQs","baseQuals",
            "jxnReadNs",
            "sizes","insSizes"
        ), keyValuePair)
    )
})
# List of 11
#  $ moleculeType:Classes ‘data.table’ and 'data.frame':  1 obs. of  24 variables:
#   ..$ sourceId  : chr "8d81de94383d6a8a55cdb4181a410f69"
#   ..$ sample    : chr "Sample_Chr1_SH_Mre11_1_IGO_08673_E_1"
#   ..$ ampliconId: int 1
#   ..$ path      : chr "chr1/+/33765030:chr1/+/33765125"
#   ..$ insSizes  : chr "0"
#   ..$ molKey    : chr "1:1402597"
#   ..$ molId     : int 1402597
#   ..$ isRef     : logi FALSE
#   ..$ nReadPairs: int 9
#   ..$ nMols     : int 4
#   ..$ mergeLevel: int 2
#   ..$ overlap   : int 121
#   ..$ chroms    : chr "chr1:chr1"
#   ..$ poss      : chr "33764688:33765125"
#   ..$ cigars    : chr "343M:138M"
#   ..$ readNs    : chr "1:1"
#   ..$ mapQs     : chr "60:60"
#   ..$ pathClass : chr "D"
#   ..$ nodePairs :List of 1
#   .. ..$ : chr "chr1/+/33765030:chr1/+/33765125"
#   ..$ sizes     : chr "94"
#   ..$ seq1      : chr "TCATGTCCGGTTTGGATGCCAAACGACCCTTTAAAAGGGGTTGCCTAAGAGCATCAG
# AAAGCAGATATTTATATTTTGATTCACAAGAGTAGCAAAATTACAGCTAAGAAG"| __truncated__
#   ..$ seq2      : chr "*"
#   ..$ qual1     : chr "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"| __truncated__
#   ..$ qual2     : chr "*"
#   ..- attr(*, ".internal.selfref")=<externalptr>
#  $ amplicon    :Classes ‘data.table’ and 'data.frame':  1 obs. of  16 variables:
#   ..$ sourceId  : chr "8d81de94383d6a8a55cdb4181a410f69"
#   ..$ sample    : chr "Sample_Chr1_SH_Mre11_1_IGO_08673_E_1"
#   ..$ amplicon  : int 1
#   ..$ type      : chr "expectOverlap"
#   ..$ nReadPairs: int 3179461
#   ..$ chrom1    : chr "chr1"
#   ..$ side1     : chr "R"
#   ..$ pos1      : int 33764688
#   ..$ ref1      : chr "TCATGTCCGGTTTGGATGCCAAACGACCCTTTAAAAGGGGTTGCCTAAGAGCATCAG
# AAAGCAGATATTTATATTTTGATTCACAAGAGTAGCAAAATTACAGCTAAGAAG"| __truncated__
#   ..$ primer1   : chr "TCATGTCCGGTTTGGATGCCAAACG"
#   ..$ chrom2    : chr "chr1"
#   ..$ side2     : chr "L"
#   ..$ pos2      : int 33765262
#   ..$ ref2      : chr "*"
#   ..$ primer2   : chr "ATCCTCTTCTGGAGTCTTTGATTCG"
#   ..$ size      : num 575
#   ..- attr(*, ".internal.selfref")=<externalptr>
#  $ ampliconSize: num 575
#  $ tLen        : int 481
#  $ isMerged    : logi TRUE
#  $ chroms      : chr [1:2] "chr1" "chr1"
#  $ poss        : int [1:2] 33764688 33765125
#  $ cigars      : chr [1:2] "343M" "138M"
#  $ readNs      : int [1:2] 1 1
#  $ sizes       : int 94
#  $ insSizes    : int 0

# Classes ‘data.table’ and 'data.frame':  234 obs. of  22 variables:
#  $ sourceId  : chr  "15dcd64821b17e6d392c2a58dc404ef5" "15dcd64821b17e6d392c2a58
# dc404ef5" "15dcd64821b17e6d392c2a58dc404ef5" "15dcd64821b17e6d392c2a58dc404ef5" 
# ...
#  $ sample    : chr  "Sample_Chr1_SH_Mre11_1_IGO_08673_E_1" "Sample_Chr1_SH_Mre11
# _1_IGO_08673_E_1" "Sample_Chr1_SH_Mre11_1_IGO_08673_E_1" "Sample_Chr1_SH_Mre11_1
# _IGO_08673_E_1" ...
#  $ ampliconId: int  1 1 1 1 1 1 1 1 1 1 ...
#  $ path      : chr  "" "chr1/+/33765001:chr1/+/33764688" "chr1/+/33765021:chr1/+
# /33764688" "chr1/+/33764999:chr1/+/33764688" ...
#  $ insSizes  : chr  "" "0" "0" "0" ...
#  $ molKey    : chr  "1:2948997" "1:3293777" "1:1369849" "1:3400260" ...
#  $ molId     : int  2948997 3293777 1369849 3400260 1402597 687422 1839587 16099
# 30 1254465 2716327 ...
#  $ isRef     : logi  TRUE FALSE FALSE FALSE FALSE FALSE ...
#  $ nReadPairs: int  3447982 39 24 17 12 10 10 9 9 9 ...
#  $ nMols     : int  1019895 32 21 11 7 8 9 4 9 5 ...
#  $ mergeLevel: int  2 2 2 2 2 2 2 2 2 2 ...
#  $ overlap   : int  27 79 68 125 121 72 46 278 102 97 ...
#  $ chroms
#  $ poss
#  $ cigars    : chr  "575M" "314M:209M" "334M:200M" "312M:165M" ...
#  $ readNs    : chr  "1" "1:1" "1:1" "1:1" ...
#  $ mapQs     : chr  "60" "60:60" "60:60" "60:60" ...
#  $ pathClass : chr  "" "D" "D" "D" ...
#  $ nodePairs :List of 234
#   ..$ : chr
#   ..$ : chr "chr1/+/33765001:chr1/+/33764688"
#   .. [list output truncated]
#  $ sizes     : chr  "" "52" "41" "98" ...
#  $ seq1      : chr  "TCATGTCCGGTTTGGATGCCAAACGACCCTTTAAAAGGGGTTGCCTAAGAGCATCAGAA
# AGCAGATATTTATATTTTGATTCACAAGAGTAGCAAAATTACAGCTAAGAAG"| __truncated__ "TCATGTCCGG
# TTTGGATGCCAAACGACCCTTTAAAAGGGGTTGCCTAAGAGCATCAGAAAGCAGATATTTATATTTTGATTCACAAGAGT
# AGCAAAATTACAGCTAAGAAG"| __truncated__ "TCATGTCCGGTTTGGATGCCAAACGACCCTTTAAAAGGGGT
# TGCCTAAGAGCATCAGAAAGCAGATATTTATATTTTGATTCACAAGAGTAGCAAAATTACAGCTAAGAAG"| __trunc
# ated__ "TCATGTCCGGTTTGGATGCCAAACGACCCTTTAAAAGGGGTTGCCTAAGAGCATCAGAAAGCAGATATTTAT
# ATTTTGATTCACAAGAGTAGCAAAATTACAGCTAAGAAG"| __truncated__ ...
#  $ seq2      : chr  "*" "*" "*" "*" ...
#  $ qual1     : chr  "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"| __truncated__ "CCCCCCCCCC
# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
# CCCCCCCCCCCCCCCCCCCCC"| __truncated__ "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"| __trunc
# ated__ "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"| __truncated__ ...
#  $ qual2     : chr  "*" "*" "*" "*" ...
#  - attr(*, ".internal.selfref")=<externalptr>

# Classes ‘data.table’ and 'data.frame':  1 obs. of  16 variables:
#  $ sourceId  : chr "15dcd64821b17e6d392c2a58dc404ef5"
#  $ sample    : chr "Sample_Chr1_SH_Mre11_1_IGO_08673_E_1"
#  $ amplicon  : int 1
#  $ type      : chr "expectOverlap"
#  $ nReadPairs: int 3179461
#  $ chrom1    : chr "chr1"
#  $ side1     : chr "R"
#  $ pos1      : int 33764688
#  $ ref1      : chr "TCATGTCCGGTTTGGATGCCAAACGACCCTTTAAAAGGGGTTGCCTAAGAGCATCAGAAA
# GCAGATATTTATATTTTGATTCACAAGAGTAGCAAAATTACAGCTAAGAAG"| __truncated__
#  $ primer1   : chr "TCATGTCCGGTTTGGATGCCAAACG"
#  $ chrom2    : chr "chr1"
#  $ side2     : chr "L"
#  $ pos2      : int 33765262
#  $ ref2      : chr "*"
#  $ primer2   : chr "ATCCTCTTCTGGAGTCTTTGATTCG"
#  $ size      : num 575
#  - attr(*, ".internal.selfref")=<externalptr>



# filteredMoleculeTypes <- reactive({
#     moleculeTypes <- moleculeTypes()
#     moleculeTypes[
#         getAmpliconKeys(moleculeTypes) %in% getAmpliconKeys(selectedAmplicons()) & 
#         mergeLevel > 0 & 
#         grepl(CONSTANTS$svTypes$Deletion, jxnsKey, ignore.case = TRUE) # TODO: deploy dynamic SV type filtering
#     ]
# })
# tableFilteredMoleculeTypes <- reactive({
#     junctionI <- junctionsTable$rows_selected()
#     if(!isTruthy(junctionI)) return(filteredMoleculeTypes())
#     jxnKey <- filteredJunctions()[junctionI, jxnKey]
#     filteredMoleculeTypes()[sapply(jxnKeys, function(x) jxnKey %in% x)]
# })
# junctions <- reactive({
#     assembleDataTable("junctions", readRDS)
# })
# filteredJunctions <- reactive({
#     junctions <- junctions()
#     junctions[
#         getAmpliconKeys(junctions) %in% getAmpliconKeys(selectedAmplicons()) &
#         hasMerged == TRUE &
#         svType == CONSTANTS$svTypes$Deletion &  # TODO: deploy dynamic SV type filtering
#         TRUE
#     ]
# })

# moleculeTypesTable <- bufferedTableServer(
#     "moleculeTypesTable",
#     id,
#     input,
#     tableData = reactive({
#         tableFilteredMoleculeTypes()[, .SD, .SDcols = c(
#             "sample",
#             "amplicon","molTypeId",
#             "nMols","nReadPairs",        
#             "isRef","mergeLevel","overlap","tLen",
#             "molClass","jxnsKey","avgQual","maxIntMapQ",
#             "nAlns","nJxns"
#         )]
#     }),
#     selection = 'single'
# )
# junctionsTable <- bufferedTableServer(
#     "junctionsTable",
#     id,
#     input,
#     tableData = reactive({
#         filteredJunctions()[, .SD, .SDcols = c(
#             "sample",
#             "amplicon","nMolTypes","nMols","nReadPairs",
#             "hasMerged","tLens","overlap","jxnBases","callTypes",
#             "svType","svSize",
#             "chrom1","side1","pos1",
#             "chrom2","side2","pos2"
#         )]
#     }),
#     selection = 'single'
# )
