# examine all called TT and TA class SVs for the presence of concomitant SNVs/indels

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(data.table)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'DATA_NAME',
        'ACTION_DIR',
        'MODULES_DIR',
        'FIND_PREFIX',
        'GENOTYPE_PREFIX',
        'GENOME',
        'GENOME_CHROMS',
        'GENOME_FASTA',
        'CHROM_FASTA_DIR'        
    ),
    integer = c(
        'N_CPU',
        'REGION_PADDING'
    )
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) ########################
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'),   c('general', 'IUPAC', 'smith_waterman'))
# sourceScripts(file.path(rUtilDir, 'genome'), c('general', 'chroms', 'targets', 'faidx'))
#-------------------------------------------------------------------------------------
nodeClasses <- list(
    GAP         = 0, # SV evidence type codes, i.e., node classes
    SPLIT       = 1,
    OUTER_CLIP  = 2
)
clipBases <- c('a', 'c', 'g', 't', 'n')
unpackNodeNames <- function(nodeNames){ 
    x <- as.data.frame(t(sapply(nodeNames, function(x){
        if(x == "*") c("0", "*", "0") else strsplit(x, '\\W')[[1]]
    })), stringsAsFactors = FALSE)
    x[[1]] <- as.integer(x[[1]])
    x[[3]] <- as.integer(x[[3]])
    x
}
#=====================================================================================

#=====================================================================================
# load required data
#-------------------------------------------------------------------------------------
message("loading unphased haplotype map")
hapFile <- paste(env$GENOTYPE_PREFIX, 'unphased_haplotypes', 'gz', sep = ".")
hapMap <- fread(
    hapFile,
    header = TRUE,
    sep = "\t"
)

message("loading structural variants")
svsFile  <- paste(env$FIND_PREFIX, 'structural_variants', 'rds', sep = ".")
molsFile <- paste(env$FIND_PREFIX, 'junction_molecules', 'rds', sep = ".")
svCalls <- readRDS(svsFile)
jxnMols <- readRDS(molsFile)

message("filtering structural variants to sequenced TT and TA target classes")
svCalls <- svCalls[TARGET_CLASS %in% c("TT", "TA") & JXN_BASES != "*"]
jxnMols[, isOuterClip := NODE_CLASS == nodeClasses$OUTER_CLIP]
jxnMols[, c('chrom1', 'side1', 'pos1') := unpackNodeNames(NODE_1)]
jxnMols[, c('chrom2', 'side2', 'pos2') := unpackNodeNames(NODE_2)]
#=====================================================================================

#=====================================================================================
# convert an SV and its set of evidence molecules to a base map
# ----------------------------------------------------------------------

getJunctionMap <- function(sv, mols){

    # initialize map and dimensions
    refWidth <- nchar(sv$GEN_REF_1)
    faidxPadding <- (refWidth - 1) / 2
    microhomologyLength <- sv[, MICROHOM_LEN]    
    leftRefI  <- faidxPadding + 1
    rightRefI <- leftRefI + 1 - microhomologyLength
    mapWidth <- refWidth + 1 - microhomologyLength
    nAln <- mols[, .N + sum(!isOuterClip)]
    yOffset <- 1
    baseMap <- matrix(NA, nrow = mapWidth, ncol = nAln) # for consensus making

    # parse the genome reference lines
    GEN_REF_1 <- sv$GEN_REF_1
    GEN_REF_2 <- sv$GEN_REF_2
    if(sv$SIDE_1 == sv$SIDE_2){ # flip one reference strand for inversion
        if(sv$SIDE_1 == "L") GEN_REF_2 <- tolower(rc(GEN_REF_2))
                        else GEN_REF_1 <- tolower(rc(GEN_REF_1))
    }
    GEN_REF_1 <- c(    strsplit(GEN_REF_1, "")[[1]],  rep(" ", mapWidth))
    GEN_REF_2 <- c(rev(strsplit(GEN_REF_2, "")[[1]]), rep(" ", mapWidth))
    GEN_REF_1 <- GEN_REF_1[1:mapWidth]
    GEN_REF_2 <- GEN_REF_2[1:mapWidth]
    GEN_REF_2 <- rev(GEN_REF_2)

    # adjust actual node positions based on microhomology / insertion
    if(microhomologyLength > 0){ 
        left  <- leftRefI - microhomologyLength + 1
        right <- leftRefI
    } else if(microhomologyLength < 0){
        left  <- leftRefI + 1
        right <- leftRefI - microhomologyLength
    }    
    
    # add all alignments to map
    processNode <- function(nodeN, isRef, jxnSide, jxnPos, CIGAR, SEQ, pos){
        if(CIGAR == "*") return(NULL) # the missing node of an outer clip
        aln <- cigarToRefAln(CIGAR, SEQ)           
        nodePosDelta <- if(jxnSide == "L") jxnPos - pos else pos - jxnPos
        mapPos <- if(nodeN == 1){
            leftRefI - nodePosDelta - aln$lengthAlignedOut - aln$leftClip + 1
        } else {
            leftRefI + 1 - microhomologyLength + nodePosDelta - aln$leftClip
        }
        j <- mapPos:(mapPos + aln$lengthOut - 1)
        baseMap[j, yOffset] <<- aln$seq                
        yOffset <<- yOffset + 1 
    }
    mols[
        order( CLIP_LEN_1, POS_1), 
        mapply(processNode, 1, IS_REFERENCE, sv$SIDE_1, sv$POS_1, CIGAR_1, SEQ_1, pos1)
    ]
    mols[
        order(-CLIP_LEN_2, POS_2 + nchar(SEQ_2)), 
        mapply(processNode, 2, IS_REFERENCE, sv$SIDE_2, sv$POS_2, CIGAR_2, SEQ_2, pos2)
    ]

    # assemble the final maps
    usedPos <- rle(apply(baseMap[, ], 1, function(v) any(!is.na(v))))$lengths
    usedPos <- (usedPos[1] + 1):(mapWidth - usedPos[length(usedPos)]) # row indices = base positions in SV allele

    # finish and return the results
    list(
        text = baseMap, 
        leftRefI = leftRefI,
        rightRefI = rightRefI,
        usedPos = usedPos,
        GEN_REF_1 = GEN_REF_1,
        GEN_REF_2 = GEN_REF_2
    )
}
# ----------------------------------------------------------------------
# create an alignment betweeen junction molecules, their consensus, and the genome reference
# ----------------------------------------------------------------------
getJunctionConsensus <- function(text){
    apply(text, 1, function(x){
        x <- x[!is.na(x)]
        if(length(x) == 0) return("~")
        xx <- x[!(x %in% clipBases)] # insertions are also lower case
        if(length(xx) == 0) xx <- x
        agg <- aggregate(xx, list(xx), length)
        agg[which.max(agg[[2]]), 1]
    })
}
#=====================================================================================

#=====================================================================================
# load required data
#-------------------------------------------------------------------------------------
x <- do.call(rbind, lapply(svCalls$SV_ID, function(svId){
    sv <- svCalls[SV_ID == svId]
    mols <- jxnMols[SV_ID == svId]
    jxnMap <- getJunctionMap(sv, mols)
    consensus <- getJunctionConsensus(jxnMap$text)

    str(sv)
    str(mols)
    str(jxnMap)
    str(consensus)

    stop(1111)



}))
#=====================================================================================

# #=====================================================================================
# # xxxxxxxxxxxxxxxxxxxxxxxx
# #-------------------------------------------------------------------------------------
# sv_del_dup <- function(hf1_file, startPos, index){
#     row_map_1 <- svJunctions[[index]]$assembly$consensusMap$refPos1 - startPos + 1
#     gene_map1 <- hf1_file[row_map_1,]
#     row_map_2 <- svJunctions[[index]]$assembly$consensusMap$refPos2 - startPos + 1
#     gene_map2 <- hf1_file[row_map_2,]

#     node1_hap1All <- gene_map1$hap1
#     node1_hap2All <- gene_map1$hap2
#     node2_hap1All <- gene_map2$hap1
#     node2_hap2All <- gene_map2$hap2

#     node1_row <- match(svJunctions[[index]]$sv$POS_1, svJunctions[[index]]$assembly$consensusMap$junctionPos1)
#     node2_row <- match(svJunctions[[index]]$sv$POS_2, svJunctions[[index]]$assembly$consensusMap$junctionPos2)

#     junction_length <- length(svJunctions[[index]]$assembly$consensusMap$refPos1)
#     node1_hap1Masked <- replace(node1_hap1All, seq(node1_row +1, junction_length), NA)
#     node1_hap2Masked <- replace(node1_hap2All, seq(node1_row +1, junction_length), NA)
#     node2_hap1Masked <- replace(node2_hap1All, seq(1, node2_row-1), NA)
#     node2_hap2Masked <- replace(node2_hap2All, seq(1, node2_row-1), NA)

#     new_consensus_map <- data.table(
#     refPos1 = svJunctions[[index]]$assembly$consensusMap$refPos1,
#     refPos2 = svJunctions[[index]]$assembly$consensusMap$refPos2,
#     junctionPos1 = svJunctions[[index]]$assembly$consensusMap$junctionPos1,
#     junctionPos2 = svJunctions[[index]]$assembly$consensusMap$junctionPos2,
#     nMols = svJunctions[[index]]$assembly$consensusMap$nMols,
#     consensus = svJunctions[[index]]$assembly$consensusMap$consensus,
#     node1_refAll = svJunctions[[index]]$assembly$consensusMap$ref1All,
#     node1_refMasked = svJunctions[[index]]$assembly$consensusMap$ref1Masked,
#     node2_refAll = svJunctions[[index]]$assembly$consensusMap$ref2All,
#     node2_refMasked = svJunctions[[index]]$assembly$consensusMap$ref2Masked,
#     node1_hap1All = node1_hap1All,
#     node1_hap1Masked = node1_hap1Masked,
#     node1_hap2All = node1_hap2All,
#     node1_hap2Masked = node1_hap2Masked,
#     node2_hap1All = node2_hap1All,
#     node2_hap1Masked = node2_hap1Masked,
#     node2_hap2All = node2_hap2All,
#     node2_hap2Masked = node2_hap2Masked,
#     isRef1 = svJunctions[[index]]$assembly$consensusMap$consensus == svJunctions[[index]]$assembly$consensusMap$ref1Masked,
#     isRef2 = svJunctions[[index]]$assembly$consensusMap$consensus == svJunctions[[index]]$assembly$consensusMap$ref2Masked,
#     isNode1Hap1 = svJunctions[[index]]$assembly$consensusMap$consensus == node1_hap1Masked,
#     isNode1Hap2 = svJunctions[[index]]$assembly$consensusMap$consensus == node1_hap2Masked,
#     isNode2Hap1 = svJunctions[[index]]$assembly$consensusMap$consensus == node2_hap1Masked,
#     isNode2Hap2 =svJunctions[[index]]$assembly$consensusMap$consensus == node2_hap2Masked
#     )
#     return(new_consensus_map)
# }


# #=====================================================================================



# sv_inversion <- function(hf1_file, startPos, index){
#   row_map_1 <- svJunctions[[index]]$assembly$consensusMap$refPos1 - startPos + 1
#   gene_map1 <- hf1_file[row_map_1,]
#   row_map_2 <- svJunctions[[index]]$assembly$consensusMap$refPos2 - startPos + 1
#   gene_map2 <- hf1_file[row_map_2,]
#   junction_length <- length(svJunctions[[index]]$assembly$consensusMap$refPos1)
  
#   if(gene_map1[1,1] < gene_map1[2,1]){
#     node1_hap1All <- gene_map1$hap1
#     node1_hap2All <- gene_map1$hap2
#     node1_row <- match(svJunctions[[index]]$sv$POS_1, svJunctions[[index]]$assembly$consensusMap$junctionPos1)
#     node1_hap1Masked <- replace(node1_hap1All, seq(node1_row +1, junction_length), NA)
#     node1_hap2Masked <- replace(node1_hap2All, seq(node1_row +1, junction_length), NA)
#   }
#   else{
#     node1_hap1_normal <- gene_map1$hap1
#     node1_hap1_CGswap <- swap(node1_hap1_normal, "C", "G")
#     node1_hap1All <- swap(node1_hap1_CGswap, "A", "T")
    
#     node1_hap2_normal <- gene_map1$hap2
#     node1_hap2_CGswap <- swap(node1_hap2_normal, "C", "G")
#     node1_hap2All <- swap(node1_hap2_CGswap, "A", "T")

#     node1_row <- match(svJunctions[[index]]$sv$POS_1, svJunctions[[index]]$assembly$consensusMap$junctionPos1)
#     node1_hap1Masked <- replace(node1_hap1All, seq(node1_row +1, junction_length), NA)
#     node1_hap2Masked <- replace(node1_hap2All, seq(node1_row +1, junction_length), NA)
#   }
  
#   if(gene_map2[1,1] < gene_map2[2,1]){
#     node2_hap1All <- gene_map2$hap1
#     node2_hap2All <- gene_map2$hap2
#     node2_row <- match(svJunctions[[index]]$sv$POS_2, svJunctions[[index]]$assembly$consensusMap$junctionPos2)
#     node2_hap1Masked <- replace(node2_hap1All, seq(1, node2_row-1), NA)
#     node2_hap2Masked <- replace(node2_hap2All, seq(1, node2_row-1), NA)
#   }
#   else{
#     node2_hap1_normal <- gene_map2$hap1
#     node2_hap1_CGswap <- swap(node2_hap1_normal, "C", "G")
#     node2_hap1All <- swap(node2_hap1_CGswap, "A", "T")
    
#     node2_hap2_normal <- gene_map2$hap2
#     node2_hap2_CGswap <- swap(node2_hap2_normal, "C", "G")
#     node2_hap2All <- swap(node2_hap2_CGswap, "A", "T")
    
#     node2_row <- match(svJunctions[[index]]$sv$POS_2, svJunctions[[index]]$assembly$consensusMap$junctionPos2)
#     node2_hap1Masked <- replace(node2_hap1All, seq(1, node2_row-1), NA)
#     node2_hap2Masked <- replace(node2_hap2All, seq(1, node2_row-1), NA)
#   }
  
#   new_consensus_map <- data.table(
#     refPos1 = svJunctions[[index]]$assembly$consensusMap$refPos1,
#     refPos2 = svJunctions[[index]]$assembly$consensusMap$refPos2,
#     junctionPos1 = svJunctions[[index]]$assembly$consensusMap$junctionPos1,
#     junctionPos2 = svJunctions[[index]]$assembly$consensusMap$junctionPos2,
#     nMols = svJunctions[[index]]$assembly$consensusMap$nMols,
#     consensus = svJunctions[[index]]$assembly$consensusMap$consensus,
#     node1_refAll = svJunctions[[index]]$assembly$consensusMap$ref1All,
#     node1_refMasked = svJunctions[[index]]$assembly$consensusMap$ref1Masked,
#     node2_refAll = svJunctions[[index]]$assembly$consensusMap$ref2All,
#     node2_refMasked = svJunctions[[index]]$assembly$consensusMap$ref2Masked,
#     node1_hap1All = node1_hap1All,
#     node1_hap1Masked = node1_hap1Masked,
#     node1_hap2All = node1_hap2All,
#     node1_hap2Masked = node1_hap2Masked,
#     node2_hap1All = node2_hap1All,
#     node2_hap1Masked = node2_hap1Masked,
#     node2_hap2All = node2_hap2All,
#     node2_hap2Masked = node2_hap2Masked,
#     isRef1 = svJunctions[[index]]$assembly$consensusMap$consensus == svJunctions[[index]]$assembly$consensusMap$ref1Masked,
#     isRef2 = svJunctions[[index]]$assembly$consensusMap$consensus == svJunctions[[index]]$assembly$consensusMap$ref2Masked,
#     isNode1Hap1 = svJunctions[[index]]$assembly$consensusMap$consensus == node1_hap1Masked,
#     isNode1Hap2 = svJunctions[[index]]$assembly$consensusMap$consensus == node1_hap2Masked,
#     isNode2Hap1 = svJunctions[[index]]$assembly$consensusMap$consensus == node2_hap1Masked,
#     isNode2Hap2 =svJunctions[[index]]$assembly$consensusMap$consensus == node2_hap2Masked
#   )
#   return(new_consensus_map)
# }

# assembly_mismatch <- function(index){ #fills out the number of mismatches between the consensus sequence and the reference/haplotype
#   falseCount_Ref1 = length(svJunctions[[index]]$assembly$consensusMap$isRef1) - 
#                     sum(svJunctions[[index]]$assembly$consensusMap$isRef1, na.rm = TRUE) - 
#                     sum(is.na(svJunctions[[index]]$assembly$consensusMap$isRef1))
  
#   falseCount_Ref2 = length(svJunctions[[index]]$assembly$consensusMap$isRef2) - 
#                     sum(svJunctions[[index]]$assembly$consensusMap$isRef2, na.rm = TRUE) - 
#                     sum(is.na(svJunctions[[index]]$assembly$consensusMap$isRef2))
  
#   falseCount_Node1 = 0
#   falseCount_Node2 = 0
  
#   for (rowNum in seq(1:length(svJunctions[[index]]$assembly$consensusMap))){
#     if (is.na(svJunctions[[index]]$assembly$consensusMap$isNode1Hap1[rowNum]) & is.na(svJunctions[[index]]$assembly$consensusMap$isNode1Hap2[rowNum])){
#     }
#     else if (is.na(svJunctions[[index]]$assembly$consensusMap$isNode1Hap1[rowNum]) & svJunctions[[index]]$assembly$consensusMap$isNode1Hap2[rowNum] == FALSE){
#       falseCount_Node1 <- falseCount_Node1 + 1
#     }
#     else if (is.na(svJunctions[[index]]$assembly$consensusMap$isNode1Hap2[rowNum]) & svJunctions[[index]]$assembly$consensusMap$isNode1Hap1[rowNum] == FALSE){
#       falseCount_Node1 <- falseCount_Node1 + 1
#     }
#     else if (svJunctions[[index]]$assembly$consensusMap$isNode1Hap1[rowNum] == FALSE & svJunctions[[index]]$assembly$consensusMap$isNode1Hap2[rowNum] == FALSE){
#       falseCount_Node1 <- falseCount_Node1 + 1
#     }
    
#     if (is.na(svJunctions[[index]]$assembly$consensusMap$isNode2Hap1[rowNum]) & is.na(svJunctions[[index]]$assembly$consensusMap$isNode2Hap2[rowNum])){
#     }
#     else if (is.na(svJunctions[[index]]$assembly$consensusMap$isNode2Hap1[rowNum]) & svJunctions[[index]]$assembly$consensusMap$isNode2Hap2[rowNum] == FALSE){
#       falseCount_Node2 <- falseCount_Node2 + 1
#     }
#     else if (is.na(svJunctions[[index]]$assembly$consensusMap$isNode2Hap2[rowNum]) & svJunctions[[index]]$assembly$consensusMap$isNode2Hap1[rowNum] == FALSE){
#       falseCount_Node2 <- falseCount_Node2 + 1
#     }
#     else if (svJunctions[[index]]$assembly$consensusMap$isNode2Hap1[rowNum] == FALSE & svJunctions[[index]]$assembly$consensusMap$isNode2Hap2[rowNum] == FALSE){
#       falseCount_Node2 <- falseCount_Node2 + 1
#     }
#   }
#   new_assembly <- list(
#     success = svJunctions[[index]]$assembly$success,
#     readMap = svJunctions[[index]]$assembly$readMap,
#     consensusMap = svJunctions[[index]]$assembly$consensusMap,
#     maxNMols = svJunctions[[index]]$assembly$maxNMols,
#     nMismatchNode1_Ref = falseCount_Ref1,
#     nMismatchNode2_Ref = falseCount_Ref2,
#     nMismatchNode1_Hap = falseCount_Node1,
#     nMismatchNode2_Hap = falseCount_Node2
#   )
#   return(new_assembly)
# }

# load(file = "dna01337_6_0pt6.hg38.assemble.junctions.RData")

# for (i in seq(1:length(svJunctions))){  #main body of code, uses the above three functions to create a new consensus map and replaces the original one
#   jxn_type <- svJunctions[[i]]$sv$JXN_TYPE
#   jxn_location <- svJunctions[[i]]$sv$CHROM_1
#   if (jxn_location == "chr1"){  #negr1
#     if (jxn_type == "L" | jxn_type == "D"){ #deletion
#       svJunctions[[i]]$assembly$consensusMap <- sv_del_dup(hf1_negr1, 71430000, i)
#       svJunctions[[i]]$assembly <- assembly_mismatch(i)
#     }
#     else if (jxn_type == "I"){ #inversion
#       svJunctions[[i]]$assembly$consensusMap <- sv_inversion(hf1_negr1, 71430000, i)
#       svJunctions[[i]]$assembly <- assembly_mismatch(i)
#     }
#     else{
#       print(c("new type", jxn_type))
#     }
#   }
#   else if (jxn_location == "chr10"){ #prkg1
#     if (jxn_type == "L" | jxn_type == "D"){ #deletion
#       svJunctions[[i]]$assembly$consensusMap <- sv_del_dup(hf1_prkg1, 51290000, i)
#       svJunctions[[i]]$assembly <- assembly_mismatch(i)
#     }
#     else if (jxn_type == "I"){ #inversion
#       svJunctions[[i]]$assembly$consensusMap <- sv_inversion(hf1_prkg1, 51290000, i)
#       svJunctions[[i]]$assembly <- assembly_mismatch(i)
#     }
#     else{
#       print(c("new type", jxn_type))
#     }
#   }
#   else if (jxn_location == "chr7"){ #magi2
#     if (jxn_type == "L" | jxn_type == "D"){ #deletion
#       svJunctions[[i]]$assembly$consensusMap <- sv_del_dup(hf1_magi2, 78210000, i)
#       svJunctions[[i]]$assembly <- assembly_mismatch(i)
#     }
#     else if (jxn_type == "I"){ #inversion
#       svJunctions[[i]]$assembly$consensusMap <- sv_inversion(hf1_magi2, 78210000, i)
#       svJunctions[[i]]$assembly <- assembly_mismatch(i)
#     }
#     else{
#       print(c("new type", jxn_type))
#     }
#   }
#   else{
#     print("gene error")
#   }
# }

# haplotype_mismatches <- list() #list of junctions with haplotype mismatch
# only_hap_mismatches <- list() #list of junctions with haplotype mismatch, but reference match 
# for (i in seq(1:length(svJunctions))){ #returns above two 
#   if (svJunctions[[i]]$assembly$nMismatchNode1_Hap > 0 | svJunctions[[i]]$assembly$nMismatchNode2_Hap > 0){
#     haplotype_mismatches <- append(haplotype_mismatches, list(i, svJunctions[[i]]$sv$JXN_TYPE, svJunctions[[i]]$sv$CHROM_1))
#     if (svJunctions[[i]]$assembly$nMismatchNode1_Ref == 0 & svJunctions[[i]]$assembly$nMismatchNode2_Ref == 0){
#       only_hap_mismatches <- append(only_hap_mismatches, list(i, svJunctions[[i]]$sv$JXN_TYPE, svJunctions[[i]]$sv$CHROM_1))
#     }
#   }
# }

#save(junctionsSummary, svJunctions, file = "DNA01336_HF1_0pt6_APH.hg38.assemble.junctions.RData") #optional save command for export 

