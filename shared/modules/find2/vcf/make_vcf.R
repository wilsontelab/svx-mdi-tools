# export svCalls to simplified VCF format for upload and other programs

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
message("  initializing")
library(data.table)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'ACTION_DIR',
        'FIND_PREFIX',
        'SAMPLES',
        'GENOME',
        'GENOME_CHROMS',
        'GENOME_FASTA',
        'CHROM_FASTA_DIR'       
    ),
    integer = c(
        'N_CPU',
        'MIN_MAPQ_ONE',
        'MIN_MAPQ_BOTH',
        'MIN_SV_SIZE',
        'SV_SIZE_FACTOR',
        'PURGE_DISTANCE',
        'PURGE_LIMIT',
        'MIN_COVERAGE',
        'MIN_MERGE_OVERLAP'
    ),
    double = c(
        'MIN_MERGE_DENSITY'
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
sourceScripts(file.path(rUtilDir, 'sequence'), c('general'))
sourceScripts(file.path(rUtilDir, 'genome'),   c('faidx'))
loadFaidx()
#-------------------------------------------------------------------------------------
# initalize targets
targetsBed <- if(is.null(env$TARGETS_BED) || 
                 env$TARGETS_BED == "" || 
                 env$TARGETS_BED == "NA" || 
                 env$TARGETS_BED == "null") NA else basename(env$TARGETS_BED)
targets <- if(is.na(targetsBed)) data.table() else fread(env$TARGETS_BED)
nTargets <- nrow(targets)
targets <- if(nTargets > 100) ">100" else {
    paste(apply(targets, 1, function(v) paste0(v[4], "::", v[1], ":", v[2], "-", v[3])), collapse = ",")
}
#-------------------------------------------------------------------------------------
# initialize samples
env$SAMPLES <- sort(strsplit(env$SAMPLES, "\\s+")[[1]])
#=====================================================================================

#-------------------------------------------------------------------------------------
# write a file with most information from the svCalls object in VCF breakend format
# thus, every called junction occupies two VCF output lines, one for each breakend
#-------------------------------------------------------------------------------------

# initialize VCF
message("  parsing VCF header")
headerLine <- function(key, value){
    if(is.list(value)) value <- {
        x <- sapply(seq_along(value), function(i) paste(names(value)[i], value[[i]], sep = "=") )
        paste0("<", paste(x, collapse = ","), ">")
    }
    writeLines(paste0("##", key, "=", value))
}

# add general header lines
headerLine("fileformat", "VCFv4.3")
headerLine("fileDate",   format(Sys.time(), "%Y%m%d"))
headerLine("reference",  env$GENOME) # not using file://<local path> to conform to repository upload standards
headerLine("source",     "https://github.com/wilsontelab/svx-mdi-tools/tree/main/pipelines/svCapture")
headerLine("samples",    paste(env$SAMPLES, collapse=","))
headerLine("nSamples",   length(env$SAMPLES))
headerLine("targets",    targets)
headerLine("nTargets",   nTargets)

# add OPTION header lines
options <- list(
    REGION_PADDING      = "this many bp on each side of target regions are considered adjacent",
    MIN_MERGE_DENSITY   = "fraction of matching overlap bases to accept an alignment-guided merge",
    MIN_MERGE_OVERLAP   = "No. of bases that must overlap between read ends for alignment-guided merging",
    PURGE_DISTANCE      = "presume molecule duplicates when OUT_POS euclidean distance < --purge-distance",
    PURGE_LIMIT         = "randomly downsample evidence molecules to --purge-limit prior to duplicate purging",
    MIN_MAPQ_ONE        = 'at least one alignment flanking an SV junction must have at least this MAPQ',
    MIN_MAPQ_BOTH       = 'both alignments flanking an SV junction must have at least this MAPQ',
    MIN_SV_SIZE         = 'SVs smaller than this predicted size will not be processed',
    SV_SIZE_FACTOR      = 'if >0, override --min-sv-size to be --sv-size-factor * maxTLen',
    MIN_COVERAGE        = 'require at least this many independent molecules over all samples to call an SV',
    ON_TARGET           = 'if applicable, at least this many SV ends must be in a padded target region'
)
for(i in seq_along(options)){
    headerLine("option", list(
        ID      = names(options[i]),
        Value   = paste0('"', env[names(options[i])], '"'),
        Description = paste0('"', options[i], '"')
    ))
}

# add FILTER header lines
filters <- list(
    "MIN_MAPQ_ONE",
    "MIN_MAPQ_BOTH",
    "MIN_SV_SIZE",
    "MIN_COVERAGE",
    "ON_TARGET"
)
for(filter in filters){
    headerLine("FILTER", list(
        ID = filter,
        Description = paste0('"', 'see matching option', '"')
    ))
}

# add INFO header lines
info <- fread(file.path(env$ACTION_DIR, "vcf/INFO.csv"))
for(i in seq_along(info[[1]])){
    headerLine("INFO", list(
        ID      = info[i, VCF_ID],
        Number  = info[i, Number],
        Type    = info[i, Type],
        Description = paste0('"', info[i, Description], '"')
    ))
}

# add FORMAT header lines
headerLine("FORMAT", list(
    ID      = "AD",
    Number  = "R",
    Type    = "Integer",
    Description = paste0('"', "Number of source DNA molecules for each allele", '"')
))

# add contig header lines
fai <- fread(paste(env$GENOME_FASTA, "fai", sep = "."))
for(i in seq_along(fai[[1]])){
    headerLine("contig", list(
        ID       = fai[i, 1],
        length   = fai[i, 2],
        assembly = env$GENOME
    ))
}

# add the header line for variant call records
variantColumns <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")
allColumns <- paste(c(variantColumns, env$SAMPLES))
writeLines(paste0("#", paste(allColumns, collapse = "\t")))

# parse INFO for two breakend lines per variant record
message("  parsing variant records")
inFile  <- paste(env$FIND_PREFIX, 'structural_variants', 'gz',  sep = ".")
svCalls <- fread(inFile)
# svCalls <- svCalls[JXN_TYPE == "I" & MICROHOM_LEN < 0 & SIDE_1 == "R" & SIDE_2 == "R" & CHROM_1 == "chr3"][1:3]
callInfo1 <- data.table(TEMPORARY = svCalls$SV_ID)
getMQ <- function(x) {
    x <- as.integer(strsplit(x, ",")[[1]])
    sqrt(mean(x^2)) # root mean square MAP for all alignments on a side
}
for(i in seq_along(info$VCF_ID)) callInfo1[[info$VCF_ID[i]]] <- switch( # fill junction side 1
    info$VCF_ID[i],    
        MATEID = paste(svCalls$SV_ID, 2, sep = ":"), # these values need format conversion
        MQ = getMQ(svCalls$MAPQ_1),
        SVTYPE = "BND",
        SVLEN = ifelse(svCalls$JXN_TYPE == "L", -svCalls$SV_SIZE, svCalls$SV_SIZE), # DEL neg SVLEN per VCF format
        HOMLEN = as.integer(pmax(0, svCalls$MICROHOM_LEN)), # only micromology appropriate here, not de novo insertions
        HOMSEQ = ifelse(svCalls$MICROHOM_LEN <= 0, ".", svCalls$JXN_BASES),
        JXNSEQ = ifelse(svCalls$JXN_SEQ == "*", ".", svCalls$JXN_SEQ),
        TGT = gsub(",", "|", svCalls$TARGET_REGION),
            svCalls[[info$SVX_ID[i]]] # simple transfer of all other values
)
callInfo2 <- copy(callInfo1) # copy to junction side 2
callInfo2[, ":="( # modify columns that differ between junction sides 1 and 2
    MATEID = paste(svCalls$SV_ID, 1, sep = ":"),
    MQ = getMQ(svCalls$MAPQ_2)
)]
callInfo <- rbind(callInfo1, callInfo2) # merge to table of all records, two per junction
rm(callInfo1, callInfo2)
for(i in seq_along(info$VCF_ID)) if(info$Number[i] == "R") { # add reference allele coverage
    # TODO: update to add coverage data for reference allele, probably just for svWGS
    callInfo[[info$VCF_ID[i]]] <- paste(".", callInfo[[info$VCF_ID[i]]], sep = ",")
}

# parse the remainder of the variant record lines
records1 <- data.table(TEMPORARY = svCalls$SV_ID)
nodeBases <- mapply(getRefSeq_padded, svCalls$CHROM_1, svCalls$POS_1, 0)
getALT <- function(THIS_SIDE_I, THIS_SIDE, 
                   OTHER_CHROM, OTHER_POS, OTHER_SIDE, 
                   MICROHOM_LEN, JXN_BASES){
    otherNode <- paste(OTHER_CHROM, OTHER_POS, sep = ":")
    otherDir <- ifelse(OTHER_SIDE == "L", "]", "[")
    otherNode <- paste0(otherDir, otherNode, otherDir)
    insertionSeq <- ifelse(MICROHOM_LEN < 0, JXN_BASES, "") # microhomologies are implicit in the node positions
    ifelse(
        THIS_SIDE == OTHER_SIDE,

        # inversion type along conjoined chromosomes, side 1 = leftmost position
        # must rc the insertion sequence on the non-canonical strand
        if(THIS_SIDE_I == 1) ifelse(
            THIS_SIDE == "L",
            # t]p] reverse comp piece extending left of p is joined after t
            paste0(paste0(nodeBases, insertionSeq), otherNode),    # this side is the canonical strand
            # s [p[t reverse comp piece extending right of p is joined before t
            paste0(otherNode, paste0(rc(insertionSeq), nodeBases)) # other side is the canonical strand

        ) else ifelse(
            THIS_SIDE == "L",
            # t]p] reverse comp piece extending left of p is joined after t
            paste0(paste0(nodeBases, rc(insertionSeq)), otherNode), # other side is the canonical strand
            # s [p[t reverse comp piece extending right of p is joined before t
            paste0(otherNode, paste0(insertionSeq, nodeBases))      # this side is the canonical strand
        ), 

        # PDL type along conjoined chromosomes, side 1 == leftward alignment
        # both sides are the canonical strand, rc not needed
                             # s t[p[ piece extending to the right of p is joined after t
        if(THIS_SIDE_I == 1) paste0(paste0(nodeBases, insertionSeq), otherNode) # THIS_SIDE == L, OTHER_SIDE = R
                             # s ]p]t piece extending to the left of p is joined before t
                        else paste0(otherNode, paste0(insertionSeq, nodeBases)) # THIS_SIDE == R, OTHER_SIDE = L
    )
}
records1[, ":="( # parse breakend 1
    CHROM = svCalls$CHROM_1,
    POS = svCalls$POS_1,
    ID = paste(svCalls$SV_ID, 1, sep = ":"),
    REF = nodeBases,
    ALT = getALT(1, svCalls$SIDE_1, 
                 svCalls$CHROM_2, svCalls$POS_2, svCalls$SIDE_2, 
                 svCalls$MICROHOM_LEN, svCalls$JXN_BASES)
)]
for(sample in env$SAMPLES) records1[[sample]] <- svCalls[[sample]] # add sample molecule counts
records2 <- data.table(TEMPORARY = svCalls$SV_ID)
nodeBases <- mapply(getRefSeq_padded, svCalls$CHROM_2, svCalls$POS_2, 0)
records2[, ":="( # parse breakend 2
    CHROM = svCalls$CHROM_2,
    POS = svCalls$POS_2,
    ID = paste(svCalls$SV_ID, 2, sep = ":"),
    REF = nodeBases,
    ALT = getALT(2, svCalls$SIDE_2, 
                 svCalls$CHROM_1, svCalls$POS_1, svCalls$SIDE_1, 
                 svCalls$MICROHOM_LEN, svCalls$JXN_BASES)
)]
for(sample in env$SAMPLES) records2[[sample]] <- svCalls[[sample]] # add sample molecule counts

# merge to final table of all records, two per junction
records <- rbind(records1, records2)
rm(records1, records2)
records[, ":="(
    QUAL = ".",
    FILTER = "PASS",
    INFO = apply(callInfo[, .SD, .SDcols = info$VCF_ID], 1, function(v){
        paste(paste(info$VCF_ID, trimws(v), sep = "="), collapse = ";")
    }),
    FORMAT = "AD"
)]

# print for handling by bgzip/tabix
write.table(
    records[, .SD, .SDcols = allColumns][order(CHROM, POS)], 
    file = "", 
    quote = FALSE, 
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE
)
#=====================================================================================

# ##fileformat=VCFv4.3
# ##fileDate=20090805
# ##source=myImputationProgramV3.1
# ##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
# ##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
# ##phasing=partial
# ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
# ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
# ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
# ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
# ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
# ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
# ##FILTER=<ID=q10,Description="Quality below 10">
# ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
# ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
# ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
# #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
# 20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
# 20 17330 . T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
# 20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
# 20 1230237 . T . 47 PASS NS=3;DP=13;AA=T GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
# 20 1234567 microsat1 GTC G,GTCT 50 PASS NS=3;DP=9;AA=G GT:GQ:DP 0/1:35:4 0/2:17:2 1/1:40:3

# The REF field of a breakend record indicates a base or
# sequence s of bases beginning at position POS, as in all VCF records. The ALT field of a breakend record indicates
# a replacement for s. This “breakend replacement” has three parts:
#   1. The string t that replaces places s. The string t may be an extended version of s if some novel bases are inserted
#      during the formation of the novel adjacency.
#   2. The position p of the mate breakend, indicated by a string of the form “chr:pos”. This is the location of the
#      first mapped base in the piece being joined at this novel adjacency.
#   3. The direction that the joined sequence continues in, starting from p. This is indicated by the orientation of
#      square brackets surrounding p.
# These 3 elements are combined in 4 possible ways to create the ALT. In each of the 4 cases, the assertion is that s
# is replaced with t, and then some piece starting at position p is joined to t. The cases are:
#   REF ALT  Meaning
#   s   t[p[ piece extending to the right of p is joined after t
#   s   t]p] reverse comp piece extending left of p is joined after t
#   s   ]p]t piece extending to the left of p is joined before t
#   s   [p[t reverse comp piece extending right of p is joined before t

#  $ SV_ID                    : chr  "10000:1" "10001:1" "10002:1" "10003:1" ...
#  $ MAPQ_1                   : chr  "60" "60" "60" "60" ...
#  $ UMI_1                    : chr  "1" "1" "1" "1" ...
#  $ N_TOTAL                  : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ N_GAPS                   : int  1 0 1 1 0 0 0 0 1 0 ...
#  $ N_SPLITS                 : int  0 1 0 0 1 1 1 1 0 1 ...
#  $ N_OUTER_CLIPS            : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ JXN_TYPE                 : chr  "L" "L" "L" "L" ...
#  $ N_DUPLEX                 : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ N_DUPLEX_GS              : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ STRAND_COUNT             : int  17 3 9 14 8 2 7 19 16 1 ...
#  $ STRAND_COUNT_GS          : int  17 3 9 14 8 2 7 19 16 1 ...
#  $ STRAND_COUNT1            : int  17 3 0 14 8 0 0 19 0 0 ...
#  $ STRAND_COUNT2            : int  0 0 9 0 0 2 7 0 16 1 ...
#  $ TARGET_CLASS             : chr  "TA" "TA" "TA" "TA" ...
#  $ SHARED_PROPER            : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ SHARED_PROPER_GS         : num  1 1 1 1 1 1 1 1 1 1 ...
#  $ SAMPLES                  : chr  "HCT_0.2APH_Ro3_Colch_M_a" "HCT_0.2APH_Ro3_Colch_M_b" "HCT_0.2APH_Ro3_Colch_M_b" "HCT_0.2APH_Ro3_Colch_M_b" ...
#  $ N_SAMPLES                : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ MAPQ_2                   : chr  "60" "60" "60" "60" ...
#  $ UMI_2                    : chr  "1" "1" "1" "1" ...
#  $ CHROM_1                  : chr  "chr16" "chr16" "chr16" "chr16" ...
#  $ SIDE_1                   : chr  "L" "L" "L" "L" ...
#  $ POS_1                    : int  78366477 78368021 78368689 78369692 78374024 78373165 78375043 78375734 78376263 78377504 ...
#  $ CHROM_2                  : chr  "chr16" "chr16" "chr16" "chr16" ...
#  $ SIDE_2                   : chr  "R" "R" "R" "R" ...
#  $ POS_2                    : int  78751195 78480682 78539034 78739838 78463200 78522163 78542833 78705043 78734570 78718865 ...
#  $ JUNCTION_NAME            : chr  "16:L:78366477,16:R:78751195" "16:L:78368021,16:R:78480682" "16:L:78368689,16:R:78539034" "16:L:78369692,16:R:78739838" ...
#  $ JUNCTION_NAMES           : chr  "16:L:78366477,16:R:78751195" "16:L:78368021,16:R:78480682" "16:L:78368689,16:R:78539034" "16:L:78369692,16:R:78739838" ...
#  $ N_AMBIGUOUS              : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ N_DOWNSAMPLED            : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ N_COLLAPSED              : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ JXN_SEQ                  : chr  "*" "GATTACAGGTGCGTGACACCACACCCCGCTAATTTTTGTATTTTTAGTGGAGACGGGATTTCACCATGTTGGTCAGGCTGGTGTCAAACTCCTGATCTTAGGTGATCTGCC"| __truncated__ "*" "*" ...
#  $ MERGE_LEN                : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ MICROHOM_LEN             : int  0 3 0 0 0 -3 4 -2 0 1 ...
#  $ JXN_BASES                : chr  "*" "AGT" "*" "*" ...
#  $ SV_SIZE                  : int  384718 112661 170345 370146 89176 148998 167790 329309 358307 341361 ...
#  $ GEN_REF_1                : chr  "TGGCACACAGAGCATCCCACCAATAAACTGTGGTTCCCAACAAGGCTTGCAACTGAATCTGGGCAACAATGAACATGACTTCAATAATGACAAAAAGATTCAGATTCCTGC"| __truncated__ "TTACATTTTTTTTTTTTTTTTTTTTTTTTTTTTAAGACGTAGTCTTGCTCTGTTGCCCAGGCTGGAGTGCAGTGGTGTGATTTCAGCTCACAGCAAGCTCTGCCTCCCAGG"| __truncated__ "CTCAGGGATAATTTCCAGGCTGGGGAGAAGGTTGTGTCATGGAGACAGCAGTACTGTTTTGGAAGGAGAATTGACTTGTCACAATGGACGGAGAACTACCATCCCATTTGG"| __truncated__ "ATGCCATCGTGAAGTGGCTTCTAATGGACATTTGTAAGGTGGTGTCCCCTCCCTTCCCCTCCCCCCTCTAATGCTGCCTGGGTTGCTCTTTGCAACACGTGCAGTCTGCAA"| __truncated__ ...
#  $ GEN_REF_2                : chr  "GTGAGCAACACCTCGGTGAAAGTAGAGGTGTTTGTACAGTTCAGCTGCTGCCTTATTCAAAAACAGGAAGAAAAACCCACTGCCCAGGTCAGCTGGTGGGTGATCTCCACT"| __truncated__ "CCAAGGAGGAAACTGTTATAGGAATGATTGACATGGAGGTTCATCTGTAGAATATGGACAGCCCTCTACAAATTGGGAACAAGGGCCAGGACACGAAATATAGTGACAGCC"| __truncated__ "TTCACCTTTCCCCGATTCTTGAAAACGGAAAAGATGAGGGTTGAATTATAGGGCAAATCATGAACACTGTCTTTTTTCCAGCATGGTATTTTGAGCGGCACAGAAGGGGGT"| __truncated__ "AGAGCCATTGCCTTGACCACAAAGCCACAGTCTCCTAGTCCTAATATTAAGAGGTAGGAGGTAGATAGAGATTGTCCCTTAAGCCACATTCCATTCTTCATATTTTTGTAC"| __truncated__ ...
#  $ GEN_COV_1                : logi  NA NA NA NA NA NA ...
#  $ GEN_COV_2                : logi  NA NA NA NA NA NA ...
#  $ TARGET_REGION            : chr  "WWOX" "WWOX" "WWOX" "WWOX" ...
#  $ TARGET_POS_1             : int  2710437 2711981 2712649 2713652 2717984 2717125 2719003 2719694 2720223 2721464 ...
#  $ TARGET_POS_2             : int  3095155 2824642 2882994 3083798 2807160 2866123 2886793 3049003 3078530 3062825 ...
#  $ HCT_0.2APH_Ro3_Colch_G2_a: int  0 0 0 0 0 0 0 0 0 0 ...
#  $ HCT_0.2APH_Ro3_Colch_G2_b: int  0 0 0 0 0 0 0 1 0 0 ...
#  $ HCT_0.2APH_Ro3_Colch_M_a : int  1 0 0 0 0 0 0 0 0 0 ...
#  $ HCT_0.2APH_Ro3_Colch_M_b : int  0 1 1 1 1 1 1 0 1 1 ...
#  $ HCT_Ro3_Colch_G2_a       : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ HCT_Ro3_Colch_M_a        : int  0 0 0 0 0 0 0 0 0 0 ...
