# assemble and tabulate SVs across multiple previously executed find actions

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
    library(yaml)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
myEnvVars <- list(
    string = c(
        'ACTION_DIR',
        'GENOMEX_MODULES_DIR',
        'INPUT_DIR',
        'SAMPLES_TABLE',
        'OUTPUT_DIR',
        'DATA_NAME',
        'DATA_FILE_PREFIX',
        'ASSEMBLE_PREFIX'
    ),
    integer = c(
        'N_CPU',
        'MIN_FAMILY_SIZE', 
        'MIN_SV_SIZE',
        'MAX_SV_SIZE',
        'MIN_SAMPLES',
        'MAX_SAMPLES',
        'MIN_MOLECULES',
        'MAX_MOLECULES',
        'MIN_SPLIT_READS',
        'MIN_INSERTION_SIZE',
        'MAX_INSERTION_SIZE',
        'FLANKING_MICROHOMOLOGY',
        'INSERTION_SEARCH_SPACE'
    ),
    double = c(
        'MAX_SHARED_PROPER'
    ),
    logical = c(
        'WARN_ONLY',
        'INCLUDE_CLIPS',
        'FIND_INSERTION_TEMPLATES'
    )
)
checkEnvVars(myEnvVars)
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(rUtilDir, 'utilities')
rUtilDir <- file.path(env$GENOMEX_MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'sequence'), c('general')) # ,'IUPAC','smith_waterman'
#=====================================================================================

#=====================================================================================
# utility functions
#-------------------------------------------------------------------------------------
getProjectDir <- function(project) file.path(env$INPUT_DIR, project)
getSampleDir  <- function(project, sample) file.path(getProjectDir(project), sample)
getSampleDir2 <- function(sample) getSampleDir(sample$project, sample$sample)
getSampleExtractTask <- function(project, sample){
    taskFile <- paste0(sample, ".svCapture.extract.task.log")
    taskFile <- file.path(getSampleDir(project, sample), "svCapture/extract/logs", taskFile)
    if(!file.exists(taskFile)) stop("file not found:", taskFile)
    read_yaml(taskFile)
}
getSampleMetadata  <- function(column, project_, sample_) samples[project == project_ & sample == sample_][[column]]
getSampleMetadata2 <- function(column, sample) getSampleMetadata(column, sample$project, sample$sample)
getSampleBamFile <- function(sample){
    genome <- getSampleMetadata2("genome", sample)
    suffix <- if(getSampleMetadata2("cram", sample)) "cram" else "bam"
    bamFile <- paste(sample$sample, genome, "name.realigned", suffix, sep = ".")
    bamFile <- file.path(getSampleDir2(sample), bamFile)
    if(!file.exists(bamFile)) stop("file not found:", bamFile)
    bamFile
}
getProjectFindFile <- function(project_, sample_){
    genome <- getSampleMetadata("genome", project_, sample_)
    rdsFile <- paste(project_, genome, "find.structural_variants.rds", sep = ".")
    rdsFile <- file.path(getProjectDir(project_), rdsFile)
    rdsFile
}
#=====================================================================================

#=====================================================================================
message("loading samples table")
samples <- fread(env$SAMPLES_TABLE, header = TRUE)
samples[, i := 1:.N]

message("collecting samples metadata")
samples <- merge(
    samples, 
    samples[, {
        task <- getSampleExtractTask(project, sample)$extract
        .(
            genomesDir      = task$genome[["genomes-dir"]],
            genome          = task$genome$genome,
            cram            = as.logical(task[["bam-format"]][["use-cram"]]),
            minMapq         = as.integer(task[["read-filtering"]][["min-mapq"]]),
            targetsBed      = task[["target-region"]][["targets-bed"]],
            regionPadding   = as.integer(task[["target-region"]][["region-padding"]])
        )
    }, by = .(i, project, sample)],
    by = c("i", "project", "sample"),
    all.x = TRUE 
)

message("cross-checking samples metadata")
passedCrosscheck <- TRUE
for(column in c("genome","targetsBed","regionPadding","minMapq")){
    values <- unique(samples[[column]])
    if(length(values) > 1){
        message()
        message(paste("!!!!!!!!!! metadata mismatch between samples on parameter",  column,  "!!!!!!!!!!"))
        message("found values are:")
        print(values)
        passedCrosscheck <- FALSE
    }
}
if(!passedCrosscheck && !env$WARN_ONLY) stop("crosscheck failed, aborting since --warn-only not set")
#=====================================================================================

#=====================================================================================
message("getting N50 and adjusted target coverage by sample")
nSamples <- nrow(samples)
samples <- merge(
    samples,
    do.call(rbind, mclapply(1:nSamples, function(i){
        sample <- samples[i]
        message(paste("", sample$project, sample$sample, sep = "\t"))
        x <- system2(
            "perl", 
            args = c(
                file.path(env$ACTION_DIR, "coverage.pl"),
                sapply(c("genome", "targetsBed", "regionPadding", "minMapq"), getSampleMetadata2, sample),
                getSampleBamFile(sample)
            ),
            stdout = TRUE, 
            stderr = FALSE
        )
        cbind(sample[, .(i, project, sample)], fread(text = x))
    }, mc.cores = env$N_CPU)),
    by = c("i","project","sample"),
    all.x = TRUE
) 
#=====================================================================================

#=====================================================================================
message("tabulating and filtering SVs across all projects and samples")
applySVFilters <- function(svs){
    svs <- svs[
        (JXN_TYPE == "T" | 
            (SV_SIZE >= env$MIN_SV_SIZE & 
             SV_SIZE <= env$MAX_SV_SIZE)
        ) &
        N_SAMPLES >= env$MIN_SAMPLES &
        N_SAMPLES <= env$MAX_SAMPLES &
        N_SPLITS >= env$MIN_SPLIT_READS & 
        TARGET_CLASS %in% c("TT", "TA") & # assemble only ever returns on-target SVs
        STRAND_COUNT >= env$MIN_FAMILY_SIZE & 
        SHARED_PROPER / 2 <= env$MAX_SHARED_PROPER
    ]   
    nTotal <- if(env$INCLUDE_CLIPS) svs$N_TOTAL else svs$N_SPLITS + svs$N_GAPS
    svs[
        nTotal >= env$MIN_MOLECULES &
        nTotal <= env$MAX_MOLECULES
    ]
}
svs <- samples[, {
    message(paste("", project, sep = "\t"))
    projectFindFile <- getProjectFindFile(project, sample[1])
    isMultiSampleFind <- file.exists(projectFindFile)
    if(!isMultiSampleFind) stop("assembly of single-sample find inputs is not yet implemented")
    applySVFilters(readRDS(projectFindFile))[, .(
        SV_ID,
        TARGET_REGION,
        TARGET_CLASS,
        JXN_TYPE,
        SV_SIZE,
        CHROM_1,
        SIDE_1,
        POS_1,
        CHROM_2,
        SIDE_2,
        POS_2,
        N_MOLECULES   = if(env$INCLUDE_CLIPS) N_TOTAL else N_SPLITS + N_GAPS,
        N_SPLITS,
        N_DUPLEX      = if(env$INCLUDE_CLIPS) N_DUPLEX else N_DUPLEX_GS,
        STRAND_COUNT  = if(env$INCLUDE_CLIPS) STRAND_COUNT else STRAND_COUNT_GS, # i.e., family size
        SHARED_PROPER = if(env$INCLUDE_CLIPS) SHARED_PROPER else SHARED_PROPER_GS, 
        N_SAMPLES,
        SAMPLES = paste0(",", SAMPLES, ","),
        MICROHOM_LEN,
        JXN_BASES
    )]
}, by = .(project)]
names(svs)[1] <- "PROJECT"

message("counting filtered SVs by sample and type")
samples <- merge(
    samples,
    samples[, {
        jxnTypes <- svs[PROJECT == project & grepl(paste0(",", sample, ","), SAMPLES), JXN_TYPE]
        .(
            deletion      = sum(jxnTypes == "L"),
            duplication   = sum(jxnTypes == "D"),
            inversion     = sum(jxnTypes == "I"),
            translocation = sum(jxnTypes == "T")
        )
    }, by = .(i, project, sample)],
    by = c("i", "project", "sample"),
    all.x = TRUE
)
#=====================================================================================

#=====================================================================================
message("assembling target region sequences")
targets <- samples[, {
    targets <- fread(targetsBed)[, 1:4]
    names(targets) <- c("chrom","regionStart","regionEnd","name")
    regionPadding <- max(regionPadding)
    genomeFasta <- file.path(genomesDir[1], "iGenomes/*/UCSC", genome, "Sequence/WholeGenomeFasta/genome.fa")
    targets[, ":="(
        regionPadding = regionPadding,
        paddedStart = regionStart - regionPadding, # all still half-open like the source BED
        paddedEnd = regionEnd + regionPadding
    )]
    targets[, ":="(
        paddedSequence = {
            x <- system2(
                "samtools",
                args = c(
                    "faidx",
                    genomeFasta,
                    paste0(chrom, ":", paddedStart + 1, "-", paddedEnd)
                ),
                stdout = TRUE
            )
            paste(toupper(x[2:length(x)]), collapse = "")
        }
    ), by = .(chrom, regionStart, regionEnd, name)]
    targets
}, by = .(genome, targetsBed)]
#=====================================================================================

#=====================================================================================
if(env$FIND_INSERTION_TEMPLATES) message("attempting to find insertion templates near junctions")
nullInsertion <- data.table(
    templateBreakpointN = NA_integer_,
    templateBreakpointSide = NA_character_,
    templateType = NA_character_, # "not_found", "foldback", or "cross_junction"
    templateDistance = NA_integer_
)
notFoundInsertion <- data.table(
    templateBreakpointN = NA_integer_,
    templateBreakpointSide = NA_character_,
    templateType = "not_found", 
    templateDistance = NA_integer_
)
isUsablePosition <- function(pos, target) pos - target$paddedStart > env$INSERTION_SEARCH_SPACE && 
                                           target$paddedEnd - pos >= env$INSERTION_SEARCH_SPACE
analyzeInsertion <- function(project_, sv){

    # parse the genome reference target
    targetsBed_ <- samples[project == project_ & sample == strsplit(sv$SAMPLES, ",")[[1]][2], targetsBed] # all samples in a project were analyzed with the same targets
    target <- targets[targetsBed == targetsBed_ & name == sv$TARGET_REGION]
    if(!isUsablePosition(sv$POS_1, target) || 
       !isUsablePosition(sv$POS_2, target)) return(nullInsertion)

    # get the four regions to search as the retained and lost (to this junction) segments on each side of each breakpoint
    # matches must land wholly within one of these four regions
    ss <- env$INSERTION_SEARCH_SPACE
    fm <- env$FLANKING_MICROHOMOLOGY
    tps <- target$paddedStart
    if(sv$SIDE_1 == "L"){
        genRefRetained1 <-    substr(target$paddedSequence, sv$POS_1 - ss + 1 - tps, sv$POS_1 - tps)
        genRefLost1     <-    substr(target$paddedSequence, sv$POS_1 + 1 - tps, sv$POS_1 + ss - tps)
    } else {
        genRefRetained1 <- rc(substr(target$paddedSequence, sv$POS_1 - tps, sv$POS_1 + ss - 1 - tps)) # flip one reference strand for inversions to match jxnSeq assembly
        genRefLost1     <- rc(substr(target$paddedSequence, sv$POS_1 - ss - tps, sv$POS_1 - 1 - tps))
    }
    if(sv$SIDE_2 == "R"){
        genRefRetained2 <-    substr(target$paddedSequence, sv$POS_2 - tps, sv$POS_2 + ss - 1 - tps)
        genRefLost2     <-    substr(target$paddedSequence, sv$POS_2 - ss - tps, sv$POS_2 - 1 - tps)
    } else {
        genRefRetained2 <- rc(substr(target$paddedSequence, sv$POS_2 - ss + 1 - tps, sv$POS_2 - tps))
        genRefLost2     <- rc(substr(target$paddedSequence, sv$POS_2 + 1 - tps, sv$POS_2 + ss - tps))
    }

    # parse the insertion search sequence as (microhomology)(insertion)(microhomology)
    microhomology1 <- substr(genRefRetained1, ss - fm + 1, ss)
    microhomology2 <- substr(genRefRetained2, 1, fm)
    searchSeq <- paste0(microhomology1, sv$JXN_BASES, microhomology2)
    rcSearchSeq <- rc(searchSeq)

    # find all matches
    matches <- list(
        match1R   = gregexpr(  searchSeq, genRefRetained1)[[1]],
        match1Rrc = gregexpr(rcSearchSeq, genRefRetained1)[[1]],
        match1L   = gregexpr(  searchSeq, genRefLost1)[[1]],
        match1Lrc = gregexpr(rcSearchSeq, genRefLost1)[[1]],
        match2R   = gregexpr(  searchSeq, genRefRetained2)[[1]],
        match2Rrc = gregexpr(rcSearchSeq, genRefRetained2)[[1]],
        match2L   = gregexpr(  searchSeq, genRefLost2)[[1]],
        match2Lrc = gregexpr(rcSearchSeq, genRefLost2)[[1]]
    )
    if(all(unlist(matches) == -1)) return(notFoundInsertion)

    # characterize all template matches
    searchSeqLength <- nchar(searchSeq)
    parseMatches <- function(matchStartPositions, breakpointN, templateType, breakpointSide, orientation){
        if(all(matchStartPositions == -1)) return(NULL) 
        matchEndPositions <- matchStartPositions + searchSeqLength - 1
        data.table(
            templateBreakpointN = breakpointN,
            templateBreakpointSide = breakpointSide,           
            templateType = templateType,
            templateDistance = as.integer(switch(
                orientation,
                left  = ss - matchEndPositions,
                right = matchStartPositions - 1
            ))
        )
    }
    d <- do.call(rbind, list(
        parseMatches(matches$match1R,   1L, "crossJxn", "retained", "left"),
        parseMatches(matches$match1Rrc, 1L, "foldback", "retained", "left"),
        parseMatches(matches$match1L,   1L, "crossJxn", "lost",     "right"),
        parseMatches(matches$match1Lrc, 1L, "foldback", "lost",     "right"),
        parseMatches(matches$match2R,   2L, "crossJxn", "retained", "right"),
        parseMatches(matches$match2Rrc, 2L, "foldback", "retained", "right"),
        parseMatches(matches$match2L,   2L, "crossJxn", "lost",     "left"),
        parseMatches(matches$match2Lrc, 2L, "foldback", "lost",     "left")
    ))
    
    # for SVs with multiple possible templates, prefer the one closest to the junction
    d[order(templateDistance)][1]
}
svs <- merge(
    svs, 
    if(env$FIND_INSERTION_TEMPLATES) svs[, {
        insertionSize <- -MICROHOM_LEN
        if(insertionSize >= env$MIN_INSERTION_SIZE && insertionSize <= env$MAX_INSERTION_SIZE) analyzeInsertion(PROJECT, .SD)
        else nullInsertion
    }, by = .(PROJECT, SV_ID)] 
    else nullInsertion,
    by = c("PROJECT", "SV_ID"),
    all.x = TRUE
)
#=====================================================================================

#=====================================================================================
message("printing output")
samples[, i := NULL]
results <- list(
    samples = samples,
    targets = targets,
    svs = svs, 
    env = env[unlist(myEnvVars)]
)
saveRDS(results, paste(env$ASSEMBLE_PREFIX, "rds", sep = "."))
#=====================================================================================

# Classes ‘data.table’ and 'data.frame':  34255 obs. of  50 variables:
#  $ SV_ID                     : chr  "10000:1" "10001:1" "10002:1" "10002:10" ...
#  $ TARGET_REGION             : chr  "WWOX" "WWOX" "WWOX" "WWOX" ...
#  $ TARGET_CLASS              : chr  "TT" "TT" "TT" "TT" ...
#  $ JXN_TYPE                  : chr  "I" "I" "I" "I" ...
#  $ SV_SIZE                   : num  1501 6331 1558 3609 1217 ...
#  $ CHROM_1                   : chr  "chr16" "chr16" "chr16" "chr16" ...
#  $ SIDE_1                    : chr  "R" "R" "R" "R" ...
#  $ POS_1                     : int  78494093 78494096 78498346 78498233 78498244 78498393 78499200 78498321 78497952 78496939 ...
#  $ CHROM_2                   : chr  "chr16" "chr16" "chr16" "chr16" ...
#  $ SIDE_2                    : chr  "R" "R" "R" "R" ...
#  $ POS_2                     : int  78495594 78500427 78499904 78501842 78499461 78501864 78501078 78501418 78501773 78499332 ...
    #  $ N_TOTAL                   : int  2 1 1 2 2 2 1 1 4 1 ...
    #  $ N_GAPS                    : int  1 0 0 2 2 1 0 0 1 0 ...
    #  $ N_SPLITS                  : int  1 1 1 0 0 1 1 1 1 1 ...
    #  $ N_OUTER_CLIPS             : int  0 0 0 0 0 0 0 0 2 0 ...
    #  $ N_DUPLEX                  : int  0 0 0 0 0 0 0 0 0 0 ...
    #  $ N_DUPLEX_GS               : int  0 0 0 0 0 0 0 0 0 0 ...
    #  $ STRAND_COUNT              : int  8 9 25 14 3 29 9 6 20 1 ...
    #  $ STRAND_COUNT_GS           : int  8 9 25 14 3 29 9 6 10 1 ...
    #  $ SHARED_PROPER             : num  2 2 2 2 2 2 2 2 0.75 2 ...
    #  $ SHARED_PROPER_GS          : num  2 2 2 2 2 2 2 2 1.5 2 ...
#  $ N_SAMPLES                 : int  2 1 1 1 2 2 1 1 3 1 ...
#  $ SAMPLES                   : chr  "RO3306_Colch_2APH_M,RO3306_G2" "RO3306_Colch_M" "RO3306_G2" "RO3306_Colch_M" ...
#  $ MICROHOM_LEN              : int  5 7 7 0 0 11 4 11 33 10 ...
#  $ JXN_BASES                 : chr  "CACCT" "CCCTTCA" "CTGGCCT" "*" ...

#  $ MAPQ_1                    : chr  "60,60" "60" "60" "60,60" ...
#  $ UMI_1                     : chr  "1,1" "1" "1" "1,1" ...
#  $ MAPQ_2                    : chr  "60,60" "60" "60" "60,60" ...
#  $ UMI_2                     : chr  "1,1" "1" "1" "1,1" ...
#  $ STRAND_COUNT1             : int  0 9 0 0 1 29 9 6 9 1 ...
#  $ STRAND_COUNT2             : int  8 0 25 14 2 0 0 0 11 0 ...
#  $ JUNCTION_NAME             : chr  "16:R:78494093,16:R:78495594" "16:R:78494096,16:R:78500427" "16:R:78498346,16:R:78499904" "16:R:78498233,16:R:78501842" ...
#  $ JUNCTION_NAMES            : chr  "16:R:78494153,16:R:78495957::16:R:78494093,16:R:78495594" "16:R:78494096,16:R:78500427" "16:R:78498346,16:R:78499904" "16:R:78498233,16:R:78501842::16:R:78498244,16:R:78501660" ...
#  $ N_AMBIGUOUS               : int  0 0 0 1 1 0 0 0 1 0 ...
#  $ N_DOWNSAMPLED             : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ N_COLLAPSED               : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ MERGE_LEN                 : int  0 0 0 0 0 0 0 0 0 0 ...
#  $ TARGET_POS_1              : int  2838053 2838056 2842306 2842193 2842204 2842353 2843160 2842281 2841912 2840899 ...
#  $ TARGET_POS_2              : int  2839554 2844387 2843864 2845802 2843421 2845824 2845038 2845378 2845733 2843292 ...
#  $ JXN_SEQ                   : chr  "CCATCCCCCTAGTGCTGTCTTGTGAAAGAGTTCTCCTGAGATCTGTTTGTTTAAAAGCATCTACCACCTCCCTACTTTCTCTCTTCCTCCAGTGCCAGCCATGTGAAGACA"| __truncated__ "TCCTCCAGTGCCAGCCATGTGAAGACATTACCTGCTTCCCCTTCATCTATCTTTGTTCCGCTTTTGTTTGTTTTTTTTAGATCGCTGCATGCTGCCGCGTTGGGTTAACCA"| __truncated__ "GGGGTGGGAAGAGATGTGCCTGGTCTATCAACTGTGAAATCACTGGCCTTGCCTAAATAGGATTAAGTAGTACCAAATGGTCCTTACCGTTCTGAGGCCTCACTGAGACTC"| __truncated__ "*" ...
#  $ GEN_REF_1                 : chr  "CCCTGTACAGCACACAACTCGGGGAATTTTAGGTTTGGTGATTTACTGCTGTTGGAGCCTGGACTTCTTCCTGGGAGAAGGCTCGAAAAACTTCATGTCAGTTACTTAATG"| __truncated__ "TGTACAGCACACAACTCGGGGAATTTTAGGTTTGGTGATTTACTGCTGTTGGAGCCTGGACTTCTTCCTGGGAGAAGGCTCGAAAAACTTCATGTCAGTTACTTAATGGAA"| __truncated__ "ATTGGAGACTTAAACATAAGATCTGACAGTGTAAAATTAGTAGAAGAAAACAAGGAGGAAAGCTCCATGACACTGGTCTGCACAATGAGTGTTTTTGTATATGAACCGAAA"| __truncated__ "GATTCCACAGCACATACTCATTTTTTTGCTCTGATTCTAGCTTTTCCTCCAATTTTATTCAGAATTTGTTCTTTTCTTACCTCACACCATATAGAAAAATCAACTCAAAGT"| __truncated__ ...
#  $ GEN_REF_2                 : chr  "AGATCTTTTTCATTTAAAGTCCAGGAGTCTTCTGTGTGAACCCCACACAGGTACCTTGGAAAGACTGCATCTCCTGGATTTTGCATGGTTCTGGCTTAATGGACAACATAT"| __truncated__ "TGTCTTTTGGCTTGTATCCCACCTGTTCACTTTTTACCATAGCTAACCTTTCTAAGTGGTCAGGTGGGACTCATGGAGATTTTCCGTTCTTTGGCAACTAGCAGCGACTGC"| __truncated__ "ACATAGGCAACTTCCTCGATTGTCTCTCCTTTATTGGCATTGAGAGATAAGCACAATAATGATTGTTAAGATGGAATATTATGAAGGCTCAAGGAATACTTCTGAAGGCTT"| __truncated__ "AGAGCATGTACAACAGCGCCTGTGCTGGGCCTGCTAGGAGTCGATATTTTCTGGCTCAAAGAAGTTCCAAAAGAATTCAAAGCCAGTGTGTCTCTTTCTCTGGACGGCCAA"| __truncated__ ...
#  $ GEN_COV_1                 : chr  "" "" "" "" ...
#  $ GEN_COV_2                 : chr  "" "" "" "" ...
#  $ RO3306_0.2APH_Colch_2APH_M: int  0 0 0 0 1 1 0 1 1 0 ...
#  $ RO3306_0.2APH_Colch_M     : int  0 0 0 0 0 0 0 0 1 0 ...
#  $ RO3306_0.2APH_G2          : int  0 0 0 0 0 1 0 0 0 0 ...
#  $ RO3306_Colch_2APH_M       : int  1 0 0 0 1 0 0 0 0 0 ...
#  $ RO3306_Colch_M            : int  0 1 0 2 0 0 0 0 2 1 ...
#  $ RO3306_G2                 : int  1 0 1 0 0 0 1 0 0 0 ...
