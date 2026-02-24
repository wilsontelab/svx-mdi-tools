# assemble and tabulate SVs across multiple previously executed find actions

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
    library(yaml)
    library(stringi)
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
        'ASSEMBLE_PREFIX',
        'FEATURE_FILES',
        'FEATURE_TYPES'
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
        'COVERAGE_BIN_SIZE',
        'MIN_INSERTION_SIZE',
        'MAX_INSERTION_SIZE',
        'MIN_TEMPLATE_SIZE',
        'INSERTION_SEARCH_SPACE',
        'BASE_USAGE_SPAN',
        'FLEXIBILITY_SPAN',
        'FEATURES_SPAN',
        'MIN_INVERSION_SIZE'
    ),
    double = c(
        'MAX_SHARED_PROPER'
    ),
    logical = c(
        'WARN_ONLY',
        'INCLUDE_CLIPS',
        'FORCE_COVERAGE',
        'FIND_INSERTION_TEMPLATES',
        'PROFILE_BASE_USAGE',
        'PROFILE_FLEXIBILITY',
        'PROFILE_FEATURES'
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
inputDirs <- strsplit(env$INPUT_DIR, ",")[[1]]
getProjectDir <- function(project) { # support multiple input directories to aggregate samples across higher order data sets
    for(inputDir in inputDirs){
        dir <- file.path(inputDir, project)
        if(dir.exists(dir)) return(dir)
    }
    stop(paste("project directory not found:", project))
}
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
getSampleCoverageFile <- function(sample){
    paste(getSampleBamFile(sample), "coverage", "rds", sep = ".")
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
if("skip_assembly" %in% names(samples)) {
    skipped <- samples[, skip_assembly != "-"] 
    message("removing skipped samples from assembly")
    message(paste0("  ", paste(samples[skipped, sample], collapse = ", ")))
    samples <- samples[!skipped]
    samples$skip_assembly <- NULL
}
samples[, i := 1:.N]

message("collecting samples metadata")
samples <- merge(
    samples, 
    samples[, {
        task <- getSampleExtractTask(project, sample)$extract
        .(
            sampleKey     = paste(project, sample, sep = "::"),
            genomesDir    = task$genome[["genomes-dir"]],
            genome        = task$genome$genome,
            cram          = as.logical(task[["bam-format"]][["use-cram"]]),
            minMapq       = as.integer(task[["read-filtering"]][["min-mapq"]]),
            targetsBed    = task[["target-region"]][["targets-bed"]],
            regionPadding = as.integer(task[["target-region"]][["region-padding"]])
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
if(!passedCrosscheck) if(env$WARN_ONLY) message() else stop("crosscheck failed, aborting since --warn-only not set")
#=====================================================================================

#=====================================================================================
message("assembling target region sequences")
refFlats <- list()
getGeneCoordinate <- function(genome, geneNames, column, aggFn){
    if(is.null(refFlats[[genome]])) NA else sapply(geneNames, function(geneName_){
        x <- refFlats[[genome]][geneName == geneName_, unlist(.SD), .SDcols = column]
        if(length(x) > 0) aggFn(x) else NA
    })
}
targets <- samples[, { # one row per target per BED file, can technically have the same target in two rows if duplicated in BED
    targets <- fread(targetsBed)[, 1:4]
    names(targets) <- c("chrom","regionStart","regionEnd","name")
    regionPadding <- max(regionPadding)
    genomeDir   <- Sys.glob(file.path(genomesDir[1], "iGenomes/*/UCSC", genome))
    refFlat     <- file.path(genomeDir, "Annotation/Genes/refFlat.txt.gz")
    genomeFasta <- file.path(genomeDir, "Sequence/WholeGenomeFasta/genome.fa")
    if(is.null(refFlats[[genome]])) refFlats[[genome]] <<- tryCatch({
        x <- fread(refFlat, sep = "\t")
        names(x) <- c("geneName","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds")
        x
    }, error = function(e) NULL)
    targets[, ":="(
        regionKey     = paste0(genome, "::", chrom, ":", regionStart, "-", regionEnd),
        regionPadding = regionPadding,
        paddedStart   = regionStart - regionPadding, # all still half-open like the source BED
        paddedEnd     = regionEnd   + regionPadding,
        geneStart     = getGeneCoordinate(genome, name, "txStart", min),
        geneEnd       = getGeneCoordinate(genome, name, "txEnd",   max),
        geneStrand    = getGeneCoordinate(genome, name, "strand",  function(x) x[1])
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
uniqueTargets <- targets[, # one row per unique target region, defined by regionKey, even if repeated between targets BED files
    {
        targetsBeds <- targetsBed
        .(
            genome = genome[1], 
            chrom = chrom[1], 
            regionStart = regionStart[1], 
            regionEnd = regionEnd[1], 
            name = paste(sort(unique(name)), sep = ","),
            regionPadding = max(regionPadding),
            paddedStart   = min(paddedStart), # all still half-open like the source BED
            paddedEnd     = max(paddedEnd),
            geneStart = geneStart[1], 
            geneEnd = geneEnd[1], 
            geneStrand = geneStrand[1],
            paddedSequence = paddedSequence[which.max(regionPadding)],
            targetsBeds = list(targetsBed),
            sampleKeys = list(samples[targetsBed %in% targetsBeds, sort(sampleKey)])
        )
    },
    by = .(regionKey)
]
#=====================================================================================

#=====================================================================================
message("getting N50 and adjusted target coverage by sample")
samples <- merge(
    samples,
    do.call(rbind, mclapply(1:nrow(samples), function(i){
        sample <- samples[i]
        message(paste("", sample$project, sample$sample, sep = "\t"))
        covFile <- getSampleCoverageFile(sample)
        sampleCoverage <- if(file.exists(covFile) && !env$FORCE_COVERAGE) readRDS(file = covFile) else {
            dt <- fread(text = system2(
                "perl", 
                args = c(
                    file.path(env$ACTION_DIR, "coverage.pl"),
                    sapply(c("genome", "targetsBed", "regionPadding", "minMapq"), getSampleMetadata2, sample),
                    getSampleBamFile(sample)
                ),
                stdout = TRUE, 
                stderr = FALSE
            ))
            saveRDS(dt, file = covFile) # cache slow-step coverage determinations to speed future incremental assemblies
            dt
        }
        cbind(sample[, .(i, sampleKey)], sampleCoverage)
    }, mc.cores = env$N_CPU)),
    by = c("i","sampleKey"),
    all.x = TRUE
) 

message("parsing bin coverage matrices")
binnedCoverage <- sapply(uniqueTargets$regionKey, function(regionKey_) {
    target <- uniqueTargets[regionKey == regionKey_]
    targetSampleKeys <- unlist(target$sampleKeys)
    matrix( # a named list of sample bin coverage matrices by unique target region
        NA_real_, 
        nrow = target[, 1 + floor( (paddedEnd - (paddedStart + 1)) / env$COVERAGE_BIN_SIZE )], # rows = bins for this target region
        ncol = length(targetSampleKeys), # one named matrix column per sample that used this target
        dimnames = list(NULL, targetSampleKeys)
    )
}, simplify = FALSE, USE.NAMES = TRUE)
for(regionKey_ in uniqueTargets$regionKey){
    target <- uniqueTargets[regionKey == regionKey_]
    message(paste("", target$name, regionKey_, sep = "\t"))
    for(targetSampleKey in unlist(target$sampleKeys)){
        sample <- samples[sampleKey == targetSampleKey]
        targetI <- targets[targetsBed == sample$targetsBed, which(regionKey == regionKey_)]
        binnedCoverage[[regionKey_]][, targetSampleKey] <- sample[, {
            as.double(strsplit(strsplit(bin_coverages, "::")[[1]][targetI], ",")[[1]])
        }]
    }
}
samples[, bin_coverages := NULL]
#=====================================================================================

#=====================================================================================
message("tabulating and filtering SVs across all projects and samples")
applySVFilters <- function(project_, svs){

    # reject SVs only found in original 'find' samples that were omitted from the assembly sample list
    # strip the missing samples from the SAMPLES list when one was, one was not present
    svSamples <- svs[, 
        .(SAMPLE = strsplit(SAMPLES, ",")[[1]]),
         by = SV_ID
    ][
        SAMPLE %in% samples[project == project_, unique(sample)]
    ][, 
        .(SAMPLES = paste(SAMPLE, collapse = ",")), 
        by = SV_ID
    ]
    svs <- svs[
        SV_ID %in% svSamples$SV_ID
    ][, 
        SAMPLES := svSamples$SAMPLES
    ][
        # apply SV filters
        (JXN_TYPE == "T" | 
            (SV_SIZE >= env$MIN_SV_SIZE & 
             SV_SIZE <= env$MAX_SV_SIZE)
        ) &
        N_SAMPLES >= env$MIN_SAMPLES &
        N_SAMPLES <= env$MAX_SAMPLES &
        N_SPLITS >= env$MIN_SPLIT_READS & 
        toupper(TARGET_CLASS) %in% c("TT", "TA") & # assemble only ever returns on-target SVs, including inter-target translocations
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
    applySVFilters(project, readRDS(projectFindFile))[, .(
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
        # ,
        # JXN_SEQ ######################
    )]
}, by = .(project)]
names(svs)[1] <- "PROJECT"

message("counting filtered SVs by sample and type")
samples <- merge(
    samples,
    samples[, {
        jxnTypes <- svs[PROJECT == project & grepl(paste0(",", sample, ","), SAMPLES), JXN_TYPE]
        svSizes  <- svs[PROJECT == project & grepl(paste0(",", sample, ","), SAMPLES), SV_SIZE]
        .(
            deletion      = sum(jxnTypes == "L"),
            duplication   = sum(jxnTypes == "D"),
                # inversions apply a distinct size filter due to unique error mechanism
                # NB: only applies to count in samples table; the svs table does NOT apply env$MIN_INVERSION_SIZE (yet)
            inversion     = sum(jxnTypes == "I" & svSizes >= env$MIN_INVERSION_SIZE), 
            translocation = sum(jxnTypes == "T")
        )
    }, by = .(i, sampleKey)],
    by = c("i", "sampleKey"),
    all.x = TRUE
)
#=====================================================================================

#=====================================================================================
if(env$FIND_INSERTION_TEMPLATES) message("attempting to find insertion templates near junctions")
nullInsertion <- data.table( # NA values for junctions where template searching was not performed
    PROJECT = NA_character_,
    SV_ID   = NA_character_,
    templateStartPos = NA_integer_,     # template match start pos, between 1 and 2 * INSERTION_SEARCH_SPACE
    templateStartFm = NA_integer_,      # microhomology length on the left side, i.e., at templateStartPos
    templateStartIsRetained = NA,       # whether templateStartPos is found in the final junction sequence
    templateEndPos = NA_integer_,       # template match end pos, between 1 and 2 * INSERTION_SEARCH_SPACE
    templateEndFm = NA_integer_,        # microhomology length on the right side, i.e., at templateEndPos
    templateEndIsRetained = NA,         # whether templateEndPos is found in the final junction sequence
    templateBreakpointN = NA_integer_,  # the breakpoint side of the final junction sequence, either 1 or 2
    templateIsRc = NA,                  # TRUE if template was the reverse complement of insertion sequence
    templateType = NA_character_,       # the inferred template mechanism: notFound, foldback, strandSwitch, crossJxn, slippage, palindrome, or other
    templateDistance = NA_integer_,     # no. of bp from junction position to innermost microhomology position
    templateInstances = NA_integer_     # no. of distinct instances of found templates, including the one that is reported
)
getSvTarget <- function(project_, sv){
    targetsBed_ <- samples[project == project_ & sample == strsplit(sv$SAMPLES, ",")[[1]][2], targetsBed] # all samples in a project were analyzed with the same targets
    targets[targetsBed == targetsBed_ & name == sv$TARGET_REGION] 
}
isUsablePosition <- function(pos, target, ss) pos - target$paddedStart >  ss && 
                                              target$paddedEnd - pos   >= ss
getBreakpointSequences <- function(sv, target, ss){
    pos1 <- sv$POS_1 - target$paddedStart
    pos2 <- sv$POS_2 - target$paddedStart
    if(sv$SIDE_1 == "L"){
        genRefRetained1 <-    substr(target$paddedSequence, pos1 - ss + 1, pos1)
        genRefLost1     <-    substr(target$paddedSequence, pos1 + 1,      pos1 + ss)
    } else {
        genRefLost1     <- rc(substr(target$paddedSequence, pos1 - ss, pos1 - 1))
        genRefRetained1 <- rc(substr(target$paddedSequence, pos1,      pos1 + ss - 1)) # flip one reference strand for inversions to match jxnSeq assembly
    }
    if(sv$SIDE_2 == "R"){
        genRefLost2     <-    substr(target$paddedSequence, pos2 - ss, pos2 - 1)
        genRefRetained2 <-    substr(target$paddedSequence, pos2,      pos2 + ss - 1)
    } else {
        genRefRetained2 <- rc(substr(target$paddedSequence, pos2 - ss + 1, pos2))
        genRefLost2     <- rc(substr(target$paddedSequence, pos2 + 1,      pos2 + ss))
    } 
    list(
        genRef1         = paste0(genRefRetained1, genRefLost1),
        genRef2         = paste0(genRefLost2, genRefRetained2),
        genRefRetained1 = genRefRetained1,
        genRefLost1     = genRefLost1,
        genRefRetained2 = genRefRetained2,
        genRefLost2     = genRefLost2
    )
}
posIsRetained <- function(breakpointN, pos, ss) switch( # retained positions are found in the final assembled junction sequence
    breakpointN,
    pos <= ss,
    pos >= ss + 1
)
# debug <- data.table( # this code helps debug and explore mechanisms by writing a tmp file of junction sequences
#     templateBreakpointN = integer(),
#     templateStartPos = integer(),
#     templateEndPos = integer(),
#     templateStartFm = integer(),
#     templateEndFm = integer(),
#     templateDistance = integer(),
#     JXN_BASES = character(),
#     JXN_SEQ = character(),
#     genRef1 = character(),
#     genRef2 = character()
# )
analyzeInsertion <- function(svI){
    sv <- svs[svI]
    ss <- env$INSERTION_SEARCH_SPACE

    # parse the genome reference target
    target <- getSvTarget(sv$PROJECT, sv)
    if(!isUsablePosition(sv$POS_1, target, ss) || 
       !isUsablePosition(sv$POS_2, target, ss)) return(nullInsertion)

    # get the genome regions to search as the retained and lost (to this junction) segments on each side of each breakpoint
    seqs <- getBreakpointSequences(sv, target, ss)

    # parse the insertion search sequence as (microhomology)(insertion)(microhomology)
    # dynamically adjust required microhomology spans to ensure search sequences are long enough to be meaningful
    fm <- {
        insertionSize <- -sv$MICROHOM_LEN
        if(insertionSize >= env$MIN_TEMPLATE_SIZE) 1L # demand at least 1 bp flanking microhomology regardless of insertion size
        else as.integer(ceiling((env$MIN_TEMPLATE_SIZE - insertionSize) / 2)) # equally distribute required microhomology to achieve MIN_TEMPLATE_SIZE     
    }
    microhomology1 <- substr(seqs$genRefRetained1, ss - fm + 1, ss) # junction sequences defined by top strand of assembled junction
    microhomology2 <- substr(seqs$genRefRetained2, 1, fm)
    searchSeq <- paste0(microhomology1, sv$JXN_BASES, microhomology2)
    rcSearchSeq <- rc(searchSeq) # accomplish rc search by searchng for rc(template) in top strand of breakpoint sequences
    searchSeqLength <- nchar(searchSeq)

    # find all matches
    matches <- list(
        match1   = gregexpr(  searchSeq, seqs$genRef1)[[1]],
        match1rc = gregexpr(rcSearchSeq, seqs$genRef1)[[1]],
        match2   = gregexpr(  searchSeq, seqs$genRef2)[[1]],
        match2rc = gregexpr(rcSearchSeq, seqs$genRef2)[[1]]
    )
    if(all(unlist(matches) == -1)) return({
        x <- copy(nullInsertion)
        x[, ":="(
            PROJECT = sv$PROJECT,
            SV_ID   = sv$SV_ID,
            templateType = "notFound" # notFound insertions NA except for this templateType
        )]
    })

    # characterize all template matches
    parseMatches <- function(matchStartPositions, seq, breakpointN, isRc){
        if(all(matchStartPositions == -1)) return(NULL) 
        jxnPos <- if(breakpointN == 1) ss else ss + 1
        seqLen <- nchar(seq)

        # if possible, widen microhomology spans to include all matching bases beyond the minimal requirement
        # from here forward, microhomologies can be different lengths on the two sides of the insertion
        cbind(
            do.call(rbind, lapply(matchStartPositions, function(pos){ # adjust start positions
                for(inc in 1:1000){ # far longer than any microhomology can ever be, loop will always break
                    pos_ <- pos - inc
                    fm_  <- fm  + inc
                    if(pos_ < 1) break
                    candidate <- substr(seq, pos_, pos_ + fm_ - 1)
                    microhomology <- if(isRc){
                        rc(substr(seqs$genRefRetained2, 1, fm_))
                    } else {
                        substr(seqs$genRefRetained1, ss - fm_ + 1, ss)
                    }
                    if(candidate != microhomology) break
                }
                data.table(
                    templateStartPos = pos_ + 1L,
                    templateStartFm  = fm_  - 1L,
                    templateStartIsRetained = posIsRetained(breakpointN, pos_ + 1L, ss)
                )
            })),
            do.call(rbind, lapply(matchStartPositions + searchSeqLength - 1L, function(pos){ # adjust end positions
                for(inc in 1:1000){
                    pos_ <- pos + inc
                    fm_  <- fm  + inc
                    if(pos_ > seqLen) break
                    candidate <- substr(seq, pos_ - fm_ + 1, pos_)
                    microhomology <- if(isRc){
                        rc(substr(seqs$genRefRetained1, ss - fm_ + 1, ss))
                    } else {
                        substr(seqs$genRefRetained2, 1, fm_)
                    }
                    if(candidate != microhomology) break
                }
                data.table(
                    templateEndPos = pos_ - 1L,
                    templateEndFm  = fm_  - 1L,
                    templateEndIsRetained = posIsRetained(breakpointN, pos_ - 1L, ss)
                )
            }))
        )[, ":="( 
            PROJECT = sv$PROJECT,
            SV_ID   = sv$SV_ID,
            templateBreakpointN = breakpointN,
            templateIsRc = isRc,
            templateType = mapply(function(startIsRetained, endIsRetained) {
                if(isRc){
                    if(startIsRetained & endIsRetained) {
                        "foldback"
                    } else {
                        "strandSwitch"
                    }
                } else {
                    if(startIsRetained & endIsRetained) {
                        "crossJxn"
                    } else if(startIsRetained != endIsRetained) {
                        "slippage"
                    } else {
                        "other"
                    }
                }
            }, templateStartIsRetained, templateEndIsRetained)
        )][,
            templateDistance := mapply(function(startIsRetained, endIsRetained, startPos, endPos) as.integer(
                if(startIsRetained & endIsRetained){ # crossJxn or foldback
                    if(breakpointN == 1) jxnPos - endPos # can be zero if all the way to the retained end
                    else startPos - jxnPos
                } else if(!startIsRetained & !endIsRetained){ # some strandSwitch, all other
                    if(breakpointN == 1) startPos - (jxnPos + 1) # can be a negative number
                    else (jxnPos - 1) - endPos
                } else { # some strandSwitch, all slippage
                    if(breakpointN == 1) endPos - (jxnPos + 1)
                    else (jxnPos - 1) - startPos
                }
            ), templateStartIsRetained, templateEndIsRetained, templateStartPos, templateEndPos)
        ]
    }
    d <- do.call(rbind, list(
        parseMatches(matches$match1,   seqs$genRef1, 1L, FALSE),
        parseMatches(matches$match1rc, seqs$genRef1, 1L, TRUE),
        parseMatches(matches$match2,   seqs$genRef2, 2L, FALSE),
        parseMatches(matches$match2rc, seqs$genRef2, 2L, TRUE)
    ))

    # for SVs with multiple possible templates, report only one, preferring ...
    d <- d[, 
        templateInstances := .N
    ][order(
        -(templateEndPos - templateStartPos + 1), # ... the longest template span, including microhomologies and the insertion template
        templateDistance                          # ... the one closest to the junction
    )][1]
    d[
        templateType == "foldback" && templateDistance == 0,
        templateType := "palindrome" # create a distinct type for this special case of apparent foldbacks
    ]
    # if(sv$JXN_TYPE == "L" && 
    #    d$templateType == "foldback" && 
    #      d$templateDistance == 0) debug <<- rbind(debug, d[, .(
    #     templateBreakpointN,
    #     templateStartPos,
    #     templateEndPos,
    #     templateStartFm,
    #     templateEndFm,
    #     templateDistance,
    #     JXN_BASES = sv$JXN_BASES,
    #     JXN_SEQ = sv$JXN_SEQ,
    #     genRef1 = seqs$genRef1,
    #     genRef2 = seqs$genRef2
    # )])
    d
}
svs <- merge(
    svs, 
    if(!env$FIND_INSERTION_TEMPLATES) nullInsertion else {
        do.call(rbind, mclapply(
        # do.call(rbind, lapply( # for debug
            which(svs$N_SPLITS > 0 & 
                 -svs$MICROHOM_LEN >= env$MIN_INSERTION_SIZE & 
                 -svs$MICROHOM_LEN <= env$MAX_INSERTION_SIZE & 
                 svs$JXN_TYPE != "T"), 
            analyzeInsertion,
            mc.cores = env$N_CPU
        ))
    },
    by = c("PROJECT", "SV_ID"),
    all.x = TRUE
)
# write.table(
#     debug, 
#     file = "/nfs/turbo/path-wilsonte-turbo/mdi/wilsontelab/greatlakes/data-scripts/wilsontelab/publications/CFS-M_phase-PolQ-2023/assemble_all_projects/DEBUG.csv", 
#     sep = ",",
#     row.names = FALSE,
#     col.names = TRUE
# )
# stop("XXXXXXXXXXXXXXXXXXXXX")
#=====================================================================================

#=====================================================================================
if(env$PROFILE_BASE_USAGE) message("profiling base usage")
bases <- list(A = 1L, C = 2L, G = 3L, T = 4L, N = NA_integer_)
getMicrohomologyType <- function(sv){
    if(sv$MICROHOM_LEN > 0) "microhomology" 
    else if(sv$MICROHOM_LEN < 0) "insertion" 
    else "blunt"
}
analyzeBaseUsage <- function(svI){
    sv <- svs[svI]
    ss <- env$BASE_USAGE_SPAN

    # parse the genome reference target
    target <- getSvTarget(sv$PROJECT, sv)
    if(!isUsablePosition(sv$POS_1, target, ss) || 
       !isUsablePosition(sv$POS_2, target, ss)) return(NULL)

    # get the genome regions to search as the retained and lost (to this junction) segments on each side of each breakpoint
    seqs <- getBreakpointSequences(sv, target, ss)
    
    # profile the base usage around each breakpoint
    # orient each breakpoint in order retained -> lost, tabulating the bases on the "top" strand when approaching jxnPos from the left
    matrix(
        c(
            unlist(bases[strsplit(   seqs$genRef1,  "")[[1]]]),
            unlist(bases[strsplit(rc(seqs$genRef2), "")[[1]]])
        ), 
        byrow = TRUE, # one row per junction breakpoint, one col for each position in the profile span; jxnPos at env$BASE_USAGE_SPAN
        nrow = 2,
        dimnames = list(paste(sv$PROJECT, sv$SV_ID, 1:2, sep = "::"), NULL)
    )
}
baseUsageProfile <- if(!env$PROFILE_BASE_USAGE) NULL else {
    do.call(rbind, mclapply(
        which(svs$N_SPLITS > 0 & svs$JXN_TYPE != "T"), # sequenced junctions with known/usable junction positions
        analyzeBaseUsage,
        mc.cores = env$N_CPU
    ))
}
junctionBaseUsageProfile <- {
    x <- svs[
        abs(MICROHOM_LEN) > 0 & JXN_TYPE != "T", 
        .(
            microhomologyType = getMicrohomologyType(.SD),
            base = strsplit(JXN_BASES, "")[[1]]
        ), 
        by = .(PROJECT, SV_ID)
    ]
    dcast(
        x,
        "PROJECT + SV_ID + microhomologyType ~ base", 
        value.var = "base",
        fun.aggregate = length, 
        fill = 0
    )[, .SD, 
        .SDcols = c("PROJECT", "SV_ID", "microhomologyType", "A", "C", "G", "T")
    ][
        A + C + G + T > 0 # in case junctions bases were exclusively N
    ]
}
targetBaseUsageProfile <- {
    x <- uniqueTargets[, .(
        base = strsplit(paddedSequence, "")[[1]] # TODO: also profile just the capture target bases?
    ), by = .(regionKey)]
    dcast(
        x,
        "regionKey ~ base", 
        value.var = "base",
        fun.aggregate = length, 
        fill = 0
    )[, .SD, 
        .SDcols = c("regionKey", "A", "C", "G", "T")
    ]
}
#=====================================================================================

#=====================================================================================
if(env$PROFILE_FLEXIBILITY) message("calculating base flexibility profiles")
getBreakpointPositions <- function(sv, target, ss){
    pos1 <- sv$POS_1 - target$paddedStart
    pos2 <- sv$POS_2 - target$paddedStart
    if(sv$SIDE_1 == "L"){
        genRefRetained1 <-    (pos1 - ss + 1):pos1
        genRefLost1     <-    (pos1 + 1):(pos1 + ss)
    } else {
        genRefLost1     <- rev((pos1 - ss):(pos1 - 1))
        genRefRetained1 <- rev(pos1:(pos1 + ss - 1)) # flip one reference strand for inversions to match jxnSeq assembly
    }
    if(sv$SIDE_2 == "R"){
        genRefLost2     <-    (pos2 - ss):(pos2 - 1)
        genRefRetained2 <-    pos2:(pos2 + ss - 1)
    } else {
        genRefRetained2 <- rev((pos2 - ss + 1):pos2)
        genRefLost2     <- rev((pos2 + 1):(pos2 + ss))
    } 
    list(
        genRef1 = c(genRefRetained1, genRefLost1),
        genRef2 = c(genRefLost2, genRefRetained2)
    )
}
analyzeFlexibility <- function(svI){
    sv <- svs[svI]
    ss <- env$FLEXIBILITY_SPAN

    # parse the genome reference target
    target <- getSvTarget(sv$PROJECT, sv)
    if(!isUsablePosition(sv$POS_1, target, ss) || 
       !isUsablePosition(sv$POS_2, target, ss)) return(NULL)
    flexibility <- uniqueTargets[regionKey == target$regionKey, unlist(flexibility)]

    # get the genome regions to search as the retained and lost (to this junction) segments on each side of each breakpoint
    positions <- getBreakpointPositions(sv, target, ss)
    
    # profile the base usage around each breakpoint
    # orient each breakpoint in order retained -> lost, tabulating the bases on the "top" strand when approaching jxnPos from the left
    matrix(
        c(
                flexibility[positions$genRef1],
            rev(flexibility[positions$genRef2])
        ), 
        byrow = TRUE, # one row per junction breakpoint, one col for each position in the profile span; jxnPos at env$BASE_USAGE_SPAN
        nrow = 2,
        dimnames = list(paste(sv$PROJECT, sv$SV_ID, 1:2, sep = "::"), NULL)
    )
}
dinucleotideAngles <- list( # after http://margalit.huji.ac.il/TwistFlex/
    A = list(
        A = 7.6,
        C = 14.6, # A to C dinucleotide, etc.
        G = 8.2,
        T = 25, # A to T dinucleotide has the highest flexibility
        N = NA_real_
    ),
    C = list(
        A = 10.9,
        C = 7.2,
        G = 8.9,
        T = 8.2,
        N = NA_real_
    ),
    G = list(
        A = 8.8,
        C = 11.1,
        G = 7.2,
        T = 14.6,
        N = NA_real_
    ),
    T = list(
        A = 12.5,
        C = 8.8,
        G = 10.9,
        T = 7.6,
        N = NA_real_
    ),
    N = list(
        A = NA_real_,
        C = NA_real_,
        G = NA_real_,
        T = NA_real_,
        N = NA_real_
    )
)
flexibilityProfile <- if(!env$PROFILE_FLEXIBILITY) NULL else {
    message(paste("", "by target", sep = "\t"))
    uniqueTargets[, flexibility := {
        bases <- strsplit(paddedSequence, "")[[1]]
        nBases <- length(bases)
        angles <- sapply(1:(nBases - 1), function(i) dinucleotideAngles[[bases[i]]][[bases[i + 1]]])
        dt <- data.table(
            enteringAngles = c(NA_real_, angles),
            exitingAngles  = c(angles, NA_real_)
        )
        angles <- apply(dt, 1, mean, na.rm = TRUE)
        angles[is.nan(angles)] <- NA_real_
        list(list(angles)) # report a base's flexibility as the average of the bonds entering and leaving it
    }, by = .(regionKey)]
    message(paste("", "by SV", sep = "\t"))
    do.call(rbind, mclapply(
        which(svs$N_SPLITS > 0 & svs$JXN_TYPE != "T"), # sequenced junctions with known/usable junction positions
        analyzeFlexibility,
        mc.cores = env$N_CPU
    ))
}
#=====================================================================================

#=====================================================================================
if(env$PROFILE_FEATURES) message("calculating base matches to genome features")
getBreakpointFeatures <- function(targetPositions, target, features){
    coord1 <- targetPositions + target$paddedStart
    ends1 <- range(coord1)
    matchingFeatures <- features[
        chrom_ == target$chrom &
        ends1[1] <= end_ &
        (start_ + 1) <= ends1[2],  # STRAND?
        .(
            (start_ + 1) - target$paddedStart,
            end_ - target$paddedStart
        )
    ]
    if(nrow(matchingFeatures) == 0) return(rep(0, length(targetPositions)))
    featurePosInTarget <- unlist( apply(matchingFeatures, 1, function(v) v[1]:v[2], simplify = FALSE) )
    x <- merge(
        data.table(
            targetPosition = targetPositions
        ),
        data.table(
            targetPosition = unique(featurePosInTarget),
            hit = TRUE
        ),
        by = "targetPosition",
        all.x = TRUE
    )
    x[is.na(x)] <- FALSE
    if(targetPositions[1] < targetPositions[2]) x[, hit] # return a logical set to TRUE for the positions in breakpoint that match any feature
    else x[, rev(hit)] # since merge forced asending order, we must reverse back when positions were descending
}
analyzeFeatures <- function(svI, features){
    sv <- svs[svI]
    ss <- env$FEATURES_SPAN

    # parse the genome reference target
    target <- getSvTarget(sv$PROJECT, sv)
    if(!isUsablePosition(sv$POS_1, target, ss) || 
       !isUsablePosition(sv$POS_2, target, ss)) return(NULL)

    # get the genome regions to search as the retained and lost (to this junction) segments on each side of each breakpoint
    targetPositions <- getBreakpointPositions(sv, target, ss) # i.e., 1 is the first position in the target region

    # profile the base usage around each breakpoint
    # orient each breakpoint in order retained -> lost, tabulating the bases on the "top" strand when approaching jxnPos from the left
    matrix(
        c(
                getBreakpointFeatures(targetPositions$genRef1, target, features),
            rev(getBreakpointFeatures(targetPositions$genRef2, target, features))
        ), 
        byrow = TRUE, # one row per junction breakpoint, one col for each position in the profile span; jxnPos at env$BASE_USAGE_SPAN
        nrow = 2,
        dimnames = list(paste(sv$PROJECT, sv$SV_ID, 1:2, sep = "::"), NULL)
    )
}
genomeFeatures <- if(!env$PROFILE_FEATURES) NULL else {
    featureFiles <- strsplit(env$FEATURE_FILES, ",")[[1]]
    featureTypes <- strsplit(env$FEATURE_TYPES, ",")[[1]]
    if(length(featureFiles) != length(featureTypes)) stop("malformed request at --feature-files or --feature-types")
    x <- sapply(
        featureFiles,
        function(featureFile){
            message(paste("", basename(featureFile), sep = "\t"))
            if(!file.exists(featureFile)) stop(paste("feature file not found:", featureFile))
            features <- fread(featureFile, sep = "\t")[, ":="(
                chrom_ = .SD[[1]],
                start_ = .SD[[2]],
                end_   = .SD[[3]]
            )]
            features <- uniqueTargets[, {
                features[
                    chrom_ == chrom &
                    (start_ + 1) <= paddedEnd &
                    (paddedStart + 1) <= end_
                ]
            }, by = .(regionKey)]
            do.call(rbind, mclapply(
                which(svs$N_SPLITS > 0 & svs$JXN_TYPE != "T"), # sequenced junctions with known/usable junction positions
                analyzeFeatures,
                features,
                mc.cores = env$N_CPU
            ))
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )
    names(x) <- featureTypes
    x
}
#=====================================================================================

#=====================================================================================
message("printing output")
samples[, i := NULL]
results <- list(
    samples = samples,
    targets = targets,
    uniqueTargets = uniqueTargets,
    binnedCoverage = binnedCoverage,
    svs = svs, 
    baseUsageProfile = baseUsageProfile,
    junctionBaseUsageProfile = junctionBaseUsageProfile,
    targetBaseUsageProfile = targetBaseUsageProfile,
    flexibilityProfile = flexibilityProfile,
    genomeFeatures = genomeFeatures,
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

# table refFlat
#   "A gene prediction with additional geneName field." 
#       (
#       string  geneName;           "Name of gene as it appears in Genome Browser." 
#       string  name;               "Name of gene" 
#       string  chrom;              "Chromosome name" 
#       char[1] strand;             "+ or - for strand" 
#       uint    txStart;            "Transcription start position" 
#       uint    txEnd;              "Transcription end position" 
#       uint    cdsStart;           "Coding region start" 
#       uint    cdsEnd;             "Coding region end" 
#       uint    exonCount;          "Number of exons" 
#       uint[exonCount] exonStarts; "Exon start positions" 
#       uint[exonCount] exonEnds;   "Exon end positions" 
#       )
