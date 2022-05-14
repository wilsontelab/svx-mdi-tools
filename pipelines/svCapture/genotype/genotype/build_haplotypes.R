# make a bi-allelic, unphased haplotype map of all padded target region bases

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(parallel)
library(data.table)
library(vcfR) # for easy reading of VCF files
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
        'GENOTYPE_PREFIX',
        'SHM_DIR_WRK',
        'CONSTITUTIVE_VCF',
        'GENOME',
        'GENOME_CHROMS',
        'GENOME_FASTA',
        'CHROM_FASTA_DIR'        
    ),
    integer = c(
        'N_CPU',
        'REGION_PADDING',
        'MIN_VARIANT_QUAL',
        'MIN_CALL_DEPTH'
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
sourceScripts(file.path(rUtilDir, 'genome'), c('general', 'chroms', 'targets', 'faidx'))
#=====================================================================================

#=====================================================================================
# load and initialize reference genome sequence across all padded target positions
#-------------------------------------------------------------------------------------
message("loading target regions")
setCanonicalChroms()
loadTargetRegions()
message("loading genome reference sequence")
loadFaidx()
message("filling allele map with reference bases")
map <- unique(do.call(rbind, lapply(seq_len(nrow(targetRegions$bed)), function(i){
    range <- targetRegions$bed[i, c(start - env$REGION_PADDING, end + env$REGION_PADDING)]
    data.table(
        CHROM = targetRegions$bed[i, chrom],
        POS   = range[1]:range[2],
        REF   = targetRegions$bed[i, strsplit(getRefSeq(chrom, range), "")[[1]]]
    )
})))
map[, ":="(
    HAP1 = REF, 
    HAP2 = REF,
    posKey = paste(CHROM, POS, sep = ":")
)]
#=====================================================================================

#=====================================================================================
# initial processing of VCF
#-------------------------------------------------------------------------------------

# load and structure called variants in the source VCF
message("loading VCF file")
vcf <- suppressWarnings( vcfR2tidy( read.vcfR(env$CONSTITUTIVE_VCF, verbose = FALSE) ) )

message("parsing VCF file")
vcf <- cbind(
    as.data.table(vcf$fix[, c(
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "QUAL",
        "DP"
    )]),
    as.data.table(vcf$gt[, c(
        "gt_GT_alleles"
    ), drop = FALSE])
)

# record informativity in the base map
# uninformative bases are in low-coverage, typically adjacent, regions
message("recording base informativity in allele map")
map <- merge(
    map, 
    vcf[, .(CHROM = CHROM, POS = POS, INF = as.integer(DP >= env$MIN_CALL_DEPTH))],
    by = c("CHROM", "POS"), 
    all.x = TRUE
)
map[is.na(INF), INF := 0]

# filter the VCF to high-quality, informative variants
message("filtering low quality allelic variants")
vcf <- vcf[
    !is.na(ALT) & 
    QUAL >= env$MIN_VARIANT_QUAL & 
      DP >= env$MIN_CALL_DEPTH,
    .(CHROM, POS, REF, gt_GT_alleles)
]

message("parsing allele values")
ALT <- transpose(as.data.table(strsplit(vcf$gt_GT_alleles, "/")))
setnames(ALT, c("HAP1", "HAP2"))
vcf <- cbind(vcf, ALT)
vcf[, index := 1:.N]
#=====================================================================================

#=====================================================================================
# loop all variants in vcf and fill into the unphased haplotype base map
#-------------------------------------------------------------------------------------
message("filling allele values into map")
setkey(map, posKey)
for(alleleI in 1:2){
    hapCol <- paste0("HAP", alleleI)
    message(paste("  ", hapCol))
    x <- do.call(rbind, mclapply(vcf$index, function(varI){
        alt <- vcf[varI][[hapCol]]
        nAltBases <- nchar(alt)
        altBases <- strsplit(alt, "")[[1]]
        nRefBases <- vcf[varI, nchar(REF)]
        vcf[varI, .(
            posKey = paste(CHROM, POS:(POS + nRefBases - 1), sep = ":"),
            base = if(nRefBases == nAltBases) c( # a SNP
                altBases
            ) else if(nRefBases > nAltBases) c( # a deletion
                altBases, 
                rep("-", nRefBases - nAltBases)
            ) else c( # an insertion
                if(nRefBases > 1) altBases[1:(nRefBases - 1)] else character(),
                paste(altBases[nRefBases:nAltBases], collapse = "")
            )      
        )]
    }, mc.cores = env$N_CPU))
    map[x$posKey, hapCol] <- x$base
}
#=====================================================================================

#=====================================================================================
# print results
#-------------------------------------------------------------------------------------
message("writing SV summary table")
map <- map[, .SD, .SDcols = c("CHROM", "POS", "REF", "HAP1", "HAP2", "INF")]
outFile <- paste(env$GENOTYPE_PREFIX, 'unphased_haplotypes', 'gz', sep = ".")
fwrite(
    map, 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = TRUE, 
    compress = "gzip"
)
outFile <- paste(env$GENOTYPE_PREFIX, 'unphased_haplotypes', 'rds', sep = ".")
saveRDS(
    map, 
    file = outFile
)
#=====================================================================================

# List of 3
#  $ fix : tibble [2,436 × 25] (S3: tbl_df/tbl/data.frame)
#   ..$ ChromKey: int [1:2436] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ CHROM   : chr [1:2436] "chr1" "chr1" "chr1" "chr1" ...
#   ..$ POS     : int [1:2436] 71186045 71189969 71190798 71194437 71205625 71211410 71213252 71213321 71215694 71216543 ...
#   ..$ ID      : chr [1:2436] NA NA NA NA ...
#   ..$ REF     : chr [1:2436] "A" "T" "A" "G" ...
#   ..$ ALT     : chr [1:2436] "ATCAC" "C" "G" "A" ...
#   ..$ QUAL    : num [1:2436] 80.4 9 9 9 9 ...
#   ..$ FILTER  : chr [1:2436] NA NA NA NA ...
#   ..$ INDEL   : logi [1:2436] TRUE FALSE FALSE FALSE FALSE FALSE ...
#   ..$ IDV     : int [1:2436] 2 NA NA NA NA NA NA NA NA NA ...
#   ..$ IMF     : num [1:2436] 1 NA NA NA NA NA NA NA NA NA ...
#   ..$ DP      : int [1:2436] 2 1 1 1 1 1 1 2 1 1 ...
#   ..$ VDB     : num [1:2436] 0.14 NA NA NA NA NA NA 0.66 NA NA ...
#   ..$ RPBZ    : num [1:2436] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ MQBZ    : num [1:2436] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ BQBZ    : num [1:2436] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ MQSBZ   : num [1:2436] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ SCBZ    : num [1:2436] NA NA NA NA NA NA NA NA NA NA ...
#   ..$ FS      : num [1:2436] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ SGB     : num [1:2436] -0.454 -0.38 -0.38 -0.38 -0.38 ...
#   ..$ MQ0F    : num [1:2436] 0 0 0 0 0 0 0 0 0 0 ...
#   ..$ AC      : chr [1:2436] "2" "2" "2" "2" ...
#   ..$ AN      : int [1:2436] 2 2 2 2 2 2 2 2 2 2 ...
#   ..$ DP4     : chr [1:2436] "0,0,2,0" "0,0,1,0" "0,0,1,0" "0,0,1,0" ...
#   ..$ MQ      : int [1:2436] 60 60 60 60 60 60 60 60 60 36 ...
#  $ gt  : tibble [2,436 × 6] (S3: tbl_df/tbl/data.frame)
#   ..$ ChromKey     : int [1:2436] 1 1 1 1 1 1 1 1 1 1 ...
#   ..$ POS          : int [1:2436] 71186045 71189969 71190798 71194437 71205625 71211410 71213252 71213321 71215694 71216543 ...
#   ..$ Indiv        : chr [1:2436] "-" "-" "-" "-" ...
#   ..$ gt_PL        : chr [1:2436] "110,6,0" "38,3,0" "38,3,0" "38,3,0" ...
#   ..$ gt_GT        : chr [1:2436] "1/1" "1/1" "1/1" "1/1" ...
#   ..$ gt_GT_alleles: chr [1:2436] "ATCAC/ATCAC" "C/C" "G/G" "A/A" ...

# Classes ‘data.table’ and 'data.frame':  3750004 obs. of  7 variables:
#  $ CHROM : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ POS   : int  71184317 71184318 71184319 71184320 71184321 71184322 71184323 71184324 71184325 71184326 ...
#  $ REF   : chr  "T" "T" "A" "A" ...
#  $ ALT1  : logi  NA NA NA NA NA NA ...
#  $ ALT2  : logi  NA NA NA NA NA NA ...
#  $ posKey: chr  "chr1:71184317" "chr1:71184318" "chr1:71184319" "chr1:71184320" ...
#  $ INF   : int  NA NA NA NA NA NA NA NA NA NA ...
#  - attr(*, ".internal.selfref")=<externalptr> 
#  - attr(*, "sorted")= chr [1:2] "CHROM" "POS"

# Classes ‘data.table’ and 'data.frame':  1261 obs. of  6 variables:
#  $ CHROM: chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ POS  : int  71234437 71366195 71371514 71371562 71474517 71616708 71684185 71684369 71684600 71684641 ...
#  $ REF  : chr  "TACAC" "T" "G" "T" ...
#  $ ALT  : chr  "T,TAC" "TTGTG" "A" "C" ...
#  $ gt_GT: chr  "1/2" "1/1" "1/1" "1/1" ...
#  $ index: int  1 2 3 4 5 6 7 8 9 10 ...
