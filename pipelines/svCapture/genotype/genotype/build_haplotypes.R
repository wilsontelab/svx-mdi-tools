
# compare structural variants between samples and annotate SV calls with any matching samples

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
# load packages
library(data.table)
library(vcfR)
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
        'MIN_VARIANT_QUAL'
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
setCanonicalChroms()
loadTargetRegions()
loadFaidx()
map <- unique(do.call(rbind, lapply(seq_len(nrow(targetRegions$bed)), function(i){
    range <- targetRegions$bed[i, c(start - env$REGION_PADDING, end + env$REGION_PADDING)]
    data.table(
        chrom = targetRegions$bed[i, chrom],
        pos   = range[1]:range[2],
        ref   = targetRegions$bed[i, strsplit(getRefSeq(chrom, range), "")[[1]]]
    )
})))
map[, ":="(
    hap1 = ref, 
    hap2 = ref,
    posKey = paste(chrom, pos, sep = ":")
)]
setkey(map, posKey)
#=====================================================================================

#=====================================================================================
# load and structure called variants in the source VCF
#-------------------------------------------------------------------------------------
vcf <- vcfR2tidy( read.vcfR(env$CONSTITUTIVE_VCF, verbose = FALSE) )
vcf <- data.table(
    index    = 1:nrow(vcf$fix),
    chrom    = vcf$fix$CHROM,
    pos      = vcf$fix$POS,
    qual     = vcf$fix$QUAL,      
    ref      = vcf$fix$REF, 
    alt      = vcf$fix$ALT, 
    genotype = vcf$gt$gt_GT
)
#=====================================================================================

#=====================================================================================
# loop all variants in vcf and fill into the unphased haplotype base map
#-------------------------------------------------------------------------------------
commitVariant <- function(chrom, pos, ref, alt, alleleI, genotype){
    nRefBases <- nchar(ref)
    nAltBases <- nchar(alt)
    altBases <- strsplit(alt, "")[[1]]
    bases <- if(nRefBases == nAltBases) c( # a SNP
        altBases
    ) else if(nRefBases > nAltBases) c( # a deletion
        altBases, 
        rep("-", nRefBases - nAltBases)
    ) else c( # an insertion
        if(nRefBases > 1) altBases[1:(nRefBases - 1)] else character(),
        paste(altBases[nRefBases:nAltBases], collapse = "")
    )
    posKeys <- paste(chrom, pos:(pos + nRefBases - 1), sep - ":")
    switch(
        genotype,
        "0/1" = {
            if(alleleI == 1) map[posKeys, ":="(
                hap1 = bases
            )]
        },
        "1/1" = {
            if(alleleI == 1) map[posKeys, ":="(
                hap1 = bases,
                hap2 = bases
            )]
        },
        "1/2" = {
            if(alleleI == 1) map[posKeys, ":="(
                hap1 = bases
            )] else if(alleleI == 2) map[posKeys, ":="(
                hap2 = bases
            )] 
        },
        "2/1" = {
            if(alleleI == 1) map[posKeys, ":="(
                hap1 = bases
            )] else if(alleleI == 2) map[posKeys, ":="(
                hap2 = bases
            )] 
        },
        "2/2" = {
            if(alleleI == 2) map[posKeys, ":="(
                hap1 = bases,
                hap2 = bases
            )]
        }
    )
}
vcf[qual >= env$MIN_VARIANT_QUAL][, {
    alt <- strsplit(alt, ",")[[1]]
    commitVariant(chrom, pos, ref, alt[1], 1, genotype)
    if(!is.null(alt[2])) {
        commitVariant(chrom, pos, ref, alt[2], 2, genotype)
    }
}, by = index]
#=====================================================================================

#=====================================================================================
# print results
#-------------------------------------------------------------------------------------
message("writing SV summary table")
outFile <- paste(env$GENOTYPE_PREFIX, 'unphased_haplotypes', 'gz', sep = ".")
fwrite(
    map[, .SD, .SDcols = c("chrom", "pos", "ref", "hap1", "hap2")], 
    file = outFile, 
    quote = FALSE, 
    sep = "\t",    
    row.names = FALSE,   
    col.names = TRUE, 
    compress = "gzip"
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
