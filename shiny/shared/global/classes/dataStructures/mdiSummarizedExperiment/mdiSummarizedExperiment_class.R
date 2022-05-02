#----------------------------------------------------------------------
# object class for reading and handling data that comes to us in the
# form of a flat table of values where rows are either
#     (overlapping) genome features
#     non-overlapping fixed width bins
#     non-overlapping variable width bins (rowType = bins)
# and columns are either
#     bin/feature attributes (in leading columns), or
#     sample values (in trailing columns)
#----------------------------------------------------------------------
# requires a column name header and a static sample manifest to help
# identify the header format; typically this is source$unique for
# Illumina sequencing data
#----------------------------------------------------------------------
# the input dataObject might be a filename, data.table or data.frame,
# all coerced to a data.table for internal use by the class, or a
# pre-existing compatible RangedSummarizedExperiment with a 'raw' assay
#----------------------------------------------------------------------
# all tables/files must have chromosome/chrom/chr, start and end columns
# typically, an extended BED (half-open) format is expected, but
# you can also set the BED argument to FALSE if start is already 1-referenced
#----------------------------------------------------------------------
# this class returns a reactiveVal object of the Bioconductor S4 class 
# RangedSummarizedExperiment, along with static methods (we do not
# formally extend the RangedSummarizedExperiment class as the MDI
# is not intended to be released into the Bioconductor framework;
# instead, we return the RangedSummarizedExperiment as an rse() reactive)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN CLASS
#----------------------------------------------------------------------
new_mdiSummarizedExperiment <- function(
        dataObject, 
        stringsAsFactors = FALSE, colClasses = NULL, # for when dataObject is a file path (character)
        manifest = NULL, BED = TRUE, stranded = FALSE, # for when dataObject is a file path or data.table; manifest is static data.table # nolint
        rpkm = FALSE # if TRUE
) {
    class <- 'mdiSummarizedExperiment' # for reportProgress tracing
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# load the source data as a data.table or RangedSummarizedExperiment
#----------------------------------------------------------------------
reportProgress('load data', class)
rse <- if(is.character(dataObject)){ # interpret character strings as file paths to a flat file
    fread(
        file = dataObject,
        header = TRUE, # input files must have headers
        stringsAsFactors = stringsAsFactors, 
        colClasses = colClasses,
        blank.lines.skip = TRUE
    )
} else if (is.data.frame(dataObject)) {
    data.table(dataObject)
} else dataObject

#----------------------------------------------------------------------
# parse a data.table to an initial RangedSummarizedExperiment
#----------------------------------------------------------------------
reportProgress('convert to RangedSummarizedExperiment', class)
if (is.data.table(rse)) rse <- dataTableToRSE(rse, manifest, BED, stranded)

#----------------------------------------------------------------------
# enforce data.table format of the SE assays
#----------------------------------------------------------------------
reportProgress('convert assays to data.table', class)
for(i in seq_along(SummarizedExperiment::assays(rse))){   
    if(!is.data.table(SummarizedExperiment::assay(rse, i))){      
        SummarizedExperiment::assay(rse, i) <- data.table(SummarizedExperiment::assay(rse, i))  
    }     
}

#----------------------------------------------------------------------
# enforce a strict, predictable, natural-chromosome order for the data rows
#----------------------------------------------------------------------
reportProgress('sorting by chromosome/position', class)
rse <- GenomeInfoDb::sortSeqlevels(rse) # seqnames is a factor, sort it first
rse <- SummarizedExperiment::sort(rse)  # when then determines sort of rows

#----------------------------------------------------------------------
# calculate some additional bin metadata
#----------------------------------------------------------------------
reportProgress('calculating metadata', class)
rseMetadata <- S4Vectors::metadata(rse)
rowData <- data.table(data.frame(SummarizedExperiment::rowData(rse)))
colData <- SummarizedExperiment::colData(rse)
colData <- data.table(sampleN = seq_len(nrow(colData)))

# bin width format
widths <- GenomicRanges::width(rse)
aggWidths <- aggregate(widths, list(widths), length)
rseMetadata$widthFormat <- if(nrow(aggWidths) == 1) 'fixed' else 'variable'
rseMetadata$width <- if(rseMetadata$widthFormat == 'fixed') aggWidths[1, 1] else NA

# calculate data aggregates
rseMetadata$nFeatures <- nrow(rowData)
rseMetadata$nSamples  <- nrow(colData)
raw <- SummarizedExperiment::assay(rse, 'raw')
rowData[, sum := rowSums(raw, na.rm = TRUE)]
colData[, sum := colSums(raw, na.rm = TRUE)]
rseMetadata$sum <- rowData[, sum(sum, na.rm = TRUE)]
colData[, weight := sum / rseMetadata$sum]

# number the ordered bins over all chromosomes
rowData[, chrom := as.character(GenomicRanges::seqnames(rse))]
rowData[, binN := 1:.N] 
rowData[, chromBinN := 1:.N, by = chrom]
rowData[, isLastBin := chromBinN == .N, by = chrom]

# calculate more data aggregates
if('excluded' %notin% names(rowData)) rowData[, excluded := 0]
rowData[, counted := widths]
if(rpkm){
    rpkm <- function(bp){
        kb <- bp / 1e3
        mapply(
            function(colRaw, colSum) colRaw / kb / (colSum / 1e6),
            raw,
            colData[, sum]
        )
    }
    SummarizedExperiment::assays(rse, withDimnames = FALSE)$rpkm <- data.table( rpkm(rowData[, counted]) )
    rowData[, rpkm := rowMeans(SummarizedExperiment::assays(rse)$rpkm, na.rm = TRUE)]
}

# finalize the rse
S4Vectors::metadata(rse) <- rseMetadata
SummarizedExperiment::rowData(rse) <- rowData
SummarizedExperiment::colData(rse) <- cbind(SummarizedExperiment::colData(rse), as.data.frame(colData))

#----------------------------------------------------------------------
# declare static methods added to an object's list value, called as obj$method()
#----------------------------------------------------------------------
# these are stored with serialized, cached objects and can therefore
# get stale; should be restricted to simple, stable, class-specific actions
#----------------------------------------------------------------------

# setter functions
setMetadata <- function(label, data){
    rse <- obj$rse()
    S4Vectors::metadata(rse)[[label]] <- data
    obj$rse(rse)
}
setMCols <- function(label, data){
    rse <- obj$rse()
    SummarizedExperiment::mcols(rse)[[label]] <- data
    obj$rse(rse)
}
setRowData <- function(label, data){
    rse <- obj$rse()
    SummarizedExperiment::rowData(rse)[[label]] <- data
    obj$rse(rse)
}
setColData <- function(label, data){
    rse <- obj$rse()
    SummarizedExperiment::colData(rse)[[label]] <- data
    obj$rse(rse)
}
setAssay <- function(assay, data){
    rse <- obj$rse()
    x <- SummarizedExperiment::assays(rse)$raw
    x <- x[, lapply(.SD, '*', 0.0)]
    x[] <- data
    SummarizedExperiment::assays(rse, withDimnames = FALSE)[[assay]] <- x
    obj$rse(rse)   
}
setCol <- function(assay, column, data){
    rse <- obj$rse()
    SummarizedExperiment::assays(rse)[[assay]][[column]] <- data
    obj$rse(rse)   
}

# getter functions
getMetadata <- function() S4Vectors::metadata(obj$rse())
getMCols    <- function(rows=TRUE)     SummarizedExperiment::mcols(obj$rse())[rows, ]
getRowData  <- function(rows=TRUE)     SummarizedExperiment::rowData(obj$rse())[rows, ]
getColData  <- function(cols=TRUE)     SummarizedExperiment::colData(obj$rse())[cols,] # yes, cols is i, remember it is transposed # nolint
getAssay <- function(assay, rows=TRUE) SummarizedExperiment::assay(obj$rse(), assay)[rows]
subset <- function(rows=TRUE, cols=TRUE) new_mdiSummarizedExperiment(obj$rse()[rows, cols])    
assayExists <- function(assay){
    if(is.null(assay)) return(FALSE)
    assay %in% SummarizedExperiment::assayNames(obj$rse())
}

#----------------------------------------------------------------------
# set the return value
#----------------------------------------------------------------------
obj <- structure(
    list(
        rse = reactiveVal(rse), # the SummarizedExperiment reactive
        #------------------------
        setMetadata = setMetadata,
        setMCols = setMCols,
        setRowData = setRowData,
        setColData = setColData,
        setAssay = setAssay,
        setCol = setCol,
        #------------------------
        getMetadata = getMetadata,
        getMCols = getMCols,
        getRowData = getRowData,
        getColData = getColData,
        getAssay = getAssay,
        assayExists = assayExists,
        subset = subset
    ),
    class = class
)
rm(rse, rseMetadata, rowData, colData)
obj

#----------------------------------------------------------------------
# END CLASS
#----------------------------------------------------------------------
}
#----------------------------------------------------------------------
