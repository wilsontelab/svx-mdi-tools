#----------------------------------------------------------------------
# helper functions for managing genomic data matrices
#----------------------------------------------------------------------

# make sure a data table has chrom/start/end columns
renameChromStartEnd <- function(dt){
    renameColumn <- function(new, old) for(x in old) names(dt)[which(names(dt) == x)] <<- new
    renameColumn('chrom', c('chr', 'Chr', 'CHR', 'chromosome', 'Chromosome', 'CHROMOSOME', 'Chrom', 'CHROM'))
    renameColumn('start', c('Start', 'START'))
    renameColumn('end',   c('End',   'END'))
    dt
}

# make sure a data table is sorted chr1,chr2,...,chr10,...,chrX,...,chr_xyz
naturalChromSortDataTable <- function(dt){
    chroms_ <- sub('CHR', '', dt[, toupper(chrom)])
    dt[order(
        grepl('_', chroms_),
        suppressWarnings(as.integer(chroms_)),
        chroms_,
        start,
        end
    )]    
}

#----------------------------------------------------------------------
# helper functions called by mdiSummarizedExperiment generic methods
# not generally intended to be called in other code, but they can be
#----------------------------------------------------------------------

# the name of assay #1 (should always be 'raw', but just in case)
getBaseAssayName_mse <- function(mse){
    SummarizedExperiment::assayNames(mse$rse())[1]
}

# get a single column of values from the mse's row data or mcols as a vector
getRowDataCol_mse <- function(mse, col, rows=TRUE){
    switch(col,
        chrom = as.character(GenomicRanges::seqnames(mse$rse())),
        start = GenomicRanges::start(mse$rse()),
        end   = GenomicRanges::end(mse$rse()),        
        width = GenomicRanges::width(mse$rse()),
        if(isRowDataCol_mse(mse, col)) SummarizedExperiment::rowData(mse$rse())[[col]]
        else SummarizedExperiment::mcols(mse$rse())[[col]]
    )[rows]
}

# get a single column of one assay's values as a vector
getValueCol_mse <- function(mse, col, assay=NULL, rows=TRUE){
    if(is.null(assay)) assay <- getBaseAssayName_mse(mse)
    SummarizedExperiment::assay(mse$rse(), assay)[[col]][rows]
}

# get a data.table of all or just a subset of the available assays for a given sample
getValuesCols_mse <- function(mse, col, assays=NULL, rows=TRUE){
    if(is.null(assays)) assays <- SummarizedExperiment::assayNames(mse$rse())
    as.data.table( sapply(assays, function(assay) getValueCol_mse(mse, col, assay)) )[rows]
}

# explore the source of a column that was requested by name
rowDataColNames_mse <- function(mse){
    names(SummarizedExperiment::rowData(mse$rse()))
}
sampleNames_mse <- function(mse){
    rownames(SummarizedExperiment::colData(mse$rse()))
}
isRowDataCol_mse <- function(mse, col){
    col %in% c(rowDataColNames_mse(mse), 'chrom', 'start', 'end', 'width')
}
isSampleCol_mse <- function(mse, col){
    col %in% sampleNames_mse(mse)
}

# help perform a complex rowData filter by passing rowData to a callback
# such a query is applied to the metadata of each data row and passed as 'rows' above
filterRowData_mse <- function(mse, filterFn, as.dt=TRUE){
    rowData <- SummarizedExperiment::rowData(mse$rse())
    filterFn( if(as.dt) as.data.table(rowData) else rowData )
}

# help perform a complex sample value filter by passing a data.table to a callback
# such a query is applied to the various assays of a given sample and passed as 'rows' above
filterSampleData_mse <- function(mse, col, assays=NULL, filterFn){
    filterFn( getValuesCols_mse(mse, col, assays) )
}
