#----------------------------------------------------------------------
# helper function to create a minimal Bioconductor RangedSummarizedExperiment
# from an input data.table; only actions required to execute that
# conversion are performed here
#----------------------------------------------------------------------
# the input data.table must have meaningful column names, typically
# acquired from the header of a flat file using fread
#----------------------------------------------------------------------

dataTableToRSE <- function(dt, manifest,
                           BED,       # BED=TRUE means start is 0-referenced
                           stranded){ # stranded=TRUE means features are strand specific
    reportProgress('dataTableToRSE')

#----------------------------------------------------------------------
# interpret the sample columns and naming format
#----------------------------------------------------------------------
reportProgress('parseColumnNames')
columnNames <- names(dt)
countMatchingColumns <- function(x) sum(x %in% columnNames)
candidateSampleColumnNames <- list( # the different ways an incoming file might be headed
    rawSampleId       = manifest$Sample_ID,
    prefixed_SampleId = paste('sample', manifest$Sample_ID, sep = '_'),
    Prefixed_SampleId = paste('Sample', manifest$Sample_ID, sep = '_'),
    PREFIXED_SampleId = paste('SAMPLE', manifest$Sample_ID, sep = '_'),
    description       = manifest$Description,
    indices           = seq_along(manifest$Sample_ID)
)
matchCounts <- sapply(candidateSampleColumnNames, countMatchingColumns)
bestMatchI <- which.max(matchCounts)
sampleNameFormat <- names(candidateSampleColumnNames)[bestMatchI] 
allowedSampleColumnNames <- candidateSampleColumnNames[[bestMatchI]] 
isSampleColumn <- columnNames %in% allowedSampleColumnNames
sampleColumnIs <- (max(which(!isSampleColumn)) + 1):max(which(isSampleColumn)) # the columns with sample values

#----------------------------------------------------------------------
# extract the sample columns and prepare colData for RangedSummarizedExperiment
#----------------------------------------------------------------------
reportProgress('extract assay raw')
assays <- list(raw = dt[, .SD, .SDcols = columnNames[sampleColumnIs]]) # raw = the original table values
sampleColumnNames <- columnNames[sampleColumnIs]
colData <- data.frame(
    Sample_ID = sapply(sampleColumnNames, function(x){
        switch(sampleNameFormat,
            rawSampleId       = x,
            prefixed_SampleId = sub('sample_', '', x),
            Prefixed_SampleId = sub('Sample_', '', x),
            PREFIXED_SampleId = sub('SAMPLE_', '', x),
            description       = manifest[manifest$Description == x, 'Sample_ID'],
            indices           = manifest$Sample_ID[as.integer(x)]  
        )
    }),    
    stringsAsFactors = FALSE
)
row.names(colData) <- colData$Sample_ID
names(assays$raw)  <- colData$Sample_ID

#----------------------------------------------------------------------
# extract the bin columns and prepare rowRanges for RangedSummarizedExperiment
# output is always 1-referenced per R convention and stranded upon request
#----------------------------------------------------------------------
reportProgress('extract row data')
rowData <- dt[, .SD, .SDcols = columnNames[-sampleColumnIs]]
dt <- renameChromStartEnd(rowData)
reservedNames <- c( # from Bioconductor; if present in data.table, rename to seqnames_ etc. in rowData
    "seqnames", "ranges",
    "seqlevels", "seqlengths", "isCircular", 
    "width", "element"
)
if(!stranded) reservedNames <- c(reservedNames, 'strand')
names(rowData) <- sapply(names(rowData), function(x) ifelse(x %in% reservedNames, paste0(x, '_'), x))
rowRanges <- GenomicRanges::makeGRangesFromDataFrame(
    rowData,
    keep.extra.columns = TRUE,
    ignore.strand = !stranded,
    starts.in.df.are.0based = BED
)
if('name' %in% names(rowData)) names(rowRanges) <- rowData[, name]

#----------------------------------------------------------------------
# assemble and return our RangedSummarizedExperiment object
#----------------------------------------------------------------------
SummarizedExperiment::SummarizedExperiment(
    assays    = assays,
    rowRanges = rowRanges,
    colData   = colData,
    metadata = list(
        inputColumnFormat = sampleNameFormat
    )
)

}
