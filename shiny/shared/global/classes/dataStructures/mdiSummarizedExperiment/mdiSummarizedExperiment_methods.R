#----------------------------------------------------------------------
# mdiSummarizedExperiment class generic methods, called as method(obj)
# these methods are updated on each session and never get stale,
# therefore, they should not make calls to the object's static methods
# instead, generics should call utility_mse functions
#----------------------------------------------------------------------

# retrieve an XY table for plotting in a scatter plot, etc.
# optionally append other named row information in addition to x and y columns
getXY.mdiSummarizedExperiment <- function(mse, xcol, ycol, assay=NULL, # basic info for required xy columns
                                          mcols=FALSE, rowData=FALSE, assays=NULL, # additional data to tack on
                                          rows=TRUE){ # row filter
    dt <- data.table(
        x = getRowDataCol_mse(mse, xcol), 
        y = getValueCol_mse(mse, ycol, assay)
    )
    if(rowData) dt <- cbind(SummarizedExperiment::rowData(mse$rse()), dt)
    if(mcols) dt <- cbind(SummarizedExperiment::mcols(mse$rse()), dt)
    if(!is.null(assays)) dt <- cbind(dt, getValuesCols_mse(mse, ycol, assays = assays))
    dt[rows] 
}

# retrieve a single column of values as a vector; work hard to find the column the user wants
getCol.mdiSummarizedExperiment <- function(mse, col = NULL, assay = NULL, rows = TRUE){
    if(isRowDataCol_mse(mse, col)) getRowDataCol_mse(mse, col, rows = rows)
    else getValueCol_mse(mse, col, assay, rows = rows)
}

# matrix dimensions
nFeatures.mdiSummarizedExperiment <- function(mse){ # nolint
    nrow(mse$getRowData())
}
nSamples.mdiSummarizedExperiment <- function(mse){ # nolint
    nrow(mse$getColData())
}
