# determine whether a sequence position is in a bad genome region
# requires data.table
BAD_REGIONS <- list()
loadBadRegions <- function(){
    if(is.null(env$BAD_REGIONS_FILE)) return(NULL)
    x <- fread(env$BAD_REGIONS_FILE, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
    x <- x[, 1:3]
    names(x) <- c('chrom', 'start', 'end')
    for(chrom_ in x$chrom){
        BAD_REGIONS[[chrom_]] <<- x[chrom == chrom_, list(start, end)]
        setkey(BAD_REGIONS[[chrom_]], start, end)
    }
}
isBadRegion <- Vectorize(function(chrom, pos){
    if(is.null(BAD_REGIONS[[chrom]])) return(FALSE)
    BAD_REGIONS[[chrom]][start <= pos & end >= pos, .N] > 0
})
    #any(pos >= BAD_REGIONS[[chrom]]$start &
    #    pos <= BAD_REGIONS[[chrom]]$end)
