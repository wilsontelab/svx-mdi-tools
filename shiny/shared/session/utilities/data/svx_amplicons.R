# parse amplicons for rapid navigation to the relevant genome location(s) (one or two per amplicon)
svx_getTrackAmpliconTargets <- function(track, parseFn = NULL, isMultiSample = TRUE, amplicons = NULL){
    if(is.null(amplicons)){
        amplicons <- track$settings$items()
        if(isMultiSample) amplicons <- getSourcesFromTrackSamples(amplicons)
    }
    nAmplicons <- length(amplicons)
    req(nAmplicons > 0) 
    amplicons <- do.call(rbind, lapply(amplicons, as.data.table))
    if(!is.null(parseFn)) amplicons <- parseFn(amplicons) # must yield at least columns ampliconKey, sampleName, chrom1/2, strand1/2, pos1/2
    navTargets <- do.call(rbind, lapply(1:nAmplicons, function(i){
        x <- amplicons[i]
        if(x$chrom1 == x$chrom2 && abs(x$pos2 - x$pos1) <= 10000){
            x[, .(ampliconTargetKey = ampliconKey, chrom = chrom1, start = min(pos1, pos2), end = max(pos1, pos2))]
        } else {
            data.table(
                ampliconTargetKey = paste(x$ampliconKey, 1:2), 
                chrom = c(x$chrom1, x$chrom2),
                start = c(x$pos1, x$pos2),
                end   = NA_integer_
            )
        }
    }))
    unique(navTargets)
}

# parse amplicon junctions into SVX compatible data.table
svx_parseTrackAmplicons <- function(x){
    setnames(x, c(
        "ampliconId","ampliconType","count",
        "chrom1","side1","pos1","ref1","primer1",
        "chrom2","side2","pos2","ref2","primer2"
    ))
    x[, ":="(
        strand1 = ifelse(side1 == "R", "+", "-"),
        strand2 = ifelse(side2 == "L", "+", "-")
    )]
    x
}
