# support functions for parsing cigar strings and positions

# throughout, read1 has POS at the start of query, read2 has POS at the end

# the number of reference bases covered by a CIGAR string
getNRefBases <- function(cigar){
    lengths    <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    operations <-            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    i <- operations %in% c("M","D")    
    sum(lengths[i])
}

# when POS is at the end of an alignment, use CIGAR to return the start position
getAlnStartPos <- function(pos, cigar){
    pos - getNRefBases(cigar) + 1
    # paste0(c(rbind(lengths[i], operations[i])), collapse = "")
}

# return the leftmost position of an alignment from its metadata
getRefPos <- function(readN, strand, pos, cigar){
    if(readN == 1){
        if(strand == "+") pos else getAlnStartPos(pos, cigar)
    } else {
        if(strand == "-") pos else getAlnStartPos(pos, cigar)
    }
}
