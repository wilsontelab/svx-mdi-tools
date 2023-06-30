# support functions for parsing cigar strings and positions

# throughout, read1 has POS at the start of query, read2 has POS at the end

# the number of reference bases covered by a CIGAR string
getNRefBases <- function(cigar){
    lengths    <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    operations <-            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    i <- operations %in% c("M","D")    
    sum(lengths[i])
}

# return the first node position of an outer alignment
# working from its second node position
getAlnStartPos <- function(node2Pos, strand, cigar){
    if(strand > 0) node2Pos - getNRefBases(cigar) + 1
    else           node2Pos + getNRefBases(cigar) - 1
}
# return the second node position of an outer alignment
# working from its first node position
getAlnEndPos <- function(node1Pos, strand, cigar){
    if(strand > 0) node1Pos + getNRefBases(cigar) - 1
    else           node1Pos - getNRefBases(cigar) + 1
}

# return the leftmost position of an alignment from its metadata
getRefPos <- function(strand, pos, cigar){
    if(strand == 1) pos else getAlnStartPos(pos, cigar)
}
