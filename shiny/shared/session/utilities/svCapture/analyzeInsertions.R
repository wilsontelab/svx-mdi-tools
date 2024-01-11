#----------------------------------------------------------------------
# analyze the insertions in a set of pre-filtered SVs
#----------------------------------------------------------------------
svx_analyzeInsertion <- function(
    microhomologyLength, minRefWidth, padding,
    JXN_BASES, GEN_REF_1, GEN_REF_2, SIDE_1, SIDE_2
){
    # flip one reference strand for inversions to match jxnSeq assembly
    if(SIDE_1 == SIDE_2){ 
        if(SIDE_1 == "L") GEN_REF_2 <- rc(GEN_REF_2)
                     else GEN_REF_1 <- rc(GEN_REF_1)
    }

    # parse the insertion search sequence as (microhomology)(insertion)(microhomology)
    thisRefWidth <- nchar(GEN_REF_1)
    thisPadding <- (thisRefWidth - 1) / 2
    microhomology1 <- substr(GEN_REF_1, 
                             thisPadding - microhomologyLength + 2,     
                             thisPadding + 1)
    microhomology2 <- substr(GEN_REF_2, 
                             thisPadding + 1, 
                             thisPadding + microhomologyLength)
    searchSeq <- paste0(microhomology1, JXN_BASES, microhomology2)
    rcSearchSeq <- rc(searchSeq)

    # trim the genome references to the mininum available sample padding
    if(padding < thisPadding){
        start <- thisPadding - padding + 1
        stop  <- start + minRefWidth - 1
        GEN_REF_1 <- substr(GEN_REF_1, start, stop)
        GEN_REF_2 <- substr(GEN_REF_2, start, stop)
    }

    # return all matches to the search sequence
    c(
        searchSeq = searchSeq,
        match1   = paste(gregexpr(  searchSeq, GEN_REF_1)[[1]], collapse = ","),
        match1rc = paste(gregexpr(rcSearchSeq, GEN_REF_1)[[1]], collapse = ","),        
        match2   = paste(gregexpr(  searchSeq, GEN_REF_2)[[1]], collapse = ","),
        match2rc = paste(gregexpr(rcSearchSeq, GEN_REF_2)[[1]], collapse = ",") 
    )
}
svx_analyzeInsertions <- function(svs, properties){
    req(svs)
    svs <- svs[
        -MICROHOM_LEN >= properties$Min_Insertion_Size$value & 
        -MICROHOM_LEN <= properties$Max_Insertion_Size$value
    ]    
    minRefWidth <- svs[, min(nchar(GEN_REF_1))]    
    padding <- (minRefWidth - 1) / 2
    svs <- cbind(svs, t(svs[, mapply(
        svx_analyzeInsertion, 
        properties$Flanking_Microhomology$value, minRefWidth, padding, 
        JXN_BASES, GEN_REF_1, GEN_REF_2, SIDE_1, SIDE_2
    )]))
    svs[, ':='(
        nCombinations = 4 ** (-MICROHOM_LEN + properties$Flanking_Microhomology$value * 2), 
        nMatches = mapply(function(match1, match1rc, match2, match2rc){
            sum(strsplit(match1,   ",")[[1]] != "-1") + # total number of times search sequence was found
            sum(strsplit(match1rc, ",")[[1]] != "-1") + 
            sum(strsplit(match2,   ",")[[1]] != "-1") + 
            sum(strsplit(match2rc, ",")[[1]] != "-1")
        }, match1, match1rc, match2, match2rc)
    )]
    svs[, found := nMatches > 0] # whether or not each SV's search sequence was found
    list(
        refWidth = minRefWidth,
        padding = padding,
        svs = svs
    )
}
