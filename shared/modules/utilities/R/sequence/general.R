
# utilities for sequence and alignment processing


# standard SAM columns and flags
SAM_columns <- list(
    'QNAME' = 'character',  
    'FLAG'  = 'integer',
    'RNAME' = 'character',  
    'POS'   = 'integer',     
    'MAPQ'  = 'integer',
    'CIGAR' = 'character',
    'RNEXT' = 'character',
    'PNEXT' = 'integer',
    'TLEN'  = 'integer',
    'SEQ'   = 'character',
    'QUAL'  = 'character'
)
SAM_flags <- list(
    'PAIRED' = 1,  
    'PROPER_PAIR'  = 2,
    'UNMAPPED' = 4,  
    'MATE_UNMAPPED'   = 8,     
    'REVERSE'  = 16,
    'MATE_REVERSE' = 32,
    'FIRST_IN_PAIR' = 64,
    'SECOND_IN_PAIR' = 128,
    'SECONDARY'  = 256,
    'FAILED_QC'   = 512,
    'DUPLICATE'  = 1024,
    'SUPPLEMENTAL' = 2048
)

# use CIGAR string and SEQ to return an alignment's query bases relative to the reference genome
# output is the same length as reference (not the aligned read)
cigarToRefAln <- function(cigar, seq){
    
    # unpack the CIGAR string
    lengths    <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    operations <-            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    bases <- strsplit(seq, '')[[1]]
    i <- 1 # index of SEQ, i.e. the query that was aligned to reference    
    lengthAlignedOut <- 0 # the length of just the aligned portion of the output (i.e., without any clips)
    nDeleted <- 0 # number of bases deleted to the left of i
    insertionPos <- c() # indices in 'out' to the left of insertions

    # return simple alignments as is
    if(length(lengths) == 1) return(list(
        seq = bases,
        leftClip = 0,
        lengthOut = lengths,
        lengthAlignedOut = lengths
    ))

    # run through all operations
    out <- unlist(mapply(function(length, operation){
        
        # deleted segments are present in reference, not in query
        if(operation == 'D'){
            nDeleted <<- nDeleted + length
            lengthAlignedOut <<- lengthAlignedOut + length
            rep('-', length)
            
        # inserted bases are present in query, but not reference
        } else if(operation == 'I'){
            insertionPos <<- c(insertionPos, i - 1 + nDeleted)
            i <<- i + length
            c()
            
        # M or m operations present in reference and query, output differs by case
        # S operations are treated like M/m bases, outside of POS
        } else {
            subSeq <- bases[i:(i + length - 1)]
            if(operation == 'S') subSeq <- tolower(subSeq)
                else lengthAlignedOut <<- lengthAlignedOut + length
            i <<- i + length            
            subSeq 
        }
    }, lengths, operations))
    
    # return the result
    out[insertionPos] <- '+' # mark the position to the left of a novel insertion
    list(
        seq = out, # map of the bases, similar to IGV in concept
        leftClip = if(operations[1] == 'S') lengths[1] else 0, # help for alignment relative to POS
        lengthOut = length(out), # length of the query on reference (not necessarily length(query))
        lengthAlignedOut = lengthAlignedOut
    )
}

# get rightmost mapped read position in reference genome from POS and CIGAR
getEnd <- Vectorize(function(start, cigar) {
    lengths    <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    operations <-            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    allowed <- operations %notin% c('S', 'I')
    start - 1 + sum(lengths[allowed])
})

# get soft clips
getLeftClip <- function(cigar){
    operations <- unlist(regmatches(cigar, gregexpr('\\D',  cigar)))     
    if(operations[1] != 'S') return(0)
    lengths <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    lengths[1]
}
getRightClip <- function(cigar){
    operations <- unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    maxI <- length(operations)
    if(operations[maxI] != 'S') return(0)
    lengths <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    lengths[maxI]
}

# reverse complement ACGTN and associated CIGAR strings
rc_bases <- list(A = "T", C = "G", G = "C", T = "A", N = "N")
rc <- Vectorize(function(SEQ){
    paste(rev(unlist(rc_bases[ strsplit(SEQ, '')[[1]] ])), collapse = "")
})
rc_cigar <- Vectorize(function(cigar){
    operations <- unlist(regmatches(cigar, gregexpr('\\d+\\D',  cigar)))
    paste(rev(operations), collapse = "")
})
