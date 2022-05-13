
# use CIGAR string and SEQ to return an alignment's query bases relative to the reference genome
# output is the same length as reference (not the aligned read)
cigarToRefAln <- function(cigar, seq){
    
    # unpack the CIGAR string
    lengths    <- as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar))))
    operations <-            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    bases <- strsplit(seq, '')[[1]]
    i <- 1 # index of SEQ, i.e. the query that was aligned to reference    
    lengthAlignedOut <- 0 # the length of just the aligned portion of the output (i.e. without any clips)
    nDeleted <- 0 # number of bases deleted to the left of i
    insertionPos <- c() # indices in 'out' to the left of insertions
    insertionSeq <- c()

    # return simple alignments as is
    if(length(lengths) == 1) return(list(
        seq = bases,
        map = bases,
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
            insertionSeq <<- c(insertionSeq, paste(bases[i:(i+length-1)], collapse=""))
            i <<- i + length
            c()
            
        # M or m operations present in reference and query, output differs by case
        # S operations are treated like M/m bases, outside of POS
        } else {
            subSeq <- bases[i:(i+length-1)]
            if(operation == 'S') subSeq <- tolower(subSeq)
                else lengthAlignedOut <<- lengthAlignedOut + length
            i <<- i + length            
            subSeq 
        }
    }, lengths, operations))
    
    # return the result
    map <- out
    map[insertionPos] <- '+' # mark the position to the left of a novel insertion and append its bases there
    for(i in seq_along(insertionPos)) out[insertionPos[i]] <- paste0(out[insertionPos[i]], insertionSeq[i])
    list(
        seq = out, # map of the bases, similar to IGV in concept
        map = map,
        leftClip = if(operations[1] == 'S') lengths[1] else 0, # help for alignment relative to POS
        lengthOut = length(out), # length of the query on reference (not necessarily length(query))
        lengthAlignedOut = lengthAlignedOut
    )
}

# I generated these examples by manually editing the reference sequence
# and then using BWA to align to hg38
examples <- list(
    c("reference","100M",            "AGGGCACGTAGACCCGGATCCCCAGCTCCCGGACCCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG"),
    c("SNP","100M",                  "AGGGCACGTAGACCCGGATCCCCAGATCCCGGACCCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG"),
    c("simple_deletion","25M1D74M",  "AGGGCACGTAGACCCGGATCCCCAGTCCCGGACCCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG"),
    c("simple_insertion","26M1I74M", "AGGGCACGTAGACCCGGATCCCCAGCATCCCGGACCCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG"),
    c("complex_deletion","24M3D73M", "AGGGCACGTAGACCCGGATCCCCACCCGGACCCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG"),
    c("complex_insertion","22M3I78M","AGGGCACGTAGACCCGGATCCCCAGCAGCTCCCGGACCCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG"),
    c("deletion_in_run","33M1D66M",  "AGGGCACGTAGACCCGGATCCCCAGCTCCCGGACCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG"),
    c("insertion_in_run","33M1I67M", "AGGGCACGTAGACCCGGATCCCCAGCTCCCGGACCCCCTTCCGACTGCTCCCTTCGCCACTCGCGCGGCTTCCACTCCTGGACTCGGCCTCCACCCCTGGG")
)

# the reference base vector
refV <- strsplit(examples[[1]][3], '')[[1]]

# show how cigarToRefAln parses different kinds of variants
for(example in examples){
    x <- cigarToRefAln(example[2],example[3])
    message()
    message(example[1])
    message(example[2])
    print(data.frame(ref=refV, seq=x$seq, map=x$map, match=ifelse(refV==x$map, "", "***")))
}

