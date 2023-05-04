fastqIndex <- list()
getReadFromSequenceIndex <- function(sourceId, qName_){
    sourceId <- sourceId()
    req(sourceId)

    # parse the input files; abort if not found (e.g., if not on analysis server)
    extract <- read_yaml(file.path(getPackageDir(sourceId), "package.yml"))$task$extract
    taskDir <- file.path(extract$output[["output-dir"]], extract$output[["data-name"]])
    sequencesFile <- paste(extract$output[["data-name"]], extract$genome$genome, "extract.sequences.txt", sep = ".")
    sequencesFile <- file.path(taskDir, sequencesFile)
    indexFile <- paste(sequencesFile, "gz.index", sep = ".") 
    req(file.exists(sequencesFile), file.exists(indexFile))

    # parse the index to the read
    if(is.null(fastqIndex[[indexFile]])) {
        x <- fread(indexFile)
        setnames(x, c("qName","offset","length"))
        fastqIndex[[indexFile]] <<- x
    }
    index <- fastqIndex[[indexFile]][qName == qName_]
    req(nrow(index) > 0)

    # pull the read data from the file
    con <- file(sequencesFile, "rb")    
    seek(con, index$offset)
    read <- as.list(strsplit(readLines(con, 1), "\t")[[1]])
    names(read) <- c("type","qName","seq","qual")
    close(con)

    # extract per-base QUAL
    read$QUAL <- as.integer(sapply(strsplit(read$qual, "")[[1]], charToRaw)) - 33L
    read
}
