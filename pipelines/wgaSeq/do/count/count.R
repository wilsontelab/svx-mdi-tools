# make a table with summary counts, alignment and duplication rates
# for each cell sample in the overall data set

# load resources
library(jsonlite)
env <- as.list(Sys.getenv())

# load counts data
d <- read.table(env$COUNTS_FILE, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
names(d) <- c('sampleName', 'filtered', 'grouped')

# collect read pair counts from fastp
d$readPairs <- sapply(d$sampleName, function(sampleName){
    json <- paste(env$LOG_FILE_PREFIX, sampleName, 'fastp', 'json', sep = ".")
    json <- read_json(json)
    json$summary$before_filtering$total_reads / 2
})

# add a few read pair rates
d$alignRate <- d$filtered / d$readPairs   # fraction of read pairs that successfully aligned
d$dupRate   <- 1 - d$grouped / d$filtered # fraction of accepted read pairs that were discarded as replicates

# load and add source cell labels, if provided with manifest-file
fillNullValues <- function(){
    d$libraryName <<- d$sampleName
    d$cellName <<- d$sampleName
    d$sampleN <<- seq_along(nrow(d)) 
}
if(env$MANIFEST_FILE == "NA"){
    fillNullValues()
} else { # parse sample numbers and cell names based on AGC patterns
    tryCatch({
        d$libraryName <- gsub('Sample_', '', d$sampleName)
        d$sampleN <- as.integer(sapply(strsplit(d$libraryName, '-'), function(v) v[length(v)]))
        manifest <- read.csv(env$MANIFEST_FILE, header = TRUE, stringsAsFactors = FALSE)
        manifest <- manifest[manifest$Lane == 1, 3:4]
        names(manifest) <- c('libraryName', 'cellName')
        d <- merge(d, manifest, by = 'libraryName', all.x = TRUE, all.y = FALSE)        
    }, error = function(e) {
        message("error attempting to parse AGC data from manifest file")
        print(e)
        fillNullValues()
    })
}

# write out results
d <- d[order(d$sampleN),
       c('sampleN', 'sampleName', 'libraryName', 'cellName',
         'readPairs', 'filtered', 'grouped',
         'alignRate', 'dupRate')]
write.table(
    d, 
    env$RATES_FILE, 
    quote = FALSE, 
    sep = "\t", 
    row.names = FALSE, 
    col.names = TRUE
)

# create a simple manifest file for the app
project <- if(env$MANIFEST_FILE == "NA"){
    basename(env$INPUT_DIR)
} else { 
    tryCatch({
        manifest <- read.csv(env$MANIFEST_FILE, header = TRUE, stringsAsFactors = FALSE)
        manifest[manifest$Lane == 1, 'Project']
    }, error = function(e) {
        basename(env$INPUT_DIR)
    })
}
write.table(
    data.frame(
        Project = project,
        Sample_ID = d$libraryName,
        Description = d$cellName,
        stringsAsFactors = FALSE
    ),
    env$SIMPLE_MANIFEST_FILE, 
    quote = FALSE, 
    sep = ",", 
    row.names = FALSE, 
    col.names = TRUE
)
