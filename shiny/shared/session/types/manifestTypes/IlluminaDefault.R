#----------------------------------------------------------------------
# standard Illumina demux file, sometimes with multiple rows per sample
#----------------------------------------------------------------------

# the columns we will add to the manifestKeyColumns in returned tables
IlluminaSampleValueTemplate <- data.frame(
    PF_Clusters  = integer(),     # reported as Yield
    Yield_Mbases = integer(),
    Mean_QualityScore = numeric() # reported as Quality
)
IlluminaSampleValueColumns <- colnames(IlluminaSampleValueTemplate)
IlluminaSampleValueMeanColumns <- c('Mean_QualityScore')

# aggregate sample statistics over all parts of an Illumina run
parseIlluminaManifest <- function(manifest){

    # fix any naming inconsistences in the manifest file header
    names(manifest) <- fixColumnNames(names(manifest), 'Sample_ID',    c('SampleID', 'SampleId'))
    names(manifest) <- fixColumnNames(names(manifest), 'Yield_Mbases', c('Yield_.Mbases.'))

    # fix any column data type inconsistencies (e.g., due to commas in numbers)
    manifest <- fixColumnDataTypes(manifest, IlluminaSampleValueTemplate)

    # correct NA values for any bad/failed libraries in the manifest
    # keeps the bad samples, with 0, not NA, for yields and qualities
    for(col in IlluminaSampleValueColumns) manifest[is.na(manifest[[col]]), col] <- 0

    # calculate the sum of the columns
    unique <- aggregate( 
        . ~ Project + Sample_ID + Description,
        data = manifest[, c(CONSTANTS$manifestKeyColumns,
                            IlluminaSampleValueColumns)],
        sum
    )

    # then finish to mean for required columns
    for(col in IlluminaSampleValueMeanColumns){
        unique[[col]] <- sapply(seq_len(nrow(unique)), function(i){
            nValues <- sum(manifest$Sample_ID == unique$Sample_ID[i])
            unique[[col]][i] / nValues
        })        
    }

    # reorder samples since aggregate changes the sort
    uniqueSampleIds <- uniqueRows(manifest, CONSTANTS$manifestKeyColumns)$Sample_ID
    unique <- unique[order(match(unique$Sample_ID, uniqueSampleIds)), ]
    rownames(unique) <- seq_len(nrow(unique))

    # set the names Yield and Quality columns    
    unique$Yield   <- unique$PF_Clusters
    unique$Quality <- round(unique$Mean_QualityScore, 1)

    # return the results
    list(
        manifest = manifest, # return our input in case we modified the column names or data types
        unique = unique
    )
}

# define the default Illumina manifest type
manifestTypes$IlluminaDefault <- list(
    
    # only accept files ending with these suffixes
    patterns = c("_DemuxStats.csv$", "_DemuxStats.txt$"),

    # function to load the manifest
    # no comments allowed, '#' is allowed character in sample names
    load = function(file){
        sep <- if(endsWith(file, 'csv')) "," else "\t"
        read.table(file, header = TRUE, sep = sep, stringsAsFactors = FALSE, comment.char = "")
    },
    
    # reduce the manifest to unique sample rows
    parse = parseIlluminaManifest

)
