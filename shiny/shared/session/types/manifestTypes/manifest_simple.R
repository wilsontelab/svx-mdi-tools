#----------------------------------------------------------------------
# manifestType 'simple' assumes
#   - file with suffix 'csv', 'txt', or 'gz'
#   - columns Project, Sample_ID, and Description already present
#   - one row per sample
#   - no further column data processing or aggregation is necessary
#----------------------------------------------------------------------

manifestTypes$simple <- list(
    patterns = c(
        'csv',
        'txt',
        'gz'
    ),
    load = fread,
    parse = function(manifest){
        list(
            manifest = manifest,
            unique   = manifest
        )
    }
)
