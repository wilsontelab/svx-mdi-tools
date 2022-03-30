#--------------------------------------------------------------
# define global constants
#--------------------------------------------------------------
# you must EXTEND not REPLACE the framework's global CONSTANTS object!
# e.g., CONSTANTS$name <- value
#--------------------------------------------------------------

# chromosomes
# CONSTANTS$CHROMS <- c(paste('chr', 1:22, sep=""), 'chrX', 'chrY')

# filter constants
CONSTANTS$allPairTypes    <- c('TT','TA','AA','tt','ta','aa','t-','a-','--')
CONSTANTS$targetPairTypes <- c('TT','TA','AA','tt','ta','aa','t-','a-') # for Targets
CONSTANTS$defaultTargetPairTypes <- c('TT','TA','AA','tt','ta','aa')    # for Targets
CONSTANTS$targetPairTypes_peaks  <- c('TT','t-') # for Peaks
CONSTANTS$defaultTargetPairTypes_peaks <- targetPairTypes_peaks # for Peaks

# color option values and other fixed plot parameters
CONSTANTS$uHomLim <- c(-30,30)
CONSTANTS$uHomLab <- "Microhomology(bp)"
CONSTANTS$plotColors <- list(
    duplex = list(
        duplex = "green4",
        single = "red3"
    ),
    svType = list(
        Dup = "red3",
        Del = "green4",
        Inv = "blue",
        Trans = "black",
        Undet = "purple"
    ),
    sample = c(
        "red3",
        "green4",
        "blue",
        "purple",
        "black"
    )
)

# SV parameters
CONSTANTS$minSVSize_peak <- 25000
