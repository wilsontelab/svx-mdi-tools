
#----------------------------------------------------------------------
# ../constants.R defines universal values
#----------------------------------------------------------------------

# chromosomes
CHROMS <- c(paste('chr', 1:22, sep=""), 'chrX', 'chrY')

# filter constants
allPairTypes    <- c('TT','TA','AA','tt','ta','aa','t-','a-','--')
targetPairTypes <- c('TT','TA','AA','tt','ta','aa','t-','a-') # for Targets server
defaultTargetPairTypes <- c('TT','TA','AA','tt','ta','aa') # for Targets server
targetPairTypes_peaks <- c('TT','t-') # for Peaks server
defaultTargetPairTypes_peaks <- targetPairTypes_peaks # for Peaks server

# color option values and other fixed plot parameters
#uHomLim <- c(-20,50)
uHomLim <- c(-30,30)
uHomLab <- "Microhomology(bp)"
plotColors <- list(
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
minSVSize_peak <- 25000

