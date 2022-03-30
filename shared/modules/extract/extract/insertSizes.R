# get environment variables
env <- as.list(Sys.getenv(c("EXTRACT_PREFIX", "PLOT_PREFIX", "READ_LEN", "MAX_TLEN"), names = TRUE))
env$READ_LEN <- as.numeric(env$READ_LEN)
env$MAX_TLEN <- as.numeric(env$MAX_TLEN)

# get input data
counts <- read.table(file = "stdin", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
countNames  <- c("untargeted", "TT", "TA", "AA", "xx")
legendNames <- c("untargeted", "TT", "TA", "AA", "--")
names(counts) <- c("insertSize", countNames)

# convert counts to frequencies
freqs <- counts
sums  <- list()
for (name in countNames) freqs[[name]] <- {
    sums[[name]] <- sum(freqs[[name]])
    if(sums[[name]] > 0) freqs[[name]] / sums[[name]] else 0
}

# set the plot scale and order
maskReadLens <- env$READ_LEN + -2:2
scaleFreq <- freqs[!(freqs$insertSize %in% maskReadLens), countNames]
maxFreq <- max(scaleFreq)
orderedNames <- names(sums)[order(unlist(sums))]

# initialize plot histogram with thresholds
pngFile <- paste(env$PLOT_PREFIX, "insertSizes", "png", sep = ".")
png(pngFile,
    width = 4, height = 4, units = "in", pointsize = 10,
    res = 600, type = c("cairo"))
plot(
    0, 0, typ = "n", 
    xlab = "Insert Size (bp)", 
    ylab = "Frequency",
    xlim = c(0, env$MAX_TLEN) * 1.1, 
    ylim = c(0, maxFreq * 1.05)
)
abline(v = c(env$READ_LEN, env$READ_LEN * 2, env$MAX_TLEN), col = "grey")

# plot all relevant data traces
minCountToPlot <- 1000
cols <- c("blue", "green4", "blue", "red3", "gray50")
plottedNames <- character()
plottedCols  <- integer()
for(name in orderedNames){
    i <- which(countNames == name)
    col  <- cols[i]
    if(sums[[name]] >= minCountToPlot){
        lines(freqs$insertSize, freqs[[name]], col = col)
        plottedNames <- c(plottedNames, legendNames[i])
        plottedCols  <- c(plottedCols,  col)
    }
}

# add a legend
legend("topright", plottedNames, col = plottedCols, lty = 1, lwd = 2)
