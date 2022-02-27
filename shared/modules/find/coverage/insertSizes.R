# get environment variables
env <- as.list(Sys.getenv(c('EXTRACT_PREFIX', 'PLOT_PREFIX', 'READ_LEN', 'MAX_TLEN'), names=TRUE))
env$READ_LEN <- as.numeric(env$READ_LEN)
env$MAX_TLEN <- as.numeric(env$MAX_TLEN)

# get input data
d <- read.table(paste(env$EXTRACT_PREFIX, "insertSizes", "txt", sep="."),
                header=TRUE, sep="\t", stringsAsFactors=FALSE)

# initialize plot histogram with thresholds
pngFile <- paste(env$PLOT_PREFIX, "insertSizes", "png", sep = ".")
png(pngFile,
    width = 4, height = 4, units = "in", pointsize = 10,
    res = 600, type = c("cairo"))
maskReadLens <- env$READ_LEN + -2:2
maxFreq <- max( d$frequency[!(d$insertSize %in% maskReadLens)] )
plot(
    0, 0, typ = "n", 
    xlab = "Insert Size (bp)", 
    ylab = "Frequency",
    xlim = c(0, env$MAX_TLEN) * 1.1, 
    ylim = c(0, maxFreq * 1.05)
)
abline(v = c(env$READ_LEN, env$READ_LEN*2, env$MAX_TLEN), col = "grey")
lines(d$insertSize, d$frequency, col = "blue")
