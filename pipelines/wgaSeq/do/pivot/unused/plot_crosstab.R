
# use parallel processing
library(parallel)

# initialize script
env <- as.list(Sys.getenv())
source(paste(env$PIPELINE_DIR, "common", "utilities.R", sep="/"))

# load cell labels
manifest <- read.csv(env$MANIFEST_FILE, header=TRUE, stringsAsFactors=FALSE)
manifest <- manifest[manifest$Lane==1,3:4]
names(manifest) <- c('libraryName','cellName')
manifest$sampleName <- paste('Sample', manifest$libraryName, sep="_")
manifest$sampleN <- as.integer(sapply(strsplit(manifest$libraryName, '-'), function(v) v[length(v)]))
manifest <- manifest[order(manifest$sampleN),]
rownames(manifest) <- sapply(manifest$sampleName, function(sampleName){
    gsub('-','.',sampleName)
})
manifest$adjustedSampleName <- rownames(manifest)
nInputCells <- nrow(manifest)

# load crosstab csv
#chrom,start,end,binN,all_cells,strand,Sample_2251-LK-1,...
csv <- paste('crosstab', env$BIN_SIZE, 'csv', sep=".")
csv <- paste(env$OUTPUT_DIR, csv, sep="/")
ct <- read.csv(csv, header=TRUE, stringsAsFactors=FALSE)
env$BIN_SIZE <- as.integer(env$BIN_SIZE)

# filter and re-order to canonical chromosomes
chroms <- paste0('chr', c(as.character(1:22), 'X')) #, 'Y'
chroms <- chroms[chroms %in% unique(ct$chrom)]
nChroms <- length(chroms)
ct <- do.call('rbind', lapply(chroms, function(chrom) ct[ct$chrom==chrom,]))

# set cell indices
cellJs <- 7:(ncol(ct)-3)
nOutputCells <- length(cellJs)
cellKs <- cellJs - 6
adjustedSampleNames <- colnames(ct)[cellJs]
sampleNs  <- sapply(adjustedSampleNames, function(asn) manifest$sampleN[manifest$adjustedSampleName == asn])
cellNames <- sapply(adjustedSampleNames, function(asn) manifest$cellName[manifest$adjustedSampleName == asn])
cols <- rainbow(nInputCells)[sampleNs]

# enforce minimum mappability and effective bin size: TODO, expose as options
# identify and exlude a set of further unreliable (always low) bins
minMappability <- 0.8
minEffectiveSize <- env$BIN_SIZE / 2
ct$binSize <- env$BIN_SIZE - ct$excluded_bases
autosomes <- ct$chrom %notin% c('chrX', 'chrY')
toFit <- data.frame(
    x  = ct[autosomes,'gc_fraction'],
    w  = ct[autosomes,'binSize'],
    y  = rowSums(ct[autosomes,cellJs])
)
toPredict <- data.frame(
    x  = ct[,'gc_fraction'],
    y  = rowSums(ct[,cellJs])
)
fit <- loess(y ~ x, toFit, weights=toFit$w) # more prone to outlier influence
exp0 <- predict(fit, toPredict)
rpa  <- exp0 / 2     
cn   <- toPredict$y / rpa
minCN <- 0.5
ct <- ct[ct$mappability >= minMappability &
         ct$binSize >= minEffectiveSize &
         cn >= minCN,]

# done filtering bins, set additional needed parameters
nBins <- nrow(ct)
nChroms <- length(unique(ct$chrom))
autosomes <- ct$chrom %notin% c('chrX', 'chrY')
chromBoundaries <- c(0, cumsum(rle(ct$chrom)$lengths)) + 0.5
chromSizes <- diff(chromBoundaries)
chromLabelPos <- chromBoundaries[1:nChroms] +  chromSizes / 2
ct$binSizeCorrection <- env$BIN_SIZE / ct$binSize
obs <- apply(ct[,cellJs], 2, '*', ct$binSizeCorrection)

# clear prior images
PLOT_DIR <- file.path(env$PLOT_DIR, 'cell_plots')
dir.create(PLOT_DIR, showWarnings=FALSE)
unlink(file.path(PLOT_DIR, '*'))

# plot cell summed and ranked read counts
file <- paste('cellSums', 'png', sep=".")
file <- file.path(env$PLOT_DIR, file)
png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
cellSums <- colSums(ct[autosomes,cellJs])
order <- order(cellSums)
plot(cellKs, cellSums[order],
     xlab="Cell Rank", ylab='# Read Pairs (no dups)',
     main="Ranked Sample Depths", pch=16, col=cols[order])
graphics.off()

# set fit and plot parameters
toFit <- data.frame(
    x  = ct[autosomes,'gc_fraction'],
    x2 = ct[autosomes,'gc_fraction'] ^ 2,
    w  = ct[autosomes,'binSize'],
    y  = rowSums(obs[autosomes,])
)
toPredict <- data.frame(
    x  = ct[,'gc_fraction'],
    x2 = ct[,'gc_fraction'] ^ 2,
    y  = rowSums(obs)
)
gcPch <- 16
gcCex <- 0.4

# make a composite lorenz plot
plotLorenz <- function(main, type, d){
    file <- paste('lorenz', type, 'png', sep=".")
    file <- file.path(env$PLOT_DIR, file)
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
    plot(c(0,1), c(0,1), typ="l", xlab="Fraction of Genome",
         xlim=c(0,1), ylab="Fraction of Reads", ylim=c(0,1), main=main)
    for(k in cellKs){
        y <- d[autosomes,k]
        order <- order(y)     
        x <- ct[autosomes,'binSize'][order]    
        y <- y[order]
        x <- cumsum(as.numeric(x)) / sum(x)
        y <- cumsum(as.numeric(y)) / sum(y)
        lines(x, y, col=cols[k])
    }
    graphics.off()    
}
plotLorenz('Raw Read Counts', 'raw', obs)

# fit and plot each cell
obs_gc <- obs
obs_z  <- obs
plotCell <- function(k){
    
    # fit a curve to account for GC bias
    obs_k <- obs[,k]
    toFit$y <- obs_k[autosomes]
    fit <- loess(y ~ x, toFit, weights=toFit$w) # more prone to outlier influence
    #fit <- lm(y ~ x + x2, toFit, weights=toFit$w)
    
    # use fit to calculate inferred copy number
    toPredict$y <- obs_k
    exp0 <- predict(fit, toPredict)
    rpa  <- exp0 / 2     
    cn   <- obs_k / rpa
    obs_gc[,k] <<- cn

    # plot the GC fit
    sampleN <- sprintf("%03d", sampleNs[k])
    cellName <- cellNames[k]
    main <- paste(adjustedSampleNames[k], cellName,
                  format(cellSums[k],big.mark=",",scientific=FALSE), sep = " / ")
    file <- paste('gc_fit', sampleN, cellName, 'png', sep=".")
    file <- file.path(PLOT_DIR, file)
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")     
    plot(toPredict$x, toPredict$y, pch=gcPch, cex=gcCex, col="red3",
         xlab="Fraction GC", ylab="Bin Read Count", main=main) # sex chroms
    points(toFit$x, toFit$y,  pch=gcPch, cex=gcCex, col="blue") # autosomes
    points(toPredict$x, exp0, pch=gcPch, cex=gcCex, col="green3")
    graphics.off()     

    # plot the genome scatter plot
    file <- paste('bins', sampleN, cellName, 'png', sep=".")
    file <- file.path(PLOT_DIR, file)
    png(file, width = 8, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
    minY <- 0
    maxY <- 5
    cn <- ifelse(cn<minY, minY, cn)
    cn <- ifelse(cn>maxY, maxY, cn)
    plot(0, 0, typ="n",
         xlab="Bin Number", xlim=c(0,nBins+1),
         ylab="Inferred Copy Number", ylim=c(minY, maxY),
         main=main )
    abline(h=minY:maxY, col="grey")
    abline(v=chromBoundaries, col="grey")   
    points(1:nBins, cn, pch=16, cex=0.7)
    mtext(gsub('chr','',chroms), at = chromLabelPos, las=2)

    agg <- aggregate(cn, ct[,'chrom',drop=FALSE], function(v) {
        c(
            mean = mean(v),
            sd   = sd(v),
            q05  = quantile(v, 0.05),
            q95  = quantile(v, 0.95)
        )
    })
    agg <- agg[match(unique(ct$chrom), agg[[1]]),]
    names(agg) <- c('chrom','cn')
    for(i in 1:nrow(agg)){
        x   <- c(chromBoundaries[i], chromBoundaries[i+1])
        mu  <- agg[i,'cn'][,1]
        lines(x, rep(mu,2),       lwd=2, col="red3")        
        #sd2 <- agg[i,'cn'][,2] * 2
        #lines(x, rep(mu-sd2,2), lwd=1, col="blue")
        #lines(x, rep(mu+sd2,2), lwd=1, col="blue")
        lines(x, rep(agg[i,'cn'][,3],2), lwd=1, col="blue")
        lines(x, rep(agg[i,'cn'][,4],2), lwd=1, col="blue")
    }
    graphics.off()
    
    obs_z[,k] <<- unlist(sapply(agg$chrom, function(chrom){
        v <- cn[ct$chrom==chrom]
        (v - mean(v)) / sd(v) # z-score relative to the mean of one chrom for one sample
    }))
}
invisible(lapply(cellKs, plotCell))

# plot GC-corrected Lorenz
plotLorenz('GC-corrected Read Counts', 'gc', obs_gc)

# plot the genome scatter plot
medianZ <- apply(obs_z, 1, median)
file <- paste('bins', 'median_z_score', 'png', sep=".")
file <- file.path(env$PLOT_DIR, file)
png(file, width = 8, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
minZ <- -3
maxZ <-  3
medianZ <- ifelse(medianZ<minZ, minZ, medianZ)
medianZ <- ifelse(medianZ>maxZ, maxZ, medianZ)
plot(0, 0, typ="n",
     xlab="Bin Number", xlim=c(0,nBins+1),
     ylab="Median Z Score", ylim=c(minZ, maxZ))
abline(h=minZ:maxZ, col="grey")
abline(v=chromBoundaries, col="grey")   
points(1:nBins, medianZ, pch=16, cex=0.7)
mtext(gsub('chr','',chroms), at = chromLabelPos, las=2)
graphics.off()

# apply a correction for the WGA method bias at each bin, based on median Z
toFit <- data.frame(
    x  = medianZ[autosomes],
    w  = ct[autosomes,'binSize']
)
toPredict <- data.frame(
    x  = medianZ
)
plotCell2 <- function(k){
    
    # fit a curve to account for GC bias
    obs_k <- obs[,k]
    toFit$y <- obs_k[autosomes]
    fit <- loess(y ~ x, toFit, weights=toFit$w) # more prone to outlier influence
    #fit <- lm(y ~ x + x2, toFit, weights=toFit$w)
    
    # use fit to calculate inferred copy number
    toPredict$y <- obs_k
    exp0 <- predict(fit, toPredict)
    rpa  <- exp0 / 2     
    cn   <- obs_k / rpa

    # plot the GC fit
    sampleN <- sprintf("%03d", sampleNs[k])
    cellName <- cellNames[k]
    main <- paste(adjustedSampleNames[k], cellName,
                  format(cellSums[k],big.mark=",",scientific=FALSE), sep = " / ")
    file <- paste('z_fit', sampleN, cellName, 'png', sep=".")
    file <- file.path(PLOT_DIR, file)
    png(file, width = 4, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
    plot(medianZ, toPredict$y, pch=gcPch, cex=gcCex, col="red3",
         xlab="Median Z", ylab="Bin Read Count", main=main) # sex chroms
    points(toFit$x, toFit$y,  pch=gcPch, cex=gcCex, col="blue") # autosomes
    points(toPredict$x, exp0, pch=gcPch, cex=gcCex, col="green3")
    graphics.off()
    
    # plot the genome scatter plot
    file <- paste('bins2', sampleN, cellName, 'png', sep=".")
    file <- file.path(PLOT_DIR, file)
    png(file, width = 8, height = 4, units = "in", pointsize = 8, res = 600,  type = "cairo")
    minY <- 0
    maxY <- 5
    cn <- ifelse(cn<minY, minY, cn)
    cn <- ifelse(cn>maxY, maxY, cn)
    plot(0, 0, typ="n",
         xlab="Bin Number", xlim=c(0,nBins+1),
         ylab="Inferred Copy Number", ylim=c(minY, maxY),
         main=main )
    abline(h=minY:maxY, col="grey")
    abline(v=chromBoundaries, col="grey")   
    points(1:nBins, cn, pch=16, cex=0.7)
    mtext(gsub('chr','',chroms), at = chromLabelPos, las=2)

    agg <- aggregate(cn, ct[,'chrom',drop=FALSE], function(v) {
        c(
            mean = mean(v),
            sd   = sd(v),
            q05  = quantile(v, 0.05),
            q95  = quantile(v, 0.95)
        )
    })
    agg <- agg[match(unique(ct$chrom), agg[[1]]),]
    names(agg) <- c('chrom','cn')
    for(i in 1:nrow(agg)){
        x   <- c(chromBoundaries[i], chromBoundaries[i+1])
        mu  <- agg[i,'cn'][,1]
        lines(x, rep(mu,2),       lwd=2, col="red3")        
        #sd2 <- agg[i,'cn'][,2] * 2
        #lines(x, rep(mu-sd2,2), lwd=1, col="blue")
        #lines(x, rep(mu+sd2,2), lwd=1, col="blue")
        lines(x, rep(agg[i,'cn'][,3],2), lwd=1, col="blue")
        lines(x, rep(agg[i,'cn'][,4],2), lwd=1, col="blue")
    }
    graphics.off()

}
invisible(lapply(cellKs, plotCell2))



