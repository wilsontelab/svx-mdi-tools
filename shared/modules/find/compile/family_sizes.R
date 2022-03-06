# get arguments
env <- as.list(Sys.getenv(c('EXTRACT_PREFIX', 'COMPILE_PREFIX', 'PLOT_PREFIX'), names = TRUE))

# constants
IS_PROPER <- 'P' # proper and anomalous molecule codes
IS_SV     <- 'V'
IS_ORPHAN <- 'O'

# load the table of strand depths
d <- read.table(file = "stdin", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
groupCol  <- 3 # count1
countCols <- 4:ncol(d)
depths    <- 5:ncol(d) - 4

# declare the molecule types to tabulate and plot
df <- data.frame(
    targetClass = c('TT',           'TT',           'tt',           '--',         '--'),
    molClass    = c(IS_PROPER,      IS_SV,          IS_SV,          IS_PROPER,    IS_SV),
    class       = c('Proper',       'Del/Dup/Inv',  'Artifact',     'Proper',     'Artifact'),
    total   = rep(0, 5),
    duplex  = rep(0, 5),
    duplex3 = rep(0, 5),
    color      = c("green4", "purple", "red3", "blue", "gray50"),
    plotDuplex = c(TRUE, TRUE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
)
singletonDepth <- list()
duplexDepth    <- list()

# sum counts over all job threads that created the data
for(i in 1:nrow(df)){ # nolint
    
    # aggregate the data over all extract threads (previously sorted)
    rows <- d$targetClass == df[i, 'targetClass'] & d$molClass == df[i, 'molClass']
    dm <- d[rows, countCols]
    dt <- rowsum(dm, d[rows, groupCol])
    
    # tabulate molecule counts by plexity
    dx  <- dt[2:nrow(dt), 2:ncol(dt)]
    dx3 <- dt[4:nrow(dt), 4:ncol(dt)]
    df[i, 'total']   <- sum(dt)
    df[i, 'duplex']  <- sum(dx)
    df[i, 'duplex3'] <- sum(dx3)    
    
    # make a table of strand depths by reference molecule type
    key <- paste(df[i, 'targetClass'], df[i, 'molClass'], sep = "")
    sd  <- t(dt)[2:nrow(dt), 1] + dt[2:nrow(dt), 1]
    ssd <- sum(sd)
    singletonDepth[[key]] <- if(ssd > 0) sd / ssd else NA
    xd  <- rowSums(dx) + colSums(dx)
    sxd <- sum(xd)
    duplexDepth[[key]]  <- if(sxd > 0) xd / sxd else NA    
}

# print the composite table of counts and rates
df$singleton <- df$total - df$duplex
df$duplexRate <- df$duplex / df$total
print(df[, c('targetClass', 'class', 'total', 'duplex', 'duplex3', 'singleton', 'duplexRate')])

# make strand depth plots by reference molecule type
png(filename = paste(env$PLOT_PREFIX, 'family_sizes', 'png', sep = "."),
    width = 4, height = 4, units = "in", pointsize = 10,
    res = 600, type = "cairo")
plot(0, 0, type = "n", 
     xlab = "Strand Family Size", ylab = "Frequency", 
     xlim = c(1, 50), ylim = c(0, 0.15))
duplexLty    <- 1
singletonLty <- 2
for(i in 1:nrow(df)){ # nolint
    key <- paste(df[i, 'targetClass'], df[i, 'molClass'], sep = "")
    if(df[i, 'plotDuplex']){
        if(length(depths) == length(duplexDepth[[key]])){
            lines(depths, duplexDepth[[key]], lty = duplexLty, col = df[i, 'color'], lwd = 2)
        }
    }
    if(length(depths) == length(singletonDepth[[key]])){
        lines(depths, singletonDepth[[key]], lty = singletonLty, col = df[i, 'color'], lwd = 2)
    }
}
legend("topright", paste(df$targetClass, df$class, sep = " "), col = df$color, pch = 19, cex = 0.8)
