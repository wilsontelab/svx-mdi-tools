# get environment variables
env <- as.list(Sys.getenv())

# write a table to a file
file <- paste(env$DATA_FILE_PREFIX, "mtcars", "csv", sep = ".")
write.csv(mtcars, file, quote = FALSE, row.names = FALSE)

# make a plot
file <- paste(env$DATA_FILE_PREFIX, "test", "png", sep = ".")
png(file, width = 3, height = 3, units = "in", 
    pointsize = 8, res = 600, type = "cairo")
plot(1:10)
dev.off()
