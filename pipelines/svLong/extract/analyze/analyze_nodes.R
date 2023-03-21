# sort and group molecules based on the SVs they carry

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("  initializing")
library(data.table)
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'OUTPUT_DIR',
        'DATA_NAME',
        'NODES_FILE'
    ),
    integer = c(
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# set some options
setDTthreads(env$N_CPU)
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
# nodeClasses <- c(
#     GAP   = 0, # SV evidence type codes, i.e. node classes
#     SPLIT = 1    
# )
#=====================================================================================


#=====================================================================================
nodes <- fread(
    env$NODES_FILE,
    col.names = c(
        "QNAME", # PAF fields
        "NODE1",
        "NODE2",
        "EDGE_TYPE",
        "MAPQ",
        "SV_SIZE",
        "INSERT_SIZE",
        "N_STRANDS"
    ),
    colClasses = c(
        "character", # PAF fields
        "integer",
        "integer",
        "character",
        "integer",
        "integer",
        "integer",
        "integer"
    )
)
str(nodes)

molecules <- nodes[, .(
    PATH_TYPE = paste(EDGE_TYPE, collapse = ""),
    N_SEGMENTS = .N,
    MIN_MAPQ = min(MAPQ),
    MAX_MAPQ = max(MAPQ),
    MIN_SV_SIZE = min(SV_SIZE),
    MAX_SV_SIZE = max(SV_SIZE),
    N_STRANDS = N_STRANDS[1],
    NODE_PATH = list(NODE1, NODE2[.N])
), by = "QNAME"]
str(molecules)
str(molecules[N_SEGMENTS > 1])

# find duplex paths sequenced as two molecules (presumably rare)

# characterize individual SV junctions, count molecules

# characterize composite SV paths, counts molecules

# above required strand-aware comparison, similar to NODE_PATH1 == reverse(-NODE_PATH2)



#=====================================================================================
