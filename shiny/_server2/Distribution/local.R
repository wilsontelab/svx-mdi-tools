
#----------------------------------------------------------------------
# local.R launches the web tool on a desktop or laptop via localhost
#----------------------------------------------------------------------
# you just need to edit the proper file paths for your computer, below
#----------------------------------------------------------------------

# get the script path
actionsPath <- dirname(sys.frame(1)$ofile)
shinyPath <- paste(actionsPath, '..', '_shiny', sep="/")

# set required environment variables
Sys.setenv(
    PIPELINE_NAME   = "svCapture",
    SERVER_SUBTYPE  = "Distribution",

    #-----------------------------------------------
    # EDIT HERE: specific data paths for local mode
    #-----------------------------------------------
    DATA_PATH       = "Z:\\data\\scCNV",  # required
    GENOMES_DIR     = "Z:\\data\\genomes" # can be "", but you won't see gene annotations
    #-----------------------------------------------
)

# call app.R
source(paste(shinyPath, 'local.R', sep="/"))

