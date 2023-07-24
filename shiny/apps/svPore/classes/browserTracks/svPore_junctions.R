#----------------------------------------------------------------------
# handle unique junction loading and filtering
#----------------------------------------------------------------------

# app-specific junction loading
svPore_loadJunctions <- function(sourceId){
    readRDS(getSourceFilePath(sourceId, "junctionClustersFile")) 
}

# get a single junction to expand
svPore_getJunction <- function(x){
    svx_loadJunctions(x$sourceId, svPore_loadJunctions)[clusterN == x$clusterN]
}
