#----------------------------------------------------------------------
# handle unique junction loading and filtering
#----------------------------------------------------------------------

# app-specific junction loading
svx_loadJunctions_app <- function(sourceId){
    readRDS(getSourceFilePath(sourceId, "junctionClustersFile")) 
}

# get a single junction to expand
svPore_getJunction <- function(x){
    svx_loadJunctions(x$sourceId)[clusterN == x$clusterN]
}
