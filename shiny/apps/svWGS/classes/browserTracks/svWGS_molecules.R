#----------------------------------------------------------------------
# handle junction loading and filtering, for molecule plots, etc.
#----------------------------------------------------------------------
svWGSMoleculesCreate <- "asNeeded"

# load all molecules for a specific sourceId (not filtered yet)
svWGS_loadSourceMolecules <- function(sourceId){
    req(sourceId)
    startSpinner(session, message = "loading molecules")
    sessionCache$get(
        'svWGS_molecules', 
        key = sourceId, 
        permanent = TRUE, 
        from = "ram", 
        create = svWGSMoleculesCreate, 
        createFn = function(...) {
            startSpinner(session, message = "loading molecules .")
            molecules <- readRDS(getSourceFilePath(sourceId, "junctionMolecules")) 
            molecules[, c('chrom1', 'side1', 'pos1') := unpackNodeNames(NODE_1)]
            molecules[, c('chrom2', 'side2', 'pos2') := unpackNodeNames(NODE_2)]
            molecules[, isOuterClip := NODE_CLASS == SVX$nodeClasses$OUTER_CLIP]          
            molecules[, ":="(
                molKey = paste(SAMPLE, MOL_ID, sep = ":")
            )]
            setkey(molecules, SV_ID) 
            molecules
        }
    )$value
}

# pull the molecules for a specific junction
svWGS_loadMolecules <- function(sourceId, svId_){
    sessionCache$get(
        'svWGS_molecules', 
        keyObject = list(sourceId = sourceId, svId = svId_), 
        permanent = FALSE, 
        from = "ram", 
        create = svWGSMoleculesCreate, 
        createFn = function(...) {
            molecules <- svWGS_loadSourceMolecules(sourceId)
            startSpinner(session, message = "loading junction")
            molecules[svId_]
        }
    )$value
}
