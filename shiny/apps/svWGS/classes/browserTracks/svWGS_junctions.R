#----------------------------------------------------------------------
# handle unique junction loading and filtering
#----------------------------------------------------------------------

# app-specific junction loading
# here, coerce svWGS format toward svPore format (result is cached upstream)
svWGS_loadJunctions <- function(targetId){
    startSpinner(session, message = "loading jxns")
    jxns <- readRDS(getSourceFilePath(targetId, "structuralVariants")) 
    chroms <- readRDS(getSourceFilePath(targetId, "chromosomesFile"))
    samples_ <- read_yaml(getSourceFilePath(targetId, "metadata"))$SAMPLES
    samples_ <- unlist(strsplit(samples_, "\\s"))
    startSpinner(session, message = "loading jxns .")
    # NOTE: several derivative column values are set downstream by svx_loadJunctions
    # this function must set values whose input columns might differ between data sources
    jxns[, ":="(
        edgeType = svx_jxnType_altCodeToX(JXN_TYPE, "code"),
        chromIndex1 = unlist(chroms$chromIndex[CHROM_1]),
        chromIndex2 = unlist(chroms$chromIndex[CHROM_2]),
        strand1 = ifelse(SIDE_1 == "L", "+", "-"),
        strand2 = ifelse(SIDE_2 == "R", "+", "-"),
        insertSize = -MICROHOM_LEN,
        nSamples = N_SAMPLES,
        nInstances = N_TOTAL, 
        nSequenced = N_SPLITS,  
        mapQ = pmin(
            sapply(MAPQ_1, function(m) max(as.integer(strsplit(m, ",")[[1]]))),
            sapply(MAPQ_2, function(m) max(as.integer(strsplit(m, ",")[[1]])))
        ), 
        flankLength = pmin(FLANK_LEN1, FLANK_LEN2),
        nLinkedJunctions = N_CLUSTERED_JUNCTIONS
    )]
    startSpinner(session, message = "loading jxns ..")
    jxns[, samples := paste(
        samples_[sapply(samples_, function(x) .SD[[x]] > 0)], 
        collapse = ","
    ), by = .(SV_ID)]
    jxns[, ":="(
        node1 = getSignedNode(chroms$chromSizes, chromIndex1, POS_1, strand1),
        node2 = getSignedNode(chroms$chromSizes, chromIndex2, POS_2, strand2) 
    )]
    startSpinner(session, message = "loading jxns ...")
    jxns[, isCanonical := isCanonicalStrand_vectorized(node1, node2)]
    jxns[, ":="(
        cChromIndex1 = ifelse(isCanonical, chromIndex1, chromIndex2),
        cChromIndex2 = ifelse(isCanonical, chromIndex2, chromIndex1),
        cRefPos1     = ifelse(isCanonical, POS_1, POS_2),
        cRefPos2     = ifelse(isCanonical, POS_2, POS_1),
        cStrand1     = ifelse(ifelse(isCanonical, strand1, strand2) == "+", 1, -1),
        cStrand2     = ifelse(ifelse(isCanonical, strand2, strand1) == "+", 1, -1)
    )] 
    jxns
}

# get a single junction cluster for the object table
svWGS_getJunction <- function(x){
    svx_loadJunctions(x$targetId, svWGS_loadJunctions)[SV_ID == x$SV_ID]
}
