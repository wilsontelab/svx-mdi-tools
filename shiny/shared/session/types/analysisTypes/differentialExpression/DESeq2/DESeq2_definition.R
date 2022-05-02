#----------------------------------------------------------------------
# differential expression analysis by DESeq2
#----------------------------------------------------------------------

# setJobParameters generates the jobParameters object sent to executeJob
# a main purpose is to convert reactives to static variables for use by promises
setJobParameters.DESeq2 <- function(job){
    job$parameters <- reactiveToStatic_commmon()
    job
}

# executeJob does the actual analysis work
#   called by runJobXXX family of functions, often within a promise
#   if jobType is not 'immediate', must not depend on any reactives
executeJob.DESeq2 <- function(job){

    # load the count data
    composite <- getSampleCountData(job)

    # set the proper design, as per "analyze by"
    design <- if(job$schema$Analyze_By == "Category1") ~Category1
        else if (job$schema$Analyze_By == "Category2") ~Category2
        else stop(safeError("unknown Analyze By"))

    # run DESeq2
    reportProgress('running DESeq2')
    tryCatchJob({
        #installAndLoad_Bioconductor("DESeq2")
        dds <- DESeq2::DESeqDataSetFromMatrix(
            countData = composite$countData,
            colData   = composite$colData,
            design    = design
        )
        dds <- DESeq2::DESeq(dds)
        list(input = composite, output = DESeq2::results(dds, tidy = TRUE))
        ## or to shrink log fold changes association with condition:
        #res <- lfcShrink(dds, coef="condition_trt_vs_untrt", type="apeglm")
        # DESeq2 has a ton of overhead and may be causing memory leaks in shiny
    }, job$schemaId) 
}

# load the results of an analysis for use in viewResults or related ui modules
loadJobOutput.DESeq2 <- function(job){
    job$output <- readRDS( getJobRdsFile(job$schemaId) )
    job
}
