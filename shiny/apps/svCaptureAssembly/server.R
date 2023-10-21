#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------
assemblyCache <- new_dataCache('assemblyCache')
loadingPlotSet <- NULL

# names of the columns used internally for plotting, etc.
# all other columns were declared by user and are offered for grouping and selection
denominatorColumns <- c("kbp_on_target","kbp_on_target_effective","kbp_on_target_filtered","kbp_on_target_filtered_effective")
svTypeColumns <- c("deletion","duplication","inversion","translocation")
internalUseSampleColumns <- c(
    "project","sample","sample_id",
    "genomesDir","genome","cram","minMapq","targetsBed","regionPadding",
    "sumTargetLens","N50",
    "n_source_molecules","n_on_target","n_on_target_filtered",
    denominatorColumns,
    svTypeColumns,
    "coverage"
)

# appServer function called after all modules are instantiated
appServer <- function(){


}
