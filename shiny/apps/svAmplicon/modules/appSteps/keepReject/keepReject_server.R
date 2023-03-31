#----------------------------------------------------------------------
# server components for the keepReject appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
keepRejectServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'keepReject'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    # settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
)

#----------------------------------------------------------------------
# get the samples and data to plot
#----------------------------------------------------------------------
sampleSet <- sampleSetServer("sampleSet", id)
samples <- samplesReactive(sampleSet)
amplicons <- ampliconsReactive(samples)
moleculeTypes <- moleculeTypesReactive(samples)

#----------------------------------------------------------------------
# construct the output tables and cascading reactives
#----------------------------------------------------------------------
ampliconsTable <- ampliconsTableServer(id, input, amplicons, selection = "single")
selectedAmplicons <- selectedAmpliconsReactive(amplicons, ampliconsTable)
ampliconMoleculeTypes <- ampliconMoleculeTypesReactive(moleculeTypes, selectedAmplicons)

#----------------------------------------------------------------------
# construct the output plots
#----------------------------------------------------------------------
moleculeMetadata <- moleculeMetadataReactive(selectedAmplicons, ampliconMoleculeTypes)
output$moleculeMetadata <- moleculeMetadataUI(moleculeMetadata)
output$moleculeQcPlot <- renderPlot({ moleculeQcPlot(moleculeMetadata) }) 

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    # settings$replace(bm$settings)
    updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
    if(!is.null(bm$outcomes)) {
        # outcomes <<- listToReactiveValues(bm$outcomes)
    }
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    # settings = settings$all_,
    outcomes = reactive({ list(
    ) }),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
