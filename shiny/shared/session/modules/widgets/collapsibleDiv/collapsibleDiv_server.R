#----------------------------------------------------------------------
# server components for the collapsibleDiv widget module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
collapsibleDivServer <- function(id) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'collapsibleDiv'
# settings <- activateMdiHeaderLinks( # uncomment as needed
#     session,
#     url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
#     dir = getAppStepDir(module), # for terminal emulator
#     envir = environment(), # for R console
#     baseDirs = getAppStepDir(module), # for code viewer/editor
#     settings = id, # for step-level settings
#     immediate = TRUE # plus any other arguments passed to settingsServer()
# )

#----------------------------------------------------------------------
# set return value, typically NULL or a list of reactives
#----------------------------------------------------------------------
NULL

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------
