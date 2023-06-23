#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------
library(bit64)
svPoreCache <- new_dataCache('svPoreCache')

# appServer function called after all modules are instantiated
appServer <- function(){


}
