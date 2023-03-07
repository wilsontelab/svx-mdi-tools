#----------------------------------------------------------------------
# if present, apps/<appName>/server.R must contain a single function 
# called 'appServer', which is called in session context and thus has 
# implicit access to:
#   input, output, session objects
#   values returned from all app step modules
#----------------------------------------------------------------------
# if not needed, simply omit file server.R from your app
#----------------------------------------------------------------------

# objects instantiated here are available to all appStep modules in a session
sourceExternalScript("genomex-mdi-tools", "shared/global/utilities/sequence.R")

# settings presets
settingsPresets <- list(
    Fragile_Site = list(
        SV_Filters = list(
            Min_SV_Size = 10000,
            Max_Samples_With_SV = 1,
            Max_Source_Molecules = 1,
            Min_Map_Quality = 20,
            SV_Type = "Del"
        ),
        Capture_Filters = list(
            Min_Read_Count = 3,
            Target_Class = "both on target"
        ),
        Junction_Properties = list(
            Min_Insertion_Size = 3,
            Max_Insertion_Size = 13
        )
    )
)

# appServer function called after all modules are instantiated
appServer <- function(){


}
