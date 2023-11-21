# ----------------------------------------------------------------------
# text representation of all an SV junction consensus sequence vs. genome references
# ----------------------------------------------------------------------

# functions used in a standard app page, e.g., in svCapture
alignmentSettingsServer <- function(id, moduleOptions){
    globalSettingsDir <- file.path(app$sources$suiteGlobalDir, "settings")
    settingsServer(
        id = 'alignmentSettings',
        parentId = id,
        templates = file.path(globalSettingsDir, "junction_alignment.yml"),
        fade = FALSE,
        title = "Junction Alignment Settings",
        immediate = TRUE
    )
}
junctionAlignmentUI <- function(ns, width){
    box(
        width = width,
        collapsible = TRUE,
        collapsed = FALSE,
        title = tagList(
            "Junction Alignment to Reference Genome",
            settingsUI(ns('alignmentSettings'))
        ),
        tags$style(HTML("
            .junction { background-color: #ddd;}
            .alignment .base_A { color: green; }
            .alignment .base_C { color: blue; }
            .alignment .base_G { color: brown; }
            .alignment .base_T { color: red; }
            .referenceGenome { font-style: oblique; }
        ")),
        div(
            htmlOutput(ns("junctionAlignment")),
            style = "font-family: monospace; font-weight: bold; width: 100; overflow: auto;"
        )
    )  
}
junctionAlignmentServer <- function(output, junctionMap, alignmentSettings){
    output$junctionAlignment <- renderText({
        getJunctionAlignment(
            junctionMap(),
            alignmentSettings$get("Alignment_Settings", "Bases_Per_Line"),
            alignmentSettings$get("Alignment_Settings", "Display_Mode")
        )
    })
}

# alternative to fill app$browser$expansionUI when reacting to browser track click
junctionAlignmentTrackExpansionUI <- function(junctionMap, Bases_Per_Line, Alignment_Mode){ 
    req(junctionMap)    
    tags$div(
        tags$style(HTML("
            .junction { background-color: #ddd;}
            .alignment .base_A { color: green; }
            .alignment .base_C { color: blue; }
            .alignment .base_G { color: brown; }
            .alignment .base_T { color: red; }
            .referenceGenome { font-style: oblique; }
        ")),
        div(
            style = "font-family: monospace; font-weight: bold; width: 100; overflow: auto; margin-top: 10px;",
            HTML(getJunctionAlignment(
                junctionMap,
                Bases_Per_Line,
                Alignment_Mode
            ))
        )
    )
} 
