# ----------------------------------------------------------------------
# text representation of all an SV junction consensus sequence vs. genome references
# ----------------------------------------------------------------------
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
