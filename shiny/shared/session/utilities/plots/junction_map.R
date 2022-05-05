# ----------------------------------------------------------------------
# colored-image representation of all molecules supporting and SV junction
# ----------------------------------------------------------------------
mapSettingsServer <- function(id, moduleOptions){
    globalSettingsDir <- file.path(app$sources$suiteGlobalDir, "settings")
    settingsServer( 
        id = 'mapSettings',
        parentId = id,
        templates = file.path(globalSettingsDir, "junction_map.yml"),
        fade = FALSE,
        title = "Base Map Settings",
        immediate = TRUE
    )    
}
junctionMapUI <- function(ns, width){
    box(
        width = width,                
        collapsible = TRUE,
        collapsed = FALSE,
        title = tagList(
            "Junction Evidence Map",
            settingsUI(ns('mapSettings'))
        ),
        imageOutput(ns("junctionMapImage"), inline = TRUE)
    )  
}
junctionMapServer <- function(output, junctionMap, mapSettings){
    output$junctionMapImage <- renderImage({
        map <- junctionMap()
        req(map)
        pngFile <- file.path(sessionDirectory, "junctionMapImage.png")
        pixels <- mapSettings$get("Map_Settings", "Pixels_Per_Base")
        suppressWarnings(
            imager::as.cimg(map$image[map$usedPos, , ]) %>% 
            expandImg(h = pixels, v = pixels) %>%
            imager::save.image(pngFile)        
        )
        list(src = pngFile)
    }, deleteFile = FALSE)
}
