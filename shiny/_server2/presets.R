
# common function for resetting filters
resetCommonFilters <- function(){
    updateSelectInput(session, 'captureTarget', selected='-')
    updateRadioButtons(session, 'duplexFilter', selected='all')
    updateRadioButtons(session, 'netReadPairFilter', selected='all')
    updateRadioButtons(session, 'moleculeCountFilter', selected='all')
    updateRadioButtons(session, 'sampleCountFilter', selected='all')
    updateRadioButtons(session, 'fracProperEnds', selected='all')
    updateCheckboxGroupInput(session, 'targetTypeFilter', selected=defaultTargetPairTypes)        
    updateRadioButtons(session, 'svTypeFilter', selected='all')        
    updateRadioButtons(session, 'sequencedFilter', selected='all')
}

# presets for useful filter combinatinos
observeEvent(input$candidateSV, {
    resetCommonFilters()
    updateRadioButtons(session, 'netReadPairFilter', selected='>3')
    updateRadioButtons(session, 'fracProperEnds', selected='<1')
    updateCheckboxGroupInput(session, 'targetTypeFilter', selected=c('TT','TA'))
})
observeEvent(input$ligationArtifact, {
    resetCommonFilters()
    updateRadioButtons(session, 'duplexFilter', selected='yes')
    updateRadioButtons(session, 'moleculeCountFilter', selected='1')
    updateRadioButtons(session, 'fracProperEnds', selected='<1')
    updateCheckboxGroupInput(session, 'targetTypeFilter', selected=c('t-'))        
})
observeEvent(input$chimericPCR, {
    resetCommonFilters()
    updateRadioButtons(session, 'duplexFilter', selected='no')
    updateRadioButtons(session, 'moleculeCountFilter', selected='1')
    updateCheckboxGroupInput(session, 'targetTypeFilter', selected=c('tt'))          
})

