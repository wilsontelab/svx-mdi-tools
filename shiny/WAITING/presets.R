# common function for resetting filters
resetCommonFilters <- function(session){
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
