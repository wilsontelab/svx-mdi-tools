#----------------------------------------------------------------------
# constants used to set assembly plot appearance and layout
#----------------------------------------------------------------------
CONSTANTS$assemblyPlots <- list(
    linesPerInch = 8.571429,
    fontSize = 7,
    nullMar = 0.5,
    titleMar = 2.1,
    titleLegendMar = 2.1,
    stdAxisMar = 4.1,
    untitledAxisMar = 2.1,
    noMar = 0.1
)
CONSTANTS$assemblyPlots$svFrequencies <- list(
    groupWidth = 0.3,    
    insideHeight = 1,
    mar = c(CONSTANTS$assemblyPlots$titleLegendMar, 8.1, CONSTANTS$assemblyPlots$titleMar, CONSTANTS$assemblyPlots$nullMar)
)
CONSTANTS$assemblyPlots$microhomology <- list(
    insideWidth  = 1.5,
    insideHeight = 1,
    mar = c(CONSTANTS$assemblyPlots$stdAxisMar, CONSTANTS$assemblyPlots$stdAxisMar, CONSTANTS$assemblyPlots$titleLegendMar, CONSTANTS$assemblyPlots$nullMar)
)
CONSTANTS$assemblyPlots$endpoints <- list(
    insideWidth  = 2.5,
    insideHeight = 0.75,
    mar = c(CONSTANTS$assemblyPlots$stdAxisMar, CONSTANTS$assemblyPlots$stdAxisMar, CONSTANTS$assemblyPlots$titleLegendMar, CONSTANTS$assemblyPlots$nullMar)
)
CONSTANTS$assemblyPlots$svSizes <- list(
    insideWidth  = 1.5,
    insideHeight = 1,
    mar = c(CONSTANTS$assemblyPlots$stdAxisMar, CONSTANTS$assemblyPlots$stdAxisMar, CONSTANTS$assemblyPlots$titleLegendMar, CONSTANTS$assemblyPlots$nullMar)
)
