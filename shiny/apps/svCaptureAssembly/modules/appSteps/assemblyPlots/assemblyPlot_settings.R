#----------------------------------------------------------------------
# settings exposed for user customization of assembly plots
#----------------------------------------------------------------------
assemblyPlotFrameSettings <- list(
    Plot = list(
        Width_Inches = list( # allow size and title overrides of automated values
            type = "textInput",
            value = "auto"
        ),
        Height_Inches = list(
            type = "textInput",
            value = "auto"
        ),         
        Title = list(
            type = "textInput",
            value = ""
        )
    )
)
assemblyPlotSettings <- list(
    svFrequencies = c(assemblyPlotFrameSettings, list(
        Groups = list(
            Group_Width_Inches = list(
                type = "numericInput",
                value = CONSTANTS$assemblyPlots$svFrequencies$groupWidth,
                min = 0.1,
                max = 1,
                step = 0.05
            )
        )
    )),
    microhomology = c(assemblyPlotFrameSettings, list(
        Limits = list(
            Min_X_Value = list(
                type = "numericInput",
                value = -10,
                min = -50, 
                max = 0,
                step = 5
            ),
            Max_X_Value = list(
                type = "numericInput",
                value = 15,
                min = 0, 
                max = 50,
                step = 5
            )
        )
    )),
    endpoints = c(assemblyPlotFrameSettings, list(
        Bins = list(
            Bin_Size = list(
                type = "numericInput",
                value = 10000,
                min = 100, 
                max = 100000,
                step = 100
            )
        )
    )),
    svSizes = c(assemblyPlotFrameSettings, list(
        Bins = list(
            Bin_Per_Log = list(
                type = "numericInput",
                value = 10,
                min = 1, 
                max = 100,
                step = 1
            ),
            Min_SV_Size = list(
                type = "numericInput",
                value = 10000,
                min = 1000, 
                max = 100000,
                step = 1000
            ),
            Max_SV_Size = list(
                type = "numericInput",
                value = 1000000,
                min = 100000, 
                max = 10000000,
                step = 100000
            )
        )
    ))
)
