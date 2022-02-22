---
title: Shiny Modules
parent: Stage 2 Apps
has_children: true
nav_order: 3
---

{% include table-of-contents.md %}

## Module structure

**modules** are 
[R Shiny modules](https://shiny.rstudio.com/articles/modules.html)
that define reusable UI components.
A module folder carries scripts that define its UI elements and 
associated server code logic.

By convention, MDI module scripts follow the naming convention:

- **module_ui.R** - carries the function that returns the module's UI elements
- **module_server.R** - carries the function that defines the module's actions
- **module_utilities.R** - additional support functions for the module

Like any R function, module server functions return a single
object that can be used by the calling code. Typical return values
are a list of reactive objects and method functions to be applied to the  
component or its data.

## Module types

Modules have different roles in MDI apps:

- **appSteps** - sequential actions comprising the leftmost dashboard 
tabs of the MDI web page

- **widgets** - UI elements, panels, plots, etc. that
might be placed once or many times by an app

## Module script templates

The following code blocks show the structure of a module's UI and server scripts. 
Replace "myModule" with the name of your module. 
See existing framework modules for extended examples.

Module UI functions must use the Shiny <code>ns()</code> function to wrap
UI element names, which causes them to have properly nested and addressable
names.

```r
# module ui function, in myModule_ui.R
myModuleUI <- function(id, ...) {
    # initialize namespace
    ns <- NS(id) 
    # add Shiny UI elements
    textInput(ns('inputName'))
    textOutput(ns('outputName'))
}
```

UI elements are then addressable in the module's 
<code>input</code> object using the name provided in the UI function.

```r
# module server function, in myModule_server.R
myModuleServer <- function(id, ...) {
    moduleServer(id, function(input, output, session) {
#----------------------------------------------------------------------
# add Shiny server reactive actions here: observe, render, etc.
output$outputName <- renderText({
    input$inputName
})
# module return value (if none, use NULL instead of a list)
list(
    myValue = reactiveVal(input$inputName),
    myMethod = function(...) {}
)
#----------------------------------------------------------------------
})}
```

## Using a module

The following code blocks show how to use your module in calling scripts. 
It is common to nest modules within modules - the examples show how to 
ensure that names are assembled and accessed correctly.

```r
# within any UI script, e.g., parentModule_ui.R
myModuleUI(ns('moduleInstanceId'), ...)
```

Be sure that the name/id of the module matches between the UI and server calls.

```r
# within any server script, e.g., parentModule_server.R
# use id, NOT ns(id), within the server function!
myObject <- myModuleServer('moduleInstanceId', ...)
output$xyz <- renderText({
    myObject$myValue() # <<< accesses the modules reactive return value
})
```
