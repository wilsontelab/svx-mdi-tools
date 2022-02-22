---
title: App Step Modules
parent: Shiny Modules
grand_parent: Stage 2 Apps
has_children: false
nav_order: 1
---

{% include table-of-contents.md %}

## Overview

**appSteps** share the properties of all MDI/Shiny modules
but have a specific purpose of defining one sequential app
step listed on the dashboard sidebar.

See the framework appStep modules for working examples.

## appStep UI function 

### Arguments

In addition to the standard Shiny 'id' argument, appStep UI functions must
take object 'options' as a second argument:

```r
# appStep module ui function, in myAppStep_ui.R
myAppStepUI <- function(id, options)
```

The options object is a nested R list with the module's configuration options,
as defined by the module's **module.yml** and the calling app's **config.yml** files.

### Element definition

appStep UI functions must create a single UI element as follows:

```r
# appStep module ui function, in myAppStep_ui.R
standardSequentialTabItem(
    title,
    leaderText,
    ... # the UI elements for the appStep's tabbed page
)  
```

## appStep server function

### Arguments

In addition to the standard Shiny 'id' argument, appStep server functions must
take the following additional arguments:

```r
# appStep module server function, in myAppStep_server.R
myAppStepServer <- function(id, options, bookmark, locks)
```

where:

- **options** = the same configuration object as passed to the UI function
- **bookmark** = access to the server bookmark reactive
- **locks** = access to the server locks reactive

### Return value

appStep server functions must return a specifically structured list
as a return object to allow proper visibility control for sequential
actions and to ensure that proper values are stored in bookmarks.

```r
# appStep module server function, in myAppStep_server.R
# return value
list(
    outcomes = list(   
        myOutcome = reactive()
    ),
    isReady = reactive({ getStepReadiness(...) }),
    ...
)
```

where:

- **outcomes** - a list of reactives for use by calling scripts and stored in bookmarks
- **isReady** - a reactive, or a function that returns a logical, to indicate whether the 
appStep has been completed; later steps will be masked until isReady() returns TRUE

The list may contain other objects or methods as needed, but they
will not be included in bookmarks unless they are a subset of the 'outcomes' key
and reactive.

### Special case: first app steps

Any module used as the _first_ step of an app, i.e., replacing the 
typical first 'sourceFileUpload' module, must also have the following
_additional_ members of its returned value list:

```r
# return value
list(
    outcomes = list(
        analysisSetName = reactive(...) # used for default naming of bookmark files
    ),
    loadSourceFile = function(incomingFile, suppressUnlink) ... # data passed from the universal launch page
)
```

Most apps are encouraged to use sourceFileUpload as the first
appStep module as it fulfills common requirements in well tested code.
