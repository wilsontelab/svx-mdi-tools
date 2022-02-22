---
title: App-Specific Files
parent: Stage 2 Apps
has_children: false
nav_order: 1
---

{% include table-of-contents.md %}

## App folders

Each folder in **.../shiny/apps** has files that define 
the behavior of a single app (i.e., one Shiny data analysis interface). 

At the root level of the app, these files include
(see comments within for additional instructions):

- **config.yml**  - names the app and describes its structure
- **server.R**    - contains the function 'appServer'
- **overview.md** - text used to describe the app on its splash page

The name of the app's folder must match the name provided in config.yml.

Optionally, you can organize additional app scripts into the
following sub-folders, which will be loaded along
with config.yml and server.R when the app loads
(see other pages for a description of files in those folders):

- modules
- types
- ui
- utilities 

## Quick start

Start by copying and modifying the '_template' app, which provides 
a working boilerplate for all required code. Read the comments and go.

## Stage 2 versioning

### App versions

Individual app versioning is optional but recommended as it will
be displayed in the web page UI and help users
access legacy versions of your code to analyze their data according
to some previous standard.

To track app versions, add a semantic version
key to config.yml and update it prior to committing new app code. 
It is not necessary to create Git tags for app versions.

```yml
# shiny/apps/<app>/config.yml
name: Example
description: "Example of descriptive text"
version: v0.0.0
```

### External suite versions

If your app uses code modules from other tool suites, you may
wish to specify the required versions of those external suites.
This is useful if you don't wish to adjust your app to account for a
breaking change made in an external tool suite.  Declare such version
requirements as follows, replacing 'suiteName' with the name of the
external tool suite.

```yml
# shiny/apps/<app>/config.yml
name: Example
description: "Example of descriptive text"
version: v0.0.0
suiteVersions: # <<< add this section <<<
    suiteName: v0.0.0
```

If you do not provide a version for an external tool suite,
the latest version will be used.

If you only use app code from within your own tool suite, the 
suiteVersions dictionary can be omitted.
