---
title: R Packages
parent: Shared Files
grand_parent: Stage 2 Apps
has_children: false
nav_order: 1
---

## Global R package dependencies

Files in **global/packages** declare R packages 
that are used by one or more apps in a suite. They are read by 
MDI::install() to discover the R packages to install.

## Component R package dependencies

Any app, module, or other component can also declare its need
for a specific R package in its **module.yml** or **config.yml** file
(there is no harm in relisting a package already listed 
elsewhere - it will only be installed once):

```yml
packages: 
   R:  
    - xxx
    - yyy
   Bioconductor: null
```

## Installing vs. attaching

Packages declared in a tool suite are _installed_ but are never _attached_ using the <code>library()</code> function. Thus, you must either:

- call functions in full syntax, e.g., <code>package::function()</code> (**preferred**)
- call <code>library(package)</code> yourself (**discouraged**, you could break the framework)

In contrast, packages listed in the 'mdi-apps-framework' repository are attached and their functions can be called directly, e.g., Shiny's <code>observeEvent()</code>.
