---
title: Analysis Types
parent: Item Types
grand_parent: Stage 2 Apps
has_children: false
nav_order: 2
---

## Analysis Types

An **analysis** is a discrete unit of computational work 
applied to a data set by a user within a running app. 
Scripts in the **analysisTypes** folder
define a specific way of analyzing data that can be:
- executed as a synchronous or asynchronous job
- potentially reused by many apps

Here's how to set an analysis type's parameters and handler functions for
use by runAnalyses and related modules.

All analysisTypes must provide, in folder **.../analysisTypes/\<typeGroup\>/\<analysisType\>**

- a **config.yml** file, with members (all are optional):
    - **name** = display name for the analysisType; default = \<analysisType\>
    - **jobType** = how/where the job should be executed; see job_execution.R; default = promise
    - **options** = job execution options, analogous to module settings; see existing prototypes. There are four UI columns for options display; use 'type: empty' for a blank position
    - **packages** = R or Bioconductor packages used by the analysisType; installed but not attached to running server process
    - **classes** = optional data classes used by the analysisType
    - **modules** = optional ui modules used by the analysisType  
    
- **\<analysisType\>_methods.R** file, with the following S3 generic functions:
    - **setJobParameters.\<analysisType\>** = convert reactives to static values for job promises
    - **executeJob.\<analysisType\>** = do the actual work of the job, potentially within a promise
    - **load.\<analysisType\>** = load the output of executeJob() for use by viewer modules

- a **results viewer Shiny module** used to populate the 'results' uiOutput in module viewResults via files:
    - **\<analysisType\>_ui.R**
    - **\<analysisType\>_server.R**
