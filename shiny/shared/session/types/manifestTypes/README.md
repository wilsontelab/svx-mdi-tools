---
title: Manifest Types
parent: Item Types
grand_parent: Stage 2 Apps
has_children: false
nav_order: 1
---

## Manifest Types

A **manifest** is a file provided by a Stage 1 pipeline or data provider
in a project zip that enumerates a related set of data samples.
Scripts in the **manifestTypes** folder:
- declare expected metadata in a specific type of manifest file
- provide handler functions used by sourceFileUpload and related modules

All manifestTypes, e.g., 'manifestTypes$xyz', must be a **list** with members:

```r
manifestTypes$xyz <- list(
    patterns = character(),
    load = function(file) ...,
    parse = function(data.frame) ...
)
```

where:

- **$patterns** = file suffixes for manifest files that match the manifestType
- **$load** = a function to read the manifest file from disk
- **$parse** = a function that takes the data frame from 'load' and returns a list of two data frames:
    - **$manifest** = all input rows, potentially with modified column names/data types
    - **$unique**   = one row per unique sample, often just a repeat of $manifest

In addition to any desired sample value columns, the parsed manifest data frames must 
ALWAYS have **Project**, **Sample_ID**, and **Description** columns - if necessary, create 
and fill those columns with a constant value (e.g., 'all_samples' or 'NA').
'Project' and 'Sample_ID' are concatenated to create unique sample identifiers, while 
'Description' is a short, human readable, sample name that is the default 
value for the editable sample name shown in the UI.

Data frame <code>$parsed$unique</code> may optionally have **Yield** and **Quality** columns, 
when appropriate (they will be displayed as NA if absent).
