---
title: Stage 2 app
parent: svAmplicon
has_children: false
nav_order: 20
published: true
---

## {{page.title}}

Results from the svAmplicon pipeline are mainly intended for 
examination and curation in the associated Shiny app.

Load the data package ending with `mdi.package.zip`
into the svAmplicon R Shiny app. The best way to run the app server
is by [installing the MDI Desktop](https://midataint.github.io/mdi-desktop-app/docs/installation),
which gives you several options for running the server on your local computer,
on your HPC server, or on a public web server.

Please see the instructions in the running app for details.
The general steps are to:
- load your desired pipeline data package(s)
- rename samples and assign them into experimental groups
- explore indvidual SV calls, keeping and rejecting individual reads based on quality assessments
- download final results tables

Additionally, the app provide a genome track browser for 
visualizing the distribution of junctions within amplicon spans.

Use the `Save your Work` link at the left to save a bookmark
of your work that you can load into the app server to restore
the page to the same state.
