---
title: Stage 2 app
parent: svCapture
has_children: false
nav_order: 20
published: true
---

## {{page.title}}

Results from the svCapture pipeline can be examined directly
and used in other downstream applications via the VCF formatted output
and other files.

You can also load the data package ending with `mdi.package.zip`
into the svCapture (for the `find` action) 
or svCaptureAssembly (for the `assemble` action) 
R Shiny apps. The best way to run the app server
is by [installing the MDI Desktop](https://midataint.github.io/mdi-desktop-app/docs/installation),
which gives you several options for running the server on your local computer,
on your HPC server, or on a public web server.

Please see the instructions in the running app for details.
The general steps are to:
- load your desired pipeline data package(s)
- rename samples and assign them into experimental groups
- examine library QC
- explore indvidual SV calls
- make aggregate plots of SV formation rate and microhomologies

Common processes in many app steps are to:
- select the samples you are working with
- set filters and other page options under the 'gear' icon
- use the download links to export data table and figure images

Use the `Save your Work` link at the left to save a bookmark
of your work that you can load into the app server to restore
the page to the same state.
