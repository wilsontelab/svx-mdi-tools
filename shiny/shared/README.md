---
title: Shared Files
parent: Stage 2 Apps
has_children: true
nav_order: 2
---

{% include table-of-contents.md %}

## Shared folder structure

This **shiny/shared** folder has code components that may be 
called by multiple analysis apps.

The **global** folder has scripts that define utility functions
that do not depend on user session data.
They are sourced prior to initialization of a specific
app and must not depend on <code>app</code> or other session 
variables, unless they are passed as function arguments.

The **session** folder has scripts that define functions that are
potentially specific to a user session or needed to assemble and serve 
an app. They are sourced after a user commits to a specific app and 
therefore have <code>app</code> and other session variables 
implicitly in their scope. 

Scripts in both the **global** and **session** folders are sourced 
into the R environment of every user session, i.e.,
in <code>sessionEnv</code>. This enables hotfixes and rolling updates
without server restart, only a page reload.

The **static** folder has fixed content for populating 
pages with text, mostly via markdown rendered in R with
<code>includeMarkdown(file.path('static/xxx.md'))</code>.

## Recurring code elements

**classes** are R Shiny S3 classes that define reusable data objects
for use by developers writing MDI code. Their purpose
is to encapsulate methods and other common logic. Classes are
never used to access or create the UI.

**modules** are R Shiny modules that define reusable UI components
and associated server logic for app steps, widgets, etc.

- **appControl** - control overall framework behavior
- **appSteps** - common analysis steps used by apps, e.g., to load files
- **widgets** - UI components to embed on pages, e.g., plot boxes

**packages** holds files that declare and load the various R packages
used to implement the framework and handle data.

**types** define data classes that help structure app contents.

- **analysisTypes** - common ways that data may be analyzed across apps
- **manifestTypes** - the file formats by which data providers deliver sample metadata  

**ui** scripts offer functions to help assemble page UIs. 

**utilities** scripts offer functions to get or set various data values.

>**modules**, **types**, **ui**, and **utilities** folders can also
be created within app folders, to be loaded only with that specific 
app. However, it is often desirable to abstract such files to be reusable 
in other apps.

## Shared code versioning

Code components in the shared folder are scoped to the
suite of tools, so the relevant versioning is to create appropriate release 
tags for the tool suite repository that reflect the code changes 
in your shared components. Doing so will help people who want to use 
shared components from your tool suite when building their own apps.

All apps within a tool suite are expected to have been adjusted for 
any changes made within its shared components. You should advance
the version of your apps to reflect this change, most likely as a patch.
