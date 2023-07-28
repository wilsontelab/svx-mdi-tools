---
title: "svx-mdi-tools"
has_children: false
nav_order: 0
---

{% include mdi-project-overview.md %} 

This **svx-mdi-tools** repository contains pipelines and apps
for genome structural variant (SV) analysis by different library strategies
from the 
[Thomas Wilson laboratory](https://wilsonte-umich.github.io)
at the University of Michigan.

{% include figure.html file="SV_logic.png" %}

### Available tools

Pipelines and apps in stable release are:
- **svCapture** = find SVs in whole-genome capture sequencing ([publication](https://academic.oup.com/nargab/article/5/2/lqad042/7157526))

Pipelines and apps in beta, which have functioning code under development
but no official releases yet, are:
- **svWGS** = find SVs in whole-genome sequencing data (without capture)
- **svPore** = find SVs in nanopore long-read whole-genome sequencing data 
- **scCNV** = find CNVs in single-cell whole-genome amplification sequencing

Pipelines and apps in alpha, with exploratory code that is not considered stable:
- **svAmplicon** = find SVs in targeted amplicon sequences
- **svDJ** = analyze normal and aberrant V(D)J recombination events

{% include mdi-project-documentation.md %}
