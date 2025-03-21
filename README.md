# Wilson lab tools for structural variant analysis

The [Michigan Data Interface](https://midataint.github.io/) (MDI) 
is a framework for developing, installing and running 
Stage 1 HPC **pipelines** and Stage 2 interactive web applications 
(i.e., **apps**) in a standardized design interface.

This **svx-mdi-tools** repository contains pipelines and apps
for genome structural variant (SV) analysis by different library strategies
from the 
[Thomas Wilson laboratory](https://wilsonte-umich.github.io)
at the University of Michigan.

## Available tools

Pipelines and apps in stable release are:
- **svCapture** = find SVs in whole-genome capture sequencing ([publication](https://academic.oup.com/nargab/article/5/2/lqad042/7157526))
- **svAmplicon** = find SVs in targeted amplicon sequences
- **svWGS** = find SVs in whole-genome sequencing data (without capture)

Pipelines and apps in beta, which have functioning code under development
but no official releases yet, are:
- **svPore** = find SVs in nanopore long-read whole-genome sequencing data 
- **scCNV** = find CNVs in single-cell whole-genome amplification sequencing

Pipelines and apps in alpha, with exploratory code that is not considered stable:
- **svDJ** = analyze normal and aberrant V(D)J recombination events

# Installation and usage

Please see the [svx-mdi-tools documentation site](https://wilsontelab.github.io/svx-mdi-tools)
for detailed information on suite installation and the usage of specific pipeline and apps.
