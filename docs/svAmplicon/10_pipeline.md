---
title: Stage 1 pipeline
parent: svAmplicon
has_children: false
nav_order: 10
published: true
---

## {{page.title}}

svAmplicon data analysis first entails execution of a pipeline with primary actions:
- **align** = align input fastq files to the reference genome
- **extract** = scan aligned reads to identify amplicon boundaries and nominate unique reads with structural variants (SVs)

## Pipeline inputs

The primary inputs to the svAmplicon pipeline are:
- genome files as described under [Installation](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/installation.html)
- two FASTQ format read files per sample, obtained from a paired-end, short-read sequencing platform. 

## Pipeline outputs

The svAmplicon pipeline produces alignment bam files and an
app data package for interactive visualization, `*.mdi.package.zip`.

## Pipeline options

Options for the different pipeline actions can be listed as follows:

```sh
mdi svAmplicon <action> --help
```

## Pipeline job files

The best way to configure pipeline execution is using YAML job files.
The source publication provides comlete examples for the reported samples.
