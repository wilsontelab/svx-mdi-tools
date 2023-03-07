---
title: Stage 1 pipeline
parent: svCapture
has_children: false
nav_order: 10
published: true
---

## {{page.title}}

svCapture data analysis first entails execution of a pipeline with actions:
- **align** = align input fastq files to the reference genome
- **collate** = assemble read groups, make consensuses, and re-align to genome
- **extract** = scan name-sorted reads for anomalous molecules with alignment discontinuities
- **find** = scan anomalous molecules from one or more samples to make SV junction calls

This stepwise implementation:
- allows users to perform alignment outside of the svCapture pipeline
- supports simultaneous SV finding across multiple related samples
- permits code-sharing with other pipelines in the svx-mdi-tools code suite

As an example, a user might apply svCapture to multipe pre-existing
bam files and then combine the extracted data into a single 'find' operation.

## Pipeline inputs

In addition to the genome files described under
[Installation](https://wilsontelab.github.io/svx-mdi-tools/docs/installation/installation.html),
you need the following:

### Unaligned or pre-aligned reads

Most users will provide two FASTQ format read files obtained
from a paired-end, short-read sequencing platform. The path to these
read files are specified using options `--input-dir` and `--input-name`.

Alternatively, you can skip the `align` action and provide pre-aligned
reads as a name-sorted bam file using option `--bam-file`.

### Your capture target regions

svCapture requires that you provide a BED-format file via option `--targets-bed`
that lists all genome spans that were targeted - typically captured - 
for sequencing. These are used to categorize SV types, calculate coverage, etc.

### Unique molecular identifiers

If your library strategy included unique molecular identifiers, you
list and describe them using options `--umi-file` and `--umi-skip-bases`.
The first column in each row should be the sequence of one known UMI value
present at the end of a read sequence. As present, random UMIs are not supported.

## Parallel SV analysis across multiple samples

Actions `align`, `collate` and `extract` are executed per sample.
The `find` action can also be executed per sample, but often it is 
helpful to merge extracted anomalous reads so that SVs can be
discovered in a manner that combines information from more than one sample,
either to gain more evidence for the existence of an SV or to discover
that a candidate SV is not unique to one sample.

The two modes of SV finding are communicated by how you set options
`--output-dir` and `--data-name`. If the output directory already
contains extracted data files, a single-sample find is executed.
If the directory contains a set of sub-folders, each with extracted read files,
a merged-sample find is executed.

## Pipeline outputs

The svCapture pipeline performs extensive grouping and consensus
making to yield sets of output molecules that are deemed likely to correspond
to single, independent source DNA molecules.

The most important pipeline outputs are the lists of characterized SV junction calls:
- an R-compatible RDS file included in a date package for use in the svCapture app, *.find.structural_variants.rds
- a gzipped flat file, *.find.structural_variants.gz
- a VCF format file, *.find.structural_variants.vcf.bgz

Some additional output files are:
- the app data package for interactive visualization, *.mdi.package.zip

## Additional pipeline options

Other options for the different pipeline actions can be listed as follows:

```sh
mdi svCapture <action> --help
```
