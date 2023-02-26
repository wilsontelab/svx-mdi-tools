---
title: Stage 1 pipeline
parent: svCapture
has_children: false
nav_order: 10
published: true
---

## {{page.title}}

svCapture data analysis first entails execution of a pipeline with steps:
- **align** = align input fastq files to the reference genome
- **collate** = assemble read groups, make consensuses, and re-align to genome
- **extract** = scan name-sorted reads for molecules with alignment discontinuities
- **find** = scan anomalous molecules from one or more samples to make SV junction calls

This stepwise implementation:
- allows users to perform alignment steps outside of the svCapture pipeline
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
read files are specific using options `input-dir` and `input-name`.

Alternatively, you can skip the `align` step and provide pre-aligned
reads as a name-sorted bam file.

### Your capture target regions

svCapture requires that you provide a BED-format file that list
all of the individual genome spans that were targetet, typically captured,
for sequencing. These are used to categorize SV types, calculate coverage, etc.

### Unique molecular identifiers

If your library strategy included unique molecular identifiers, you
list and describe them using options `umi-file` and `umi-skip-bases`.

### Parallel SV analysis across multiple samples

Actions `align`, `collate` and `extract` are executed per sample.
The `find` action can also be execute per sample, but often it is 
helpful to merge extracted anomalous reads so that SVs can be
discovered in a manner that combines information from more than one sample.

The two modes of SV finding are communicated by how you set options
`output-dir` and `data-name`. If the target file output directory
contains extracted data files, a single-sample find is executed.
If the directory contains a set of sub-folders, each with extracted read files,
that a merged-sample find is executed.

## Pipeline outputs

The svCapture pipeline performs extensive grouping and consensus
making to yield sets of output molecules that are deemed likely to correspond
to single, independent source DNA molecules.


The most important pipeline outputs are the lists of characterized SV junction calls.

In addition

## Additional pipeline options

There are a number of other options for the different pipeline actions
that can be listed as follows:

```sh
mdi svCapture <action> --help
```
