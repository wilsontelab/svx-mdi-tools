---
title: svCapture
# parent: Stage 2 Apps
has_children: true
nav_order: 10
published: true
---

## svCapture

The svCapture pipeline and app offer a combination of 
functionalities - error-correction and SV junction identification -
optimized for detecting rare structural variants (SVs)
in short-read whole-genome sequencing libraries subjected to target capture, 
typically by probe hybridization.

The motivation for developing svCapture was to extend the principles
of error-corrected duplex sequencing to SVs, beyond just SNVs/indels.
Doing so requires careful attention to the distinct error modes
by which false SVs and SNVs arise, and simultaneous
bookkeeping of both source molecule DNA strands and (mis)alignments.
Once satisfied, confident detection of even single-molecule
nonhomologous SVs can then be asserted.

### Citation

svCapture and its use for SV identification is described in detail here:
- <bioRxiv>

### Usage overview

Data analysis first entails execution of a pipeline with steps:
- **align** = align input fastq files to the reference genome
- **collate** = assemble read groups, make consensuses, and re-align to genome
- **extract** = scan name-sorted reads for molecules with alignment discontinuities
- **find** = scan anomalous molecules from one or more samples to make SV junction calls

This stepwise implementation allows users to perform alignment and re-alignment
steps separately from the svCapture pipeline, supports analysis of SVs
in multiple related samples, and permits code-sharing
with other pipelines in the svx-mdi-tools code suite.

As an example, a user might apply svCapture to multipe pre-existing
bam files and then combine those into a single svCapture 'find' operation.






